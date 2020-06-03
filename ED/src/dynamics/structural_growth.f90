!==========================================================================================!
!==========================================================================================!
!    This is the module containing the routines that control growth of structural tissues, !
! and the change in population due to growth and mortality.                                !
!------------------------------------------------------------------------------------------!
module structural_growth
   contains

   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will control the structural growth of plants.                     !
   !                                                                                       !
   ! IMPORTANT: Do not change the order of operations below unless you know what you are   !
   !            doing.  Changing the order can affect the C/N budgets.                     !
   !---------------------------------------------------------------------------------------!
   subroutine dbstruct_dt(cgrid,veget_dyn_on,new_year)
      use ed_state_vars       , only : edtype                      & ! structure
                                     , polygontype                 & ! structure
                                     , sitetype                    & ! structure
                                     , patchtype                   ! ! structure
      use pft_coms            , only : seedling_mortality          & ! intent(in)
                                     , c2n_leaf                    & ! intent(in)
                                     , c2n_storage                 & ! intent(in)
                                     , c2n_recruit                 & ! intent(in)
                                     , c2n_stem                    & ! intent(in)
                                     , l2n_stem                    & ! intent(in)
                                     , negligible_nplant           & ! intent(in)
                                     , is_grass                    & ! intent(in)
                                     , agf_bs                      & ! intent(in)
                                     , q                           & ! intent(in)
                                     , storage_reflush_times       & ! intent(in)
                                     , is_liana                    & ! intent(in)
                                     , cbr_severe_stress           & ! intent(in)
                                     , h_edge                      & ! intent(in)
                                     , f_labile_leaf               & ! intent(in)
                                     , f_labile_stem               ! ! intent(in)
      use disturb_coms        , only : cl_fseeds_harvest           ! ! intent(in)
      use ed_max_dims         , only : n_pft                       & ! intent(in)
                                     , n_dbh                       ! ! intent(in)
      use ed_misc_coms        , only : igrass                      & ! intent(in)
                                     , ibigleaf                    & ! intent(in)
                                     , frqsumi                     & ! intent(in)
                                     , current_time                ! ! intent(in)
      use ed_therm_lib        , only : calc_veg_hcap               & ! function
                                     , update_veg_energy_cweh      ! ! function
      use physiology_coms     , only : ddmort_const                & ! intent(in)
                                     , iddmort_scheme              & ! intent(in)
                                     , cbr_scheme                  ! ! intent(in)
      use fuse_fiss_utils     , only : sort_cohorts                ! ! subroutine
      use stable_cohorts      , only : is_resolvable               ! ! subroutine
      use update_derived_utils, only : update_cohort_derived_props & ! subroutine
                                     , update_vital_rates          ! ! subroutine
      use allometry           , only : size2bl                     ! ! function
      use consts_coms         , only : yr_sec                      & ! intent(in)
                                     , almost_zero                 ! ! intent(in)
      use plant_hydro         , only : rwc2tw                      & ! subroutine
                                     , twi2twe                     ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      logical          , intent(in) :: veget_dyn_on
      logical          , intent(in) :: new_year
      !----- Local variables --------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      integer                       :: ipft
      integer                       :: prev_month
      integer                       :: prev_year
      integer                       :: prev_ndays
      integer                       :: imonth
      integer                       :: phenstatus_in
      integer                       :: krdepth_in
      real                          :: tor_fact
      real                          :: fgc_in_in
      real                          :: fsc_in_in
      real                          :: stgc_in_in
      real                          :: stsc_in_in
      real                          :: nplant_loss
      real                          :: pat_balive_in
      real                          :: pat_bdead_in
      real                          :: pat_bstorage_in
      real                          :: pat_mortality
      real                          :: pat_carbon_miss
      real                          :: carbon_miss
      real                          :: bleaf_in
      real                          :: broot_in
      real                          :: bsapwooda_in
      real                          :: bsapwoodb_in
      real                          :: bbarka_in
      real                          :: bbarkb_in
      real                          :: balive_in
      real                          :: bdeada_in
      real                          :: bdeadb_in
      real                          :: bevery_in
      real                          :: hite_in
      real                          :: dbh_in
      real                          :: nplant_in
      real                          :: bstorage_in
      real                          :: bstorage_reserve
      real                          :: agb_in
      real                          :: lai_in
      real                          :: wai_in
      real                          :: cai_in
      real                          :: ba_in
      real                          :: bag_in
      real                          :: bam_in
      real                          :: vm_bar_in
      real                          :: sla_in
      real                          :: psi_open_in
      real                          :: psi_closed_in
      real                          :: cb_act
      real                          :: cb_lightmax
      real                          :: cb_moistmax
      real                          :: cb_mlmax
      real                          :: cbr_light
      real                          :: cbr_moist
      real                          :: cbr_ml
      real                          :: f_bseeds
      real                          :: f_bdeada
      real                          :: f_bdeadb
      real                          :: f_growth
      real                          :: f_bstorage
      real                          :: a_bfast_mort_litter
      real                          :: b_bfast_mort_litter
      real                          :: a_bstruct_mort_litter
      real                          :: b_bstruct_mort_litter
      real                          :: a_bstorage_mort_litter
      real                          :: b_bstorage_mort_litter
      real                          :: a_bfast
      real                          :: b_bfast
      real                          :: a_bstruct
      real                          :: b_bstruct
      real                          :: a_bstorage
      real                          :: b_bstorage
      real                          :: maxh !< maximum patch height
      real                          :: mort_litter
      real                          :: bseeds_mort_litter
      real                          :: net_seed_N_uptake
      real                          :: net_stem_N_uptake
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      real                          :: old_leaf_water
      real                          :: old_wood_water
      real                          :: old_leaf_water_im2
      real                          :: old_wood_water_im2
      logical          , parameter  :: printout  = .false.
      character(len=17), parameter  :: fracfile  = 'struct_growth.txt'
      !----- Locally saved variables. -----------------------------------------------------!
      logical          , save       :: first_time = .true.
      !----- External function. -----------------------------------------------------------!
      integer          , external   :: num_days
      !------------------------------------------------------------------------------------!


      !----- First time, and the user wants to print the output.  Make a header. ----------!
      if (first_time) then

         !----- Make the header. ----------------------------------------------------------!
         if (printout) then
            open (unit=66,file=fracfile,status='replace',action='write')
            write (unit=66,fmt='(20(a,1x))')                                               &
                     '  YEAR',     '  MONTH',      '   DAY',      '   PFT',      '   ICO'  &
              ,'       BA_IN','      BAG_IN','      BAM_IN','      DBH_IN','   NPLANT_IN'  &
              ,'      BA_OUT','     BAG_OUT','     BAM_OUT','     DBH_OUT','  NPLANT_OUT'  &
              ,' TOTAL_BA_PY','TOTAL_BAG_PY','TOTAL_BAM_PY','TOTAL_BAR_PY','FIRST_CENSUS'
            close (unit=66,status='keep')
         end if
         !---------------------------------------------------------------------------------!

         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!


      !----- Previous month (for crop yield update and carbon balance mortality). ---------!
      prev_month = 1 + modulo(current_time%month-2,12)
      select case (prev_month)
      case (12)
         prev_year = current_time%year -1
      case default
         prev_year = current_time%year
      end select
      prev_ndays = num_days(prev_month,prev_year)
      !----- Correction factor for sapwood turnover rate. ---------------------------------!
      tor_fact   = real(prev_ndays) / yr_sec
      !------------------------------------------------------------------------------------!


      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)


         !----- Initialization. -----------------------------------------------------------!
         cgrid%crop_yield(prev_month,ipy) = 0.0
         !---------------------------------------------------------------------------------!



         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)
            !----- Initialization. --------------------------------------------------------!
            cpoly%basal_area       (:,:,isi) = 0.0
            cpoly%agb              (:,:,isi) = 0.0
            cpoly%crop_yield(prev_month,isi) = 0.0
            !------------------------------------------------------------------------------!

            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               !---------------------------------------------------------------------------!
               !    Here we must get the tallest tree cohort inside the patch. This will   !
               ! be then used to set the maximum height that lianas should aim at.         !
               ! Currently h_edge is set at 0.5m.  The arrested succession (only lianas    !
               ! inside the patch is a bit of a limit case: the maximum height will be set !
               ! to (0.5 + 0.5)m. It might still be reasonable.                            !
               !---------------------------------------------------------------------------!
               maxh = h_edge
               hcohortloop: do ico = 1,cpatch%ncohorts
                  if (.not. is_liana(cpatch%pft(ico)) .and. cpatch%hite(ico) > maxh) then
                     maxh = cpatch%hite(ico)
                  end if
               end do hcohortloop
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !    Add 0.5 to maxh has the effect of letting lianas grow 0.5 m above the  !
               ! tallest tree cohort. This number could be tweaked...                      !
               !---------------------------------------------------------------------------!
               maxh = maxh + h_edge
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Save patch-level litter inputs before growth balive.  We use these   !
               ! variables to check carbon conservation at the patch level.                !
               !---------------------------------------------------------------------------!
               fgc_in_in       = csite%fgc_in (ipa)
               fsc_in_in       = csite%fsc_in (ipa)
               stgc_in_in      = csite%stgc_in(ipa)
               stsc_in_in      = csite%stsc_in(ipa)
               pat_balive_in   = 0.0
               pat_bdead_in    = 0.0
               pat_bstorage_in = 0.0
               pat_mortality   = 0.0
               pat_carbon_miss = 0.0
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Loop over cohorts.                                                     !
               !---------------------------------------------------------------------------!
               cohortloop: do ico = 1,cpatch%ncohorts

                  !----- Assigning an alias for PFT type. ---------------------------------!
                  ipft    = cpatch%pft(ico)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Remember inputs in order to calculate increments (and to revert   !
                  ! back to these values later on in case vegetation dynamics is off).     !
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
                  bevery_in       = bleaf_in  + broot_in  + bsapwooda_in + bsapwoodb_in    &
                                  + bbarka_in + bbarkb_in + bdeada_in    + bdeadb_in
                  bag_in          = sum(cpoly%basal_area_growth(ipft,:,isi))
                  bam_in          = sum(cpoly%basal_area_mort(ipft,:,isi))
                  vm_bar_in       = cpatch%vm_bar          (ico)
                  sla_in          = cpatch%sla             (ico)
                  psi_open_in     = cpatch%psi_open        (ico)
                  psi_closed_in   = cpatch%psi_closed      (ico)
                  pat_balive_in   = pat_balive_in   + nplant_in * balive_in
                  pat_bdead_in    = pat_bdead_in    + nplant_in * (bdeada_in + bdeadb_in)
                  pat_bstorage_in = pat_bstorage_in + nplant_in * bstorage_in
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Save original heat capacitiy and water content for both leaves    !
                  ! and wood.  These are used to track changes in energy and water         !
                  ! storage due to vegetation dynamics.                                    !
                  !------------------------------------------------------------------------!
                  old_leaf_hcap      = cpatch%leaf_hcap     (ico)
                  old_wood_hcap      = cpatch%wood_hcap     (ico)
                  old_leaf_water     = cpatch%leaf_water    (ico)
                  old_wood_water     = cpatch%wood_water    (ico)
                  old_leaf_water_im2 = cpatch%leaf_water_im2(ico)
                  old_wood_water_im2 = cpatch%wood_water_im2(ico)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Apply mortality, and do not allow nplant < negligible_nplant (such  !
                  ! a sparse cohort is about to be terminated, anyway).                    !
                  ! NB: monthly_dlnndt will be zero or negative.                           !
                  !------------------------------------------------------------------------!
                  cpatch%monthly_dlnndt(ico) = max( cpatch%monthly_dlnndt (ico)            &
                                                  , log( negligible_nplant(ipft)           &
                                                       / cpatch%nplant    (ico) ) )
                  cpatch%nplant(ico)         = cpatch%nplant(ico)                          &
                                             * exp(cpatch%monthly_dlnndt(ico))
                  nplant_loss                = nplant_in - cpatch%nplant(ico)
                  !------------------------------------------------------------------------!



                  !----- Split biomass components that are labile or structural. ----------!
                  a_bfast    = f_labile_leaf(ipft) * cpatch%bleaf(ico)                     &
                             + f_labile_stem(ipft)                                         &
                             * ( cpatch%bsapwooda(ico) + cpatch%bbarka(ico)                &
                               + cpatch%bdeada   (ico) )
                  b_bfast    = f_labile_leaf(ipft) * cpatch%broot(ico)                     &
                             + f_labile_stem(ipft)                                         &
                             * ( cpatch%bsapwoodb(ico) + cpatch%bbarkb(ico)                &
                               + cpatch%bdeadb   (ico) )
                  a_bstruct  = (1.0 - f_labile_leaf(ipft)) * cpatch%bleaf(ico)             &
                             + (1.0 - f_labile_stem(ipft))                                 &
                             * ( cpatch%bsapwooda(ico) + cpatch%bbarka(ico)                &
                               + cpatch%bdeada   (ico) )
                  b_bstruct  = (1.0 - f_labile_leaf(ipft)) * cpatch%broot(ico)             &
                             + (1.0 - f_labile_stem(ipft))                                 &
                             * ( cpatch%bsapwoodb(ico) + cpatch%bbarkb(ico)                &
                               + cpatch%bdeadb   (ico) )
                  a_bstorage =        agf_bs(ipft)  * cpatch%bstorage(ico)
                  b_bstorage = (1.0 - agf_bs(ipft)) * cpatch%bstorage(ico)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Calculate litter owing to mortality.                               !
                  ! MLO - Use the actual change in nplant, which is derived from dlnndt.   !
                  !       The dlnndt is preferred because it linearises the demographic    !
                  !       ODE.  The litter inputs should be exactly the same as the gross  !
                  !       mortality loss (i.e. the actual change in nplant), otherwise     !
                  !       carbon is not conserved.                                         !
                  !------------------------------------------------------------------------!
                  a_bfast_mort_litter    = a_bfast    * nplant_loss
                  b_bfast_mort_litter    = b_bfast    * nplant_loss
                  a_bstruct_mort_litter  = a_bstruct  * nplant_loss
                  b_bstruct_mort_litter  = b_bstruct  * nplant_loss
                  a_bstorage_mort_litter = a_bstorage * nplant_loss
                  b_bstorage_mort_litter = b_bstorage * nplant_loss
                  mort_litter            = a_bfast_mort_litter    + b_bfast_mort_litter    &
                                         + a_bstruct_mort_litter  + b_bstruct_mort_litter  &
                                         + a_bstorage_mort_litter + b_bstorage_mort_litter
                  !------------------------------------------------------------------------!


                  !----- Integrate mortality losses. --------------------------------------!
                  pat_mortality = pat_mortality + mort_litter
                  !------------------------------------------------------------------------!


                  !----- Reset monthly mortality rate. ------------------------------------!
                  cpatch%monthly_dlnndt(ico) = 0.0
                  !------------------------------------------------------------------------!

                  !----- Calculate bstorage reserved for future refulushing needs ---------!
                  bstorage_reserve = (1.0 + q(ipft)) * storage_reflush_times(ipft)         &
                                   * size2bl(cpatch%dbh(ico),cpatch%hite(ico)              &
                                            ,cpatch%sla(ico),ipft)
                  !------------------------------------------------------------------------!


                  !----- Determine how to distribute what is in bstorage. -----------------!
                  call plant_structural_allocation(cpatch%pft(ico),cpatch%hite(ico)        &
                                                  ,cpatch%dbh(ico),cgrid%lat(ipy)          &
                                                  ,cpatch%phenology_status(ico)            &
                                                  ,cpatch%elongf(ico)                      &
                                                  ,bdeada_in,bdeadb_in, bstorage_in        &
                                                  ,bstorage_reserve, maxh                  &
                                                  ,f_bseeds,f_growth,f_bstorage)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Update bdead.   It turns out that when the fraction allocated to   !
                  ! sapwood is reasonable, structural growth causes the tree to have a     !
                  ! large debt for the following month, and GPP alone may not be           !
                  ! sufficient to make the trees with compatible balive.  This solution    !
                  ! attempts to correct for this, by assuming that f_growth is used to     !
                  ! grow all tissues, and only a fraction of f_growth goes to bdead.  The  !
                  ! remaining stays in bstorage and can be used to bring the cohort back   !
                  ! to allometry during the upcoming month (so respiration can be          !
                  ! properly accounted for).                                               !
                  !------------------------------------------------------------------------!
                  call bdead_structural_allocation(cpatch%pft(ico),cpatch%bstorage(ico)    &
                                                  ,bleaf_in,broot_in,bsapwooda_in          &
                                                  ,bsapwoodb_in,bbarka_in,bbarkb_in        &
                                                  ,bdeada_in,bdeadb_in,bevery_in           &
                                                  ,f_bstorage,f_growth,f_bdeada,f_bdeadb)
                  !------------------------------------------------------------------------!
                  cpatch%bdeada(ico) = cpatch%bdeada(ico) + f_bdeada * cpatch%bstorage(ico)
                  cpatch%bdeadb(ico) = cpatch%bdeadb(ico) + f_bdeadb * cpatch%bstorage(ico)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Assign NPP allocation.                                             !
                  !------------------------------------------------------------------------!
                  select case (ibigleaf)
                  case (0)
                     !------ NPP allocation to wood and coarse roots in KgC /m2 -----------!
                     cpatch%today_nppwood    (ico) = f_bdeada * cpatch%bstorage(ico)       &
                                                   * cpatch%nplant(ico)
                     cpatch%today_nppcroot   (ico) = f_bdeadb * cpatch%bstorage(ico)       &
                                                   * cpatch%nplant(ico)
                  end select
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Rebalance the plant nitrogen uptake considering the actual alloc- !
                  ! ation to structural growth.  This is necessary because c2n_stem does   !
                  ! not necessarily equal c2n_storage.                                     !
                  !------------------------------------------------------------------------!
                  net_stem_N_uptake = ( cpatch%bdeada(ico) + cpatch%bdeadb(ico)            &
                                      - bdeada_in          - bdeadb_in          )          &
                                    * cpatch%nplant(ico)                                   &
                                    * ( 1.0 / c2n_stem(cpatch%pft(ico)) - 1.0 / c2n_storage)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Calculate total seed production and seed litter.  The seed pool   !
                  ! gets a fraction f_bseeds of bstorage.  In case this is agriculture, we !
                  ! also set a fraction of seeds as crop yield.                            !
                  !------------------------------------------------------------------------!
                  select case (csite%dist_type(ipa))
                  case (8)
                     !---------------------------------------------------------------------!
                     !      Cropland patch.  Here we must account for harvesting of seeds. !
                     ! Leaves and non-structural carbon must be harvested at the end of    !
                     ! the cropland cycle.                                                 !
                     !---------------------------------------------------------------------!


                     !----- bseeds contains only non-harvested seeds. ---------------------!
                     cpatch%bseeds(ico) = (1.-cl_fseeds_harvest)                           &
                                        * f_bseeds * cpatch%bstorage(ico)
                     cpatch%byield(ico) = cl_fseeds_harvest                                &
                                        * f_bseeds * cpatch%bstorage(ico)
                     !---------------------------------------------------------------------!


                  case default
                     !---------------------------------------------------------------------!
                     !      Not a cropland patch. No harvesting should occur.              !
                     !---------------------------------------------------------------------!
                     cpatch%bseeds(ico) = f_bseeds * cpatch%bstorage(ico)
                     cpatch%byield(ico) = 0.
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !      Net primary productivity used for seed production.  This must     !
                  ! include seeds that have been harvested.                                !
                  !------------------------------------------------------------------------!
                  cpatch%today_NPPseeds(ico) = cpatch%nplant(ico)                          &
                                             * f_bseeds * cpatch%bstorage(ico)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Send dead seed and failed recruits to the litter pool (this time, !
                  ! exclude fraction of bseeds that was harvested.                         !
                  !------------------------------------------------------------------------!
                  bseeds_mort_litter = cpatch%bseeds(ico) * cpatch%nplant(ico)             &
                                     * seedling_mortality(ipft)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Split seedling mortality.  This pools is comprised of several     !
                  ! things: flowers, fruits, seeds, and failed recruits.  Because we do    !
                  ! not have the actual fractions that go to each pool, we 
                  !------------------------------------------------------------------------!
                  !------------------------------------------------------------------------!



                  !----- Integrate mortality losses. --------------------------------------!
                  pat_mortality = pat_mortality + bseeds_mort_litter
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Rebalance the plant nitrogen uptake considering the actual alloc- !
                  ! ation to seeds.  This is necessary because c2n_recruit does not have   !
                  ! to be equal to c2n_storage.  Here we must also account for harvested   !
                  !seeds.                                                                  !
                  !------------------------------------------------------------------------!
                  net_seed_N_uptake =  cpatch%nplant(ico)                                  &
                                    * f_bseeds * cpatch%bstorage(ico)                      &
                                    * (1.0 / c2n_recruit(ipft) - 1.0 / c2n_storage)
                  !------------------------------------------------------------------------!


                  !----- Decrement the storage pool due to allocation. --------------------!
                  cpatch%bstorage(ico) = cpatch%bstorage(ico) * f_bstorage
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Finalize litter inputs.                                            !
                  !------------------------------------------------------------------------!
                  csite%fgc_in (ipa) = csite%fgc_in(ipa) + a_bfast_mort_litter             &
                                     + a_bstorage_mort_litter + bseeds_mort_litter
                  csite%fsc_in (ipa) = csite%fsc_in(ipa) + b_bfast_mort_litter             &
                                     + b_bstorage_mort_litter
                  csite%fgn_in (ipa) = csite%fgn_in(ipa)                                   &
                                     + a_bfast_mort_litter    / c2n_leaf   (ipft)          &
                                     + a_bstorage_mort_litter / c2n_storage                &
                                     + bseeds_mort_litter     / c2n_recruit(ipft)
                  csite%fsn_in (ipa) = csite%fsn_in(ipa)                                   &
                                     + b_bfast_mort_litter    / c2n_leaf   (ipft)          &
                                     + b_bstorage_mort_litter / c2n_storage
                  csite%stgc_in(ipa) = csite%stgc_in(ipa) + a_bstruct_mort_litter
                  csite%stsc_in(ipa) = csite%stsc_in(ipa) + b_bstruct_mort_litter
                  csite%stgl_in(ipa) = csite%stgl_in(ipa)                                  &
                                     + a_bstruct_mort_litter * l2n_stem / c2n_stem(ipft)
                  csite%stsl_in(ipa) = csite%stsl_in(ipa)                                  &
                                     + b_bstruct_mort_litter * l2n_stem / c2n_stem(ipft)
                  csite%stgn_in(ipa) = csite%stgn_in(ipa)                                  &
                                     + a_bstruct_mort_litter  / c2n_stem   (ipft)
                  csite%stsn_in(ipa) = csite%stsn_in(ipa)                                  &
                                     + b_bstruct_mort_litter  / c2n_stem   (ipft)
                  csite%total_plant_nitrogen_uptake(ipa) =                                 &
                         csite%total_plant_nitrogen_uptake(ipa) + net_seed_N_uptake        &
                       + net_stem_N_uptake
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Yield will leave the patch, accumulate it in the loss term for     !
                  ! carbon budget, and the total site yield productivity.                  !
                  !------------------------------------------------------------------------!
                  csite%cbudget_loss2yield(ipa)    = csite%cbudget_loss2yield(ipa)         &
                                                   + cpatch%nplant(ico)                    &
                                                   * cpatch%byield(ico) * frqsumi
                  cpoly%crop_yield(prev_month,isi) = cpoly%crop_yield(prev_month,isi)      &
                                                   + cpatch%nplant(ico)                    &
                                                   * cpatch%byield(ico) * csite%area(ipa)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Calculate some derived cohort properties:                          !
                  ! - DBH                                                                  !
                  ! - Height                                                               !
                  ! - Recruit and census status                                            !
                  ! - Phenology status                                                     !
                  ! - Area indices                                                         !
                  ! - Basal area                                                           !
                  ! - AGB                                                                  !
                  ! - Rooting depth                                                        !
                  !------------------------------------------------------------------------!
                  call update_cohort_derived_props(cpatch,ico,cpoly%lsl(isi),new_year      &
                                                  ,cpoly%llspan_toc(ipft,isi)              &
                                                  ,cpoly%vm_bar_toc(ipft,isi)              &
                                                  ,cpoly%rd_bar_toc(ipft,isi)              &
                                                  ,cpoly%sla_toc   (ipft,isi) )
                  !------------------------------------------------------------------------!


                  !----- Update annual average carbon balances for mortality. -------------!
                  cpatch%cb          (prev_month,ico) = cpatch%cb          (13,ico)
                  cpatch%cb_lightmax (prev_month,ico) = cpatch%cb_lightmax (13,ico)
                  cpatch%cb_moistmax (prev_month,ico) = cpatch%cb_moistmax (13,ico)
                  cpatch%cb_mlmax    (prev_month,ico) = cpatch%cb_mlmax    (13,ico)
                  !------------------------------------------------------------------------!

                  !----- Update monhtly average PLC for hydraulic  mortality. -------------!
                  cpatch%plc_monthly(prev_month,ico) = cpatch%plc_monthly(13,ico)
                  cpatch%plc_monthly(13,ico) = 0.0
                  !------------------------------------------------------------------------!

                  !----- If monthly files are written, save the current carbon balance. ---!
                  if (associated(cpatch%mmean_cb)) then
                     cpatch%mmean_cb(ico)         = cpatch%cb(13,ico)
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Reset the current month integrator.  The initial value may depend !
                  ! on the storage.  By including this term we make sure that plants won't !
                  ! start dying as soon as they shed their leaves, but only when they are  !
                  ! in negative carbon balance and without storage.  This is done only     !
                  ! when iddmort_scheme is set to 1, otherwise the initial value is 0.     !
                  !------------------------------------------------------------------------!
                  select case (iddmort_scheme)
                  case (0)
                     !------ Storage is not accounted. ------------------------------------!
                     cpatch%cb          (13,ico) = 0.0
                     cpatch%cb_lightmax (13,ico) = 0.0
                     cpatch%cb_moistmax (13,ico) = 0.0
                     cpatch%cb_mlmax    (13,ico) = 0.0

                  case (1)
                   !------ Storage is accounted. ------------------------------------------!
                     cpatch%cb          (13,ico) = cpatch%bstorage(ico)
                     cpatch%cb_lightmax (13,ico) = cpatch%bstorage(ico)
                     cpatch%cb_moistmax (13,ico) = cpatch%bstorage(ico)
                     cpatch%cb_mlmax    (13,ico) = cpatch%bstorage(ico)
                  end select
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !  Set up CB/CBmax as running sums and use that in the calculate cbr     !
                  !------------------------------------------------------------------------!

                  !----- Initialize with 0 ------------------------------------------------!
                  cb_act       = 0.0
                  cb_lightmax  = 0.0
                  cb_moistmax  = 0.0
                  cb_mlmax     = 0.0

                  !------------------------------------------------------------------------!
                  !      Compute the relative carbon balance.                              !
                  !------------------------------------------------------------------------!
                  if (is_grass(ipft).and. igrass==1) then  
                     !----- Grass loop, use past month's CB only. -------------------------!
                     cb_act      =  cpatch%cb          (prev_month,ico)
                     cb_lightmax =  cpatch%cb_lightmax (prev_month,ico)
                     cb_moistmax =  cpatch%cb_moistmax (prev_month,ico)
                     cb_mlmax    =  cpatch%cb_mlmax    (prev_month,ico)
                     !---------------------------------------------------------------------!
                  else  
                     !----- Tree loop, use annual average carbon balance. -----------------!
                     do imonth = 1,12
                        cb_act      = cb_act      + cpatch%cb          (imonth,ico)
                        cb_lightmax = cb_lightmax + cpatch%cb_lightmax (imonth,ico)
                        cb_moistmax = cb_moistmax + cpatch%cb_moistmax (imonth,ico)
                        cb_mlmax    = cb_mlmax    + cpatch%cb_mlmax    (imonth,ico)
                     end do
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!


                  !----- Light-related carbon balance. ------------------------------------!
                  if (cb_lightmax > 0.0) then
                     cbr_light = min(1.0, cb_act / cb_lightmax)
                  else
                     cbr_light = cbr_severe_stress(ipft)
                  end if
                  !------------------------------------------------------------------------!


                  !----- Soil moisture-related carbon balance. ----------------------------!
                  if (cb_moistmax > 0.0) then
                     cbr_moist = min(1.0, cb_act / cb_moistmax )
                  else
                     cbr_moist = cbr_severe_stress(ipft)
                  end if
                  !------------------------------------------------------------------------!


                  !----- Soil moisture+light related carbon balance. ----------------------!
                  if (cb_mlmax > 0.0) then
                     cbr_ml    = min(1.0, cb_act / cb_mlmax )
                  else
                     cbr_ml    = cbr_severe_stress(ipft)
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !  calculate CBR according to the specified CBR_SCHEME                   !
                  !------------------------------------------------------------------------!
                  select case (cbr_scheme)
                  case (0)
                     !----- CBR from absolute max CB --------------------------------------!
                     cpatch%cbr_bar(ico) = max(cbr_ml, cbr_severe_stress(ipft))
                     !---------------------------------------------------------------------!


                  case (1)
                     !---------------------------------------------------------------------!
                     !   CBR from combination of light & moist CBR                         !
                     !   Relative carbon balance: a combination of the two factors.        !
                     !---------------------------------------------------------------------!
                     if ( cbr_light <= cbr_severe_stress(ipft) .and.                       &
                          cbr_moist <= cbr_severe_stress(ipft)       ) then
                        cpatch%cbr_bar(ico) = cbr_severe_stress(ipft)
                     else
                        cpatch%cbr_bar(ico) = cbr_severe_stress(ipft)                      &
                                            + ( cbr_light - cbr_severe_stress(ipft) )      &
                                            * ( cbr_moist - cbr_severe_stress(ipft) )      &
                                            / (        ddmort_const  * cbr_moist           &
                                              + (1.0 - ddmort_const) * cbr_light           &
                                              - cbr_severe_stress(ipft) )
                     end if
                     !---------------------------------------------------------------------!

                  case (2)
                     !----- CBR from most limiting CBR ------------------------------------!
                     cpatch%cbr_bar(ico) = max( min(cbr_moist, cbr_light)                  &
                                              , cbr_severe_stress(ipft) )
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!



                  !----- Update interesting output quantities. ----------------------------!
                  call update_vital_rates(cpatch,ico,dbh_in,nplant_in,agb_in,ba_in         &
                                         ,csite%area(ipa),cpoly%basal_area(:,:,isi)        &
                                         ,cpoly%agb(:,:,isi)                               &
                                         ,cpoly%basal_area_growth(:,:,isi)                 &
                                         ,cpoly%agb_growth(:,:,isi)                        &
                                         ,cpoly%basal_area_mort(:,:,isi)                   &
                                         ,cpoly%agb_mort(:,:,isi))
                  !------------------------------------------------------------------------!


                  !----- Record monthly diameter growth  ----------------------------------!
                  cpatch%ddbh_monthly(prev_month,ico) = cpatch%ddbh_dt(ico) !cm/yr
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     In case vegetation dynamics is off, revert back to previous        !
                  ! values:                                                                !
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
                     cpatch%vm_bar          (ico) = vm_bar_in    
                     cpatch%sla             (ico) = sla_in       
                     cpatch%psi_open        (ico) = psi_open_in  
                     cpatch%psi_closed      (ico) = psi_closed_in
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  ! MLO. We now update the heat capacity and the vegetation internal       !
                  !      energy.  Since no energy or water balance is done here, we simply !
                  !      update the energy in order to keep the same temperature and water !
                  !      as before.  Internal energy is an extensive variable, we just     !
                  !      account for the difference in the heat capacity to update it.     !
                  !------------------------------------------------------------------------!
                  call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdeada(ico)                  &
                                    ,cpatch%bsapwooda(ico),cpatch%bbarka(ico)              &
                                    ,cpatch%nplant(ico),cpatch%pft(ico)                    &
                                    ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
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
                  call is_resolvable(csite,ipa,ico,.false.,.false.,'dbstruct_dt')
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  if (printout) then
                     open (unit=66,file=fracfile,status='old',position='append'            &
                          ,action='write')
                     write (unit=66,fmt='(5(i6,1x),14(f12.6,1x),1(11x,l1,1x))')            &
                        current_time%year,current_time%month,current_time%date,ipft,ico    &
                        ,ba_in,bag_in,bam_in,dbh_in,nplant_in                              &
                        ,cpatch%basarea(ico)                                               &
                        ,sum(cpoly%basal_area_growth(ipft,:,isi))                          &
                        ,sum(cpoly%basal_area_mort(ipft,:,isi))                            &
                        ,cpatch%dbh(ico),cpatch%nplant(ico)                                &
                        ,cgrid%total_basal_area(ipy)                                       &
                        ,cgrid%total_basal_area_growth(ipy)                                &
                        ,cgrid%total_basal_area_mort(ipy)                                  &
                        ,cgrid%total_basal_area_recruit(ipy)                               &
                        ,cpatch%first_census(ico)
                     close (unit=66,status='keep')
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Before we let these plants carry on with their lives, we must check !
                  ! that we can account for all changes in carbon.  In case there is any-  !
                  ! thing missing or excessive, stop the simulation.                       !
                  !------------------------------------------------------------------------!
                  if (veget_dyn_on) then
                     call check_bstruct_cohort(csite,ipa,ico,bleaf_in,broot_in             &
                                              ,bsapwooda_in,bsapwoodb_in,bbarka_in         &
                                              ,bbarkb_in,balive_in,bstorage_in,bdeada_in   &
                                              ,bdeadb_in,f_bseeds,f_growth,f_bdeada        &
                                              ,f_bdeadb,f_bstorage,carbon_miss)
                     pat_carbon_miss = pat_carbon_miss + carbon_miss
                  end if
                  !------------------------------------------------------------------------!


               end do cohortloop
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Make sure that the patch did not try to smuggle or evade carbon.      !
               ! For this we compare the live stocks and the necromass inputs before and   !
               ! after updating structural and reproduction tissues.                       !
               !---------------------------------------------------------------------------!
               if (veget_dyn_on) then
                  call check_bstruct_patch(csite,ipa,fgc_in_in,fsc_in_in,stgc_in_in        &
                                          ,stsc_in_in,pat_balive_in,pat_bdead_in           &
                                          ,pat_bstorage_in,pat_carbon_miss,pat_mortality)
               end if
               !---------------------------------------------------------------------------!



            end do patchloop
            !------------------------------------------------------------------------------!


            !------ Update polygon-level crop yield. --------------------------------------!
            cgrid%crop_yield(prev_month,ipy) = cgrid%crop_yield(prev_month,ipy)            &
                                             + cpoly%crop_yield(prev_month,isi)            &
                                             * cpoly%area(isi)
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!
      return
   end subroutine dbstruct_dt
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will decide the partition of storage biomass into seeds and dead  !
   ! (structural) biomass.                                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine plant_structural_allocation(ipft,hite,dbh,lat,phen_status,elongf,bdeada      &
                                         ,bdeadb,bstorage,bstorage_reserve,maxh            &
                                         ,f_bseeds,f_growth,f_bstorage)
      use pft_coms      , only : phenology      & ! intent(in)
                               , repro_min_h    & ! intent(in)
                               , repro_min_dbh  & ! intent(in)
                               , hgt_max        & ! intent(in)
                               , r_bang         & ! intent(in)
                               , r_fract        & ! intent(in)
                               , r_cv50         & ! intent(in)
                               , st_fract       & ! intent(in)
                               , dbh_crit       & ! intent(in)
                               , is_grass       & ! intent(in)
                               , is_liana       ! ! intent(in)
      use ed_misc_coms  , only : current_time   & ! intent(in)
                               , igrass         & ! intent(in)
                               , ibigleaf       ! ! intent(in)
      use consts_coms   , only : r_tol_trunc    & ! intent(in)
                               , tiny_num       ! ! intent(in)
      use allometry     , only : size2bd        & ! intent(in)
                               , h2dbh          ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer, intent(in)  :: ipft
      real   , intent(in)  :: hite
      real   , intent(in)  :: dbh
      real   , intent(in)  :: lat
      real   , intent(in)  :: bdeada     !> Current dead biomass
      real   , intent(in)  :: bdeadb     !> Current dead biomass
      real   , intent(in)  :: bstorage   !> Current storage pool
      real   , intent(in)  :: bstorage_reserve !> Target bstorage reserve for reflushing
      real   , intent(in)  :: maxh       !> Height of the tallest cohort in the patch
      integer, intent(in)  :: phen_status
      real   , intent(in)  :: elongf     !> Elongation factor
      real   , intent(out) :: f_bseeds   !> Fraction to use for reproduction
      real   , intent(out) :: f_growth   !> Fraction to use for growth
      real   , intent(out) :: f_bstorage !> Fraction to keep as storage
      !----- Local variables --------------------------------------------------------------!
      real                         :: bd_target   !> Target Bd to reach maxh height
      real                         :: delta_bd    !> Target Bd - actual Bd
      real                         :: dnorm       !> Normalised DBH
      real                         :: r_fract_act !> Hgt-dependent reproduction allocation
      logical                      :: late_spring
      logical                      :: use_storage
      logical                      :: zero_growth
      logical                      :: zero_repro
      logical          , parameter :: printout  = .false.
      character(len=13), parameter :: fracfile  = 'storalloc.txt'
      !----- Locally saved variables. -----------------------------------------------------!
      logical          , save      :: first_time = .true.
      !------------------------------------------------------------------------------------!


      !----- First time, and the user wants to print the output.  Make a header. ----------!
      if (first_time) then

         !----- Make the header. ----------------------------------------------------------!
         if (printout) then
            open (unit=66,file=fracfile,status='replace',action='write')
            write (unit=66,fmt='(18(a,1x))')                                               &
               '  YEAR',' MONTH','   DAY','   PFT','PHENOL','PH_STT','SPRING',' GRASS'     &
              ,'      HEIGHT','   REPRO_HGT','         DBH','    DBH_CRIT','      ELONGF'  &
              ,'       BDEAD','    BSTORAGE','   F_STORAGE','     F_SEEDS','    F_GROWTH'
            close (unit=66,status='keep')
         end if
         !---------------------------------------------------------------------------------!

         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!

      
      !------------------------------------------------------------------------------------!
      ! Check whether plants want to reserve bstorage for reflushing leaves and roots      !
      !------------------------------------------------------------------------------------!
      if (bstorage <= bstorage_reserve) then
          ! plants want to save bstorage for potential future needs
          f_bseeds = 0.0
          f_growth = 0.0
          f_bstorage = 1.0 - f_growth - f_bseeds
          return
      endif
      !------------------------------------------------------------------------------------!



      !----- Check whether this is late spring... -----------------------------------------!
      late_spring = (lat >= 0.0 .and. current_time%month ==  6) .or.                       &
                    (lat <  0.0 .and. current_time%month == 12)
      !------------------------------------------------------------------------------------!

      !----- Use storage. -----------------------------------------------------------------!
      use_storage = (phenology(ipft) /= 2   .or.  late_spring) .and.                       &
                    phen_status == 0  .and. bstorage > 0.0
      !------------------------------------------------------------------------------------!



      !----- Find the current target for allocation to reproduction. ----------------------!
      if (r_bang(ipft)) then
         !----- "Bang" reproduction once plant reaches reproductive maturity. -------------!
         if ( hite <  ( (1.0-r_tol_trunc) * repro_min_h(ipft) ) ) then
            r_fract_act = 0.0
         else
            r_fract_act = min(r_fract(ipft), 1.0 - st_fract(ipft))
         end if
         !---------------------------------------------------------------------------------!
      else
         !----- Find normalised height to calculate asymptote. ----------------------------!
         if ( repro_min_dbh(ipft) < ( (1.0 -r_tol_trunc) * dbh_crit(ipft) ) ) then
            dnorm = (dbh - repro_min_dbh(ipft)) / (dbh_crit(ipft) - repro_min_dbh(ipft))
            dnorm = max(0.,dnorm)
         elseif ( dbh < ( (1.0 -r_tol_trunc) * repro_min_dbh(ipft) ) ) then
            dnorm = 0.
         else
            dnorm = 1.
         end if
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Find allocation to reproduction.  In case the denominator is zero, assume   !
         ! partial bang (sensu Wenk and Falster 2015).                                     !
         !---------------------------------------------------------------------------------!
         if ( (dnorm + r_cv50(ipft)) >= tiny_num ) then
            r_fract_act = (1.0 - st_fract(ipft)) * dnorm / (dnorm + r_cv50(ipft))
         elseif (dnorm >= tiny_num) then
            r_fract_act = (1.0 - st_fract(ipft))
         else
            r_fract_act = 0.0
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Check size structure to decide how to allocate carbon to growth, reproduction  !
      ! and storage.                                                                       !
      !------------------------------------------------------------------------------------!
      select case (ibigleaf)
      case (0)
         !---------------------------------------------------------------------------------!
         !      Size and age structure.  Calculate fraction of bstorage going to bdead and !
         ! reproduction.  First we must make sure that the plant should do something here. !
         ! A plant should not allocate anything to reproduction or growth if it is not the !
         ! right time of year (for cold deciduous plants), or if the plants are actively   !
         ! dropping leaves or off allometry.                                               !
         !---------------------------------------------------------------------------------!
         if (use_storage) then
            !------------------------------------------------------------------------------!
            !      Decide allocation to seeds and heartwood based on size and life form.   !
            !------------------------------------------------------------------------------!
            if (is_liana(ipft)) then
               zero_growth = hite >= ( (1.0-r_tol_trunc) * maxh )
               zero_repro  = hite <  ( (1.0-r_tol_trunc) * repro_min_h(ipft) )

               !---------------------------------------------------------------------------!
               !    Lianas: we must check height relative to the rest of the local plant   !
               ! community.                                                                !
               !---------------------------------------------------------------------------!
               if (zero_growth) then
                  f_bseeds = merge(0.0, r_fract_act, zero_repro)
                  f_growth = 0.0
               else
                  bd_target = size2bd(h2dbh(maxh,ipft),maxh,ipft)
                  delta_bd  = bd_target - bdeada - bdeadb
                  !------------------------------------------------------------------------!
                  !    If bstorage is 0 or lianas have already reached their bd_target     !
                  ! don't grow otherwise invest what is needed (or everything in case it's !
                  ! not enough) to reach bd_target.                                        !
                  !    For seeds first check if the cohort has reached repro_min_h. In     !
                  ! case so, then check that it has enough left to invest r_fract_act in   !
                  ! reproduction.  In case so, invest r_fract_act in reproduction and the  !
                  ! rest will stay in storage.  Otherwise, invest everything in            !
                  ! reproduction.  Finally in case the liana hasn't reached repro_min_h,   !
                  ! leave everything in storage.                                           !
                  !------------------------------------------------------------------------!
                  f_growth  = merge(0.0                                                    &
                                   ,min(delta_bd / bstorage, 1.0)                          &
                                   ,bstorage * delta_bd <= 0.0)
                  f_bseeds  = merge( 0.0, min(r_fract_act,1.0-f_growth),zero_repro)
               end if
               !---------------------------------------------------------------------------!
            else if (is_grass(ipft) .and. igrass == 1) then
               !---------------------------------------------------------------------------!
               !     New grasses don't growth here (they do in dbalive_dt).  Decide        !
               ! whether they may reproduce or not.                                        !
               !---------------------------------------------------------------------------!
               zero_repro = hite <  ( (1.0-r_tol_trunc) * repro_min_h(ipft) )
               if (zero_repro) then
                  f_bseeds = 0.0
               else
                  f_bseeds = 1.0 - st_fract(ipft)
               end if
               f_growth = 0.0
               !---------------------------------------------------------------------------!
            else
               !---------------------------------------------------------------------------!
               !    Trees and old grasses.  Currently the only difference is that grasses  !
               ! stop growing once they reach maximum height (as they don't have an actual !
               ! DBH).                                                                     !
               !---------------------------------------------------------------------------!
               zero_growth = is_grass(ipft) .and.                                          &
                             hite >= ( (1.0-r_tol_trunc) * hgt_max(ipft)     )
               !---------------------------------------------------------------------------!


               !----- Decide allocation based on size. ------------------------------------!
               if (zero_growth) then
                  f_bseeds = 1.0 - st_fract(ipft)
                  f_growth = 0.0
               else
                  f_bseeds = r_fract_act
                  f_growth = max(0.0,1.0 - st_fract(ipft) - r_fract_act)
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         else  
            !----- Plant should not allocate carbon to seeds or grow new biomass. ---------!
            f_growth = 0.0
            f_bseeds = 0.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      case (1)
         !---------------------------------------------------------------------------------!
         !    Big-leaf solver.  As long as it is OK to grow, everything goes into          !
         ! 'reproduction'.  This will ultimately be used to increase NPLANT of the         !
         ! 'big leaf' cohort.                                                              !
         !---------------------------------------------------------------------------------!
         if (use_storage) then
            !------------------------------------------------------------------------------!
            ! A plant should only grow if it is the right time of year (for cold deciduous !
            ! plants), or if the plants are not actively dropping leaves or off allometry. !
            !------------------------------------------------------------------------------!
            f_bseeds = 1.0 - st_fract(ipft)
            f_growth = 0.0
         else
            f_growth = 0.0
            f_bseeds = 0.0
         end if
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Make sure allocation to growth and seeds is not negligible.                    !
      !------------------------------------------------------------------------------------!
      if (f_bseeds < r_tol_trunc) f_bseeds = 0.0
      if (f_growth < r_tol_trunc) f_growth = 0.0

      ! plants will save bstorage_reserve as mobile carbon supply for potential reflushing
      ! we take the maximum of the residual of f_bseeds and f_growth and bstorage_reserve 
      ! / bstorage to maintain compatibility with previous versions (bstorage_reserve = 0.)
      f_bstorage = max(1.0 - f_bseeds - f_growth, bstorage_reserve / bstorage)

      ! we need to modify f_bseeds and f_growth accordingly if f_bseeds+f_growth is non-zero
      if (f_bseeds + f_growth > 0.0) then
          f_bseeds = f_bseeds / (f_bseeds + f_growth) * (1.0 - f_bstorage)
          f_growth = f_growth / (f_bseeds + f_growth) * (1.0 - f_bstorage)
      endif
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      if (printout) then
         open (unit=66,file=fracfile,status='old',position='append',action='write')
         write (unit=66,fmt='(6(i6,1x),2(5x,l1,1x),10(f12.6,1x))')                         &
               current_time%year,current_time%month,current_time%date,ipft,phenology(ipft) &
              ,phen_status,late_spring,is_grass(ipft),hite,repro_min_h(ipft),dbh           &
              ,dbh_crit(ipft),elongf,bdeada+bdeadb,bstorage,f_bstorage,f_bseeds,f_growth
         close (unit=66,status='keep')
      end if
      !------------------------------------------------------------------------------------!
      return
   end subroutine plant_structural_allocation
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine calculates the fraction of carbon that should be allocated to    !
   ! heartwood (bdeada and bdeadb), based on the growth allocation and the optimal size of !
   ! the enlarged plant.  When the fraction allocated to sapwood is reasonable, structural !
   ! growth causes the tree to have a large debt for the following month to go back to     !
   ! allometry, and GPP alone may not be sufficient to make the trees with compatible      !
   ! balive within one month.  When new allometry is used, f_growth becomes the fraction   !
   ! of allocation to grow all tissues, and only a fraction of f_growth goes to bdead.     !
   ! The remaining stays in bstorage and can be used to bring the cohort back to allometry !
   ! during the upcoming month (so respiration can be properly accounted for).             !
   !---------------------------------------------------------------------------------------!
   subroutine bdead_structural_allocation(ipft,bstorage_in,bleaf_in,broot_in,bsapwooda_in  &
                                         ,bsapwoodb_in,bbarka_in,bbarkb_in,bdeada_in       &
                                         ,bdeadb_in,bevery_in,f_bstorage,f_growth,f_bdeada &
                                         ,f_bdeadb)
      use physiology_coms, only : istruct_growth_scheme  ! ! intent(in)
      use allometry      , only : expand_bevery          ! ! subroutine
      use consts_coms    , only : almost_zero            & ! intent(in)
                                , r_tol_trunc            ! ! intent(in)
      use pft_coms       , only : agf_bs                 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer, intent(in)    :: ipft
      real   , intent(in)    :: bstorage_in
      real   , intent(in)    :: bleaf_in
      real   , intent(in)    :: broot_in
      real   , intent(in)    :: bsapwooda_in
      real   , intent(in)    :: bsapwoodb_in
      real   , intent(in)    :: bbarka_in
      real   , intent(in)    :: bbarkb_in
      real   , intent(in)    :: bdeada_in
      real   , intent(in)    :: bdeadb_in
      real   , intent(in)    :: bevery_in
      real   , intent(inout) :: f_bstorage
      real   , intent(inout) :: f_growth
      real   , intent(out)   :: f_bdeada
      real   , intent(out)   :: f_bdeadb
      !----- Local variables. -------------------------------------------------------------!
      real                   :: dbh_aim
      real                   :: hite_aim
      real                   :: bleaf_aim
      real                   :: broot_aim
      real                   :: bsapwooda_aim
      real                   :: bsapwoodb_aim
      real                   :: bbarka_aim
      real                   :: bbarkb_aim
      real                   :: balive_aim
      real                   :: bdeada_aim
      real                   :: bdeadb_aim
      real                   :: bevery_aim
      real                   :: tr_bleaf
      real                   :: tr_broot
      real                   :: tr_bsapwooda
      real                   :: tr_bsapwoodb
      real                   :: tr_bbarka
      real                   :: tr_bbarkb
      real                   :: tr_bdeada
      real                   :: tr_bdeadb
      real                   :: tr_every
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Plant is not supposed to go to growth, or there is so little storage that it   !
      ! is not worth growing.  Skip growing and instead wait until storage accumulates.    !
      !------------------------------------------------------------------------------------!
      if (f_growth <= almost_zero .or. (bstorage_in < (r_tol_trunc * bevery_in)) ) then
         f_bstorage = f_bstorage + f_growth
         f_growth   = 0.
         f_bdeada   = 0.
         f_bdeadb   = 0.
         return
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Allocation will depend on the structural growth option.  This option decides  !
      ! whether to allocate all storage available to growth to heartwood (option 0) or to  !
      ! growth of heartwood and living tissues (option 1).                                 !
      !------------------------------------------------------------------------------------!
      select case (istruct_growth_scheme)
      case (1)
         !----- Find the new biomass with the storage inputs. -----------------------------!
         bevery_aim = bevery_in + f_growth * bstorage_in
         call expand_bevery(ipft,bevery_aim,dbh_aim,hite_aim,bleaf_aim,broot_aim           &
                           ,bsapwooda_aim,bsapwoodb_aim,bbarka_aim,bbarkb_aim,balive_aim   &
                           ,bdeada_aim,bdeadb_aim)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Find the potential transfer rates.  It is possible that some of the        !
         ! tissues are already at the ideal or even exceeding the target biomass (e.g. in  !
         ! case of extremely low storage availability).  In this case, we only account and !
         ! scale allocation with positive increments.                                      !
         !---------------------------------------------------------------------------------!
         tr_bleaf     = max( 0.0, bleaf_aim     - bleaf_in     )
         tr_broot     = max( 0.0, broot_aim     - broot_in     )
         tr_bsapwooda = max( 0.0, bsapwooda_aim - bsapwooda_in )
         tr_bsapwoodb = max( 0.0, bsapwoodb_aim - bsapwoodb_in )
         tr_bbarka    = max( 0.0, bbarka_aim    - bbarka_in    )
         tr_bbarkb    = max( 0.0, bbarkb_aim    - bbarkb_in    )
         tr_bdeada    = max( 0.0, bdeada_aim    - bdeada_in    )
         tr_bdeadb    = max( 0.0, bdeadb_aim    - bdeadb_in    )
         tr_every     = tr_bleaf  + tr_broot  + tr_bsapwooda + tr_bsapwoodb                &
                      + tr_bbarka + tr_bbarkb + tr_bdeada    + tr_bdeadb
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Make sure growth is going to be non-negligible.  If not, wait until storage  !
         ! has more to offer.                                                              !
         !---------------------------------------------------------------------------------!
         if ( tr_every < (r_tol_trunc * bevery_in) ) then
            f_bstorage = f_bstorage + f_growth
            f_bdeada   = 0.0
            f_bdeadb   = 0.0
            f_growth   = 0.0
         else
            f_bdeada   = f_growth * tr_bdeada / tr_every
            f_bdeadb   = f_growth * tr_bdeadb / tr_every
            f_bstorage = f_bstorage + f_growth - f_bdeada - f_bdeadb
         end if
         !---------------------------------------------------------------------------------!
      case default
         !----- Allocate growth proportionally to the above/below ground fractions. -------!
         f_bdeada      = f_growth * agf_bs(ipft)
         f_bdeadb      = f_growth * (1.0 - agf_bs(ipft))
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      return
   end subroutine bdead_structural_allocation
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine checks that carbon is conserved at the cohort level.  Minor      !
   ! truncation errors ma cause slightly negative pools.  In this case, we fix them and    !
   ! account for the missed carbon source.  Otherwise, we stop any cohort that is          !
   ! attempting to smuggle or to evade carbon.                                             !
   !---------------------------------------------------------------------------------------!
   subroutine check_bstruct_cohort(csite,ipa,ico,bleaf_in,broot_in,bsapwooda_in            &
                                  ,bsapwoodb_in,bbarka_in,bbarkb_in,balive_in,bstorage_in  &
                                  ,bdeada_in,bdeadb_in,f_bseeds,f_growth,f_bdeada,f_bdeadb &
                                  ,f_bstorage,carbon_miss)
      use ed_state_vars, only : sitetype          & ! structure
                              , patchtype         ! ! structure
      use allometry    , only : size2bl           & ! function
                              , size2bd           ! ! function
      use budget_utils , only : tol_carbon_budget ! ! intent(in)
      use pft_coms     , only : agf_bs            & ! intent(in)
                              , min_dbh           & ! intent(in)
                              , hgt_min           ! ! intent(in)
      use ed_misc_coms , only : current_time      ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target        :: csite
      integer         , intent(in)    :: ipa
      integer         , intent(in)    :: ico
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
      real            , intent(in)    :: f_bseeds
      real            , intent(in)    :: f_growth
      real            , intent(in)    :: f_bdeada
      real            , intent(in)    :: f_bdeadb
      real            , intent(in)    :: f_bstorage
      real            , intent(out)   :: carbon_miss
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer    :: cpatch
      integer                      :: ipft
      real                         :: bleaf_ok_min
      real                         :: bdead_ok_min
      real                         :: bdeada_ok_min
      real                         :: bdeadb_ok_min
      real                         :: bstorage_ok_min
      real                         :: btotal_in
      real                         :: btotal_fn
      real                         :: delta_btotal
      real                         :: resid_btotal
      logical                      :: neg_biomass
      logical                      :: btotal_violation
      !----- Local constants. -------------------------------------------------------------!
      character(len=10), parameter :: fmti='(a,1x,i14)'
      character(len=13), parameter :: fmtf='(a,1x,es14.7)'
      character(len=27), parameter :: fmtt='(a,i4.4,2(1x,i2.2),1x,f6.0)'
      !------------------------------------------------------------------------------------!


      !----- Handy aliases. ---------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      ipft = cpatch%pft(ico)
      !------------------------------------------------------------------------------------!


      !----- First, find the minimum possible scale for each pool. ------------------------!
      bleaf_ok_min     = size2bl(min_dbh(ipft),hgt_min(ipft),cpatch%sla(ico),ipft)
      bdead_ok_min     = size2bd(min_dbh(ipft),hgt_min(ipft),ipft)
      bdeada_ok_min    =     agf_bs(ipft)  * bdead_ok_min
      bdeadb_ok_min    = (1.-agf_bs(ipft)) * bdead_ok_min
      bstorage_ok_min  = bleaf_ok_min
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Then, if possible, set the minimum acceptable deviation based on the input, to !
      ! avoid false alarms due to truncation errors when the pool is much greater than the !
      ! minimum size but loses all its stocks.                                             !
      !------------------------------------------------------------------------------------!
      bdeada_ok_min    = - tol_carbon_budget * max(bdeada_in   , bdeada_ok_min   )
      bdeadb_ok_min    = - tol_carbon_budget * max(bdeadb_in   , bdeadb_ok_min   )
      bstorage_ok_min  = - tol_carbon_budget * max(bstorage_in , bstorage_ok_min )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Check every tissue and storage, to make sure they have sensible numbers.  Tiny  !
      ! negative stocks will be tolerated, but don't fix anything in case any of the pools !
      ! are too negative.                                                                  !
      !------------------------------------------------------------------------------------!
      neg_biomass    = cpatch%bdeada   (ico) < bdeada_ok_min    .or.                       &
                       cpatch%bdeadb   (ico) < bdeadb_ok_min    .or.                       &
                       cpatch%bstorage (ico) < bstorage_ok_min
      if (.not. neg_biomass) then
         !----- Account for any potential violation of carbon stocks. ---------------------!
         carbon_miss = - min(cpatch%bdeada   (ico),0.0) - min(cpatch%bdeadb   (ico),0.0)   &
                       - min(cpatch%bstorage (ico),0.0)
         !---------------------------------------------------------------------------------!


         !----- Make sure that all pools are non-negative. --------------------------------!
         cpatch%bdeada   (ico) = max(cpatch%bdeada   (ico),0.0)
         cpatch%bdeadb   (ico) = max(cpatch%bdeadb   (ico),0.0)
         cpatch%bstorage (ico) = max(cpatch%bstorage (ico),0.0)
         !---------------------------------------------------------------------------------!
      else
         !----- Set missing carbon to zero so the code works with debugging. --------------!
         carbon_miss = 0.0
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Check carbon stocks before and after active tissue growth.                     !
      !------------------------------------------------------------------------------------!
      btotal_in         = balive_in + bdeada_in + bdeadb_in + bstorage_in
      btotal_fn         = cpatch%balive(ico) + cpatch%bdeada  (ico)                        &
                        + cpatch%bdeadb(ico) + cpatch%bstorage(ico)
      delta_btotal      = f_bseeds * bstorage_in
      resid_btotal      = btotal_fn - btotal_in + delta_btotal - carbon_miss
      btotal_violation  = abs(resid_btotal) > (tol_carbon_budget * btotal_in)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case we identify carbon conservation issues, print information on screen    !
      ! and stop the model.                                                                !
      !------------------------------------------------------------------------------------!
      if ( neg_biomass .or. btotal_violation ) then
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|      !!!   Cohort Bstruct budget failed   !!!      |'
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
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
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
         write(unit=*,fmt=fmtf )  ' F_BSEEDS            : ',f_bseeds
         write(unit=*,fmt=fmtf )  ' F_GROWTH            : ',f_growth
         write(unit=*,fmt=fmtf )  ' F_BDEADA            : ',f_bdeada
         write(unit=*,fmt=fmtf )  ' F_BDEADB            : ',f_bdeadb
         write(unit=*,fmt=fmtf )  ' F_BSTORAGE          : ',f_bstorage
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BLEAF_FN            : ',cpatch%bleaf    (ico)
         write(unit=*,fmt=fmtf )  ' BROOT_FN            : ',cpatch%broot    (ico)
         write(unit=*,fmt=fmtf )  ' BSAPWOODA_FN        : ',cpatch%bsapwooda(ico)
         write(unit=*,fmt=fmtf )  ' BSAPWOODB_FN        : ',cpatch%bsapwoodb(ico)
         write(unit=*,fmt=fmtf )  ' BBARKA_FN           : ',cpatch%bbarka   (ico)
         write(unit=*,fmt=fmtf )  ' BBARKB_FN           : ',cpatch%bbarkb   (ico)
         write(unit=*,fmt=fmtf )  ' BDEADA_FN           : ',cpatch%bdeada   (ico)
         write(unit=*,fmt=fmtf )  ' BDEADB_FN           : ',cpatch%bdeadb   (ico)
         write(unit=*,fmt=fmtf )  ' BSEEDS_FN           : ',cpatch%bseeds   (ico)
         write(unit=*,fmt=fmtf )  ' BYIELD_FN           : ',cpatch%byield   (ico)
         write(unit=*,fmt=fmtf )  ' BSTORAGE_FN         : ',cpatch%bstorage (ico)
         write(unit=*,fmt=fmtf )  ' BALIVE_FN           : ',cpatch%balive   (ico)
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
                         ,'check_bstruct_cohort','structural_growth.f90')
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine check_bstruct_cohort
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine checks that carbon is conserved at the patch level.  Minor       !
   ! truncation errors may cause slightly negative pools, which are accounted.  Otherwise, !
   ! we stop any patch that is attempting to smuggle or to evade carbon.                   !
   !---------------------------------------------------------------------------------------!
   subroutine check_bstruct_patch(csite,ipa,fgc_in_in,fsc_in_in,stgc_in_in,stsc_in_in      &
                                 ,pat_balive_in,pat_bdead_in,pat_bstorage_in               &
                                 ,pat_carbon_miss,pat_mortality)
      use ed_state_vars, only : sitetype           & ! structure
                              , patchtype          ! ! structure
      use budget_utils , only : tol_carbon_budget  ! ! intent(in)
      use ed_misc_coms , only : current_time       ! ! intent(in)
      use pft_coms     , only : seedling_mortality ! ! intent(in)
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
      real            , intent(in)    :: pat_carbon_miss
      real            , intent(in)    :: pat_mortality
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer       :: cpatch
      integer                         :: ico
      integer                         :: ipft
      real                            :: fgc_in_fn
      real                            :: fsc_in_fn
      real                            :: stgc_in_fn
      real                            :: stsc_in_fn
      real                            :: pat_balive_fn
      real                            :: pat_bdead_fn
      real                            :: pat_bstorage_fn
      real                            :: pat_bseeds_fn
      real                            :: pat_byield_fn
      real                            :: pat_btotal_in
      real                            :: pat_btotal_fn
      real                            :: pat_biomass_in
      real                            :: pat_biomass_fn
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
      fsc_in_fn  =  csite%fsc_in (ipa)
      stgc_in_fn =  csite%stgc_in(ipa)
      stsc_in_fn =  csite%stsc_in(ipa)
      !------------------------------------------------------------------------------------!



      !----- Count current stocks. --------------------------------------------------------!
      pat_balive_fn   = 0.0
      pat_bdead_fn    = 0.0
      pat_bstorage_fn = 0.0
      pat_bseeds_fn   = 0.0
      pat_byield_fn   = 0.0
      do ico=1,cpatch%ncohorts
         ipft            = cpatch%pft(ico)
         pat_balive_fn   = pat_balive_fn   + cpatch%nplant(ico) * cpatch%balive  (ico)
         pat_bdead_fn    = pat_bdead_fn                                                    &
                         + cpatch%nplant(ico) * ( cpatch%bdeada(ico) + cpatch%bdeadb(ico) )
         pat_bstorage_fn = pat_bstorage_fn + cpatch%nplant(ico) * cpatch%bstorage(ico)
         pat_bseeds_fn   = pat_bseeds_fn   + cpatch%nplant(ico) * cpatch%bseeds(ico)       &
                         * (1.0-seedling_mortality(ipft))
         pat_byield_fn   = pat_byield_fn   + cpatch%nplant(ico) * cpatch%byield(ico)
      end do
      !------------------------------------------------------------------------------------!


      !------ Summary of the carbon stocks. -----------------------------------------------!
      pat_biomass_in = pat_balive_in + pat_bdead_in + pat_bstorage_in
      pat_biomass_fn = pat_balive_fn + pat_bdead_fn + pat_bstorage_fn
      soilc_in_in    = fgc_in_in + fsc_in_in + stgc_in_in + stsc_in_in
      soilc_in_fn    = fgc_in_fn + fsc_in_fn + stgc_in_fn + stsc_in_fn
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Check carbon stocks before and after active tissue growth.  Yield is positive  !
      ! in the residual calculation because it is a committed permanent loss of carbon     !
      ! (shipped off-site).                                                                !
      !------------------------------------------------------------------------------------!
      pat_btotal_in        = pat_biomass_in + soilc_in_in
      pat_btotal_fn        = pat_biomass_fn + soilc_in_fn
      resid_pat_btotal     = pat_btotal_fn + pat_bseeds_fn + pat_byield_fn - pat_btotal_in &
                           - pat_carbon_miss
      pat_btotal_violation = abs(resid_pat_btotal) > (tol_carbon_budget * pat_btotal_in)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case we identify carbon conservation issues, print information on screen    !
      ! and stop the model.                                                                !
      !------------------------------------------------------------------------------------!
      if ( pat_btotal_violation ) then
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|       !!!   Patch Bstruct budget failed   !!!      |'
         write(unit=*,fmt='(a)')  '|----------------------------------------------------|'
         write(unit=*,fmt=fmtt )  ' TIME                : ',current_time%year              &
                                                           ,current_time%month             &
                                                           ,current_time%date              &
                                                           ,current_time%time
         write(unit=*,fmt=fmti )  ' PATCH               : ',ipa
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
         write(unit=*,fmt=fmtf )  ' BSEEDS_FN (alive)   : ',pat_bseeds_fn
         write(unit=*,fmt=fmtf )  ' FGC_IN_FN           : ',fgc_in_fn
         write(unit=*,fmt=fmtf )  ' FSC_IN_FN           : ',fsc_in_fn
         write(unit=*,fmt=fmtf )  ' STGC_IN_FN          : ',stgc_in_fn
         write(unit=*,fmt=fmtf )  ' STSC_IN_FN          : ',stsc_in_fn
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BIOMASS_IN          : ',pat_biomass_in
         write(unit=*,fmt=fmtf )  ' SOILC_IN_IN         : ',soilc_in_in
         write(unit=*,fmt=fmtf )  ' BIOMASS_FN          : ',pat_biomass_fn
         write(unit=*,fmt=fmtf )  ' SOILC_IN_FN         : ',soilc_in_fn
         write(unit=*,fmt=fmtf )  ' MORTALITY           : ',pat_mortality
         write(unit=*,fmt=fmtf )  ' YIELD               : ',pat_byield_fn
         write(unit=*,fmt='(a)')  ' ---------------------------------------------------- '
         write(unit=*,fmt=fmtf )  ' BTOTAL_IN           : ',pat_btotal_in
         write(unit=*,fmt=fmtf )  ' BTOTAL_FN           : ',pat_btotal_fn
         write(unit=*,fmt=fmtf )  ' CARBON_MISS         : ',pat_carbon_miss
         write(unit=*,fmt=fmtf )  ' RESIDUAL_BTOTAL     : ',resid_pat_btotal
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  ' '


         call fatal_error('Budget check has failed, see message above.'                    &
                         ,'check_bstruct_patch','structural_growth.f90')
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine check_bstruct_patch
   !=======================================================================================!
   !=======================================================================================!
end module structural_growth
!==========================================================================================!
!==========================================================================================!
