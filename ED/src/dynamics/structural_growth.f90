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
   subroutine dbstruct_dt(cgrid,veget_dyn_on)
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
                                     , is_liana                    & ! intent(in)
                                     , cbr_severe_stress           & ! intent(in)
                                     , h_edge                      & ! intent(in)
                                     , f_labile_leaf               & ! intent(in)
                                     , f_labile_stem               ! ! intent(in)
      use disturb_coms        , only : cl_fseeds_harvest           ! ! intent(in)
      use ed_max_dims         , only : n_pft                       & ! intent(in)
                                     , n_dbh                       ! ! intent(in)
      use ed_misc_coms        , only : iallom                      & ! intent(in)
                                     , igrass                      & ! intent(in)
                                     , ibigleaf                    & ! intent(in)
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
      use allometry           , only : size2bl                     & ! function
                                     , expand_bevery               ! ! subroutine
      use consts_coms         , only : yr_sec                      & ! intent(in)
                                     , almost_zero                 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      logical          , intent(in) :: veget_dyn_on
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
      real                          :: agb_in
      real                          :: lai_in
      real                          :: wai_in
      real                          :: cai_in
      real                          :: ba_in
      real                          :: bag_in
      real                          :: bam_in
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
      real                          :: dbh_aim
      real                          :: hite_aim
      real                          :: bleaf_aim
      real                          :: broot_aim
      real                          :: bsapwooda_aim
      real                          :: bsapwoodb_aim
      real                          :: bbarka_aim
      real                          :: bbarkb_aim
      real                          :: balive_aim
      real                          :: bdeada_aim
      real                          :: bdeadb_aim
      real                          :: bevery_aim
      real                          :: maxh !< maximum patch height
      real                          :: mort_litter
      real                          :: bseeds_mort_litter
      real                          :: net_seed_N_uptake
      real                          :: net_stem_N_uptake
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
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


        !----- Initialization. ------------------------------------------------------------!
        cgrid%crop_yield(prev_month,ipy) = 0.0

         

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)
            !----- Initialization. --------------------------------------------------------!
            cpoly%basal_area(:,:,isi)       = 0.0
            cpoly%agb(:,:,isi)               = 0.0
            cpoly%crop_yield(prev_month,isi) = 0.0

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

               cohortloop: do ico = 1,cpatch%ncohorts

                  !----- Assigning an alias for PFT type. ---------------------------------!
                  ipft    = cpatch%pft(ico)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      Remember inputs in order to calculate increments (and to revert   !
                  ! back to these values later on in case vegetation dynamics is off).     !
                  !------------------------------------------------------------------------!
                  balive_in     = cpatch%balive          (ico)
                  bdeada_in     = cpatch%bdeada          (ico)
                  bdeadb_in     = cpatch%bdeadb          (ico)
                  bleaf_in      = cpatch%bleaf           (ico)
                  broot_in      = cpatch%broot           (ico)
                  bsapwooda_in  = cpatch%bsapwooda       (ico)
                  bsapwoodb_in  = cpatch%bsapwoodb       (ico)
                  bbarka_in     = cpatch%bbarka          (ico)
                  bbarkb_in     = cpatch%bbarkb          (ico)
                  hite_in       = cpatch%hite            (ico)
                  dbh_in        = cpatch%dbh             (ico)
                  nplant_in     = cpatch%nplant          (ico)
                  bstorage_in   = cpatch%bstorage        (ico)
                  agb_in        = cpatch%agb             (ico)
                  ba_in         = cpatch%basarea         (ico)
                  phenstatus_in = cpatch%phenology_status(ico)
                  lai_in        = cpatch%lai             (ico)
                  wai_in        = cpatch%wai             (ico)
                  cai_in        = cpatch%crown_area      (ico)
                  krdepth_in    = cpatch%krdepth         (ico)
                  bevery_in     = bleaf_in  + broot_in  + bsapwooda_in + bsapwoodb_in      &
                                + bbarka_in + bbarkb_in + bdeada_in    + bdeadb_in
                  bag_in        = sum(cpoly%basal_area_growth(ipft,:,isi))
                  bam_in        = sum(cpoly%basal_area_mort(ipft,:,isi))
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !    Apply mortality, and do not allow nplant < negligible_nplant (such  !
                  ! a sparse cohort is about to be terminated, anyway).                    !
                  ! NB: monthly_dndt may be negative.                                      !
                  !------------------------------------------------------------------------!
                  cpatch%monthly_dndt  (ico) = max( cpatch%monthly_dndt   (ico)            &
                                                  , negligible_nplant     (ipft)           &
                                                  - cpatch%nplant         (ico) )
                  cpatch%monthly_dlnndt(ico) = max( cpatch%monthly_dlnndt (ico)            &
                                                  , log( negligible_nplant(ipft)           &
                                                       / cpatch%nplant    (ico) ) )
                  cpatch%nplant(ico)         = cpatch%nplant(ico)                          &
                                             * exp(cpatch%monthly_dlnndt(ico))
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


                  !----- Calculate litter owing to mortality. -----------------------------!
                  a_bfast_mort_litter    = - a_bfast    * cpatch%monthly_dndt(ico)
                  b_bfast_mort_litter    = - b_bfast    * cpatch%monthly_dndt(ico)
                  a_bstruct_mort_litter  = - a_bstruct  * cpatch%monthly_dndt(ico)
                  b_bstruct_mort_litter  = - b_bstruct  * cpatch%monthly_dndt(ico)
                  a_bstorage_mort_litter = - a_bstorage * cpatch%monthly_dndt(ico)
                  b_bstorage_mort_litter = - b_bstorage * cpatch%monthly_dndt(ico)
                  mort_litter            = a_bfast_mort_litter    + b_bfast_mort_litter    &
                                         + a_bstruct_mort_litter  + b_bstruct_mort_litter  &
                                         + a_bstorage_mort_litter + b_bstorage_mort_litter
                  !------------------------------------------------------------------------!



                  !----- Reset monthly_dndt. ----------------------------------------------!
                  cpatch%monthly_dndt  (ico) = 0.0
                  cpatch%monthly_dlnndt(ico) = 0.0
                  !------------------------------------------------------------------------!


                  !----- Determine how to distribute what is in bstorage. -----------------!
                  call plant_structural_allocation(cpatch%pft(ico),cpatch%hite(ico)        &
                                                  ,cpatch%dbh(ico),cgrid%lat(ipy)          &
                                                  ,cpatch%phenology_status(ico)            &
                                                  ,cpatch%elongf(ico)                      &
                                                  ,bdeada_in,bdeadb_in, bstorage_in, maxh  &
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
                  if (f_growth > almost_zero) then
                     select case (iallom)
                     case (3)
                        bevery_aim     = bevery_in + f_growth * cpatch%bstorage(ico)
                        call expand_bevery(cpatch%pft(ico),bevery_aim,dbh_aim,hite_aim     &
                                          ,bleaf_aim,broot_aim,bsapwooda_aim,bsapwoodb_aim &
                                          ,bbarka_aim,bbarkb_aim,balive_aim,bdeada_aim     &
                                          ,bdeadb_aim)
                        if (bevery_aim > bevery_in .and. bdeada_aim > bdeada_in) then
                           f_bdeada    = f_growth * ( bdeada_aim - bdeada_in )             &
                                                  / ( bevery_aim - bevery_in )
                        else
                           f_bdeada    = 0.0
                        end if
                        if (bevery_aim > bevery_in .and. bdeadb_aim > bdeadb_in) then
                           f_bdeadb    = f_growth * ( bdeadb_aim - bdeadb_in )             &
                                                  / ( bevery_aim - bevery_in )
                        else
                           f_bdeadb    = 0.0
                        end if
                        f_bstorage     = f_bstorage + f_growth - f_bdeada - f_bdeadb
                     case default
                        f_bdeada       = f_growth * agf_bs(ipft)
                        f_bdeadb       = f_growth * (1.0 - agf_bs(ipft))
                     end select
                  else
                     f_bstorage  = f_bstorage + f_growth
                     f_growth    = 0.
                     f_bdeada    = 0.
                     f_bdeadb    = 0.
                  end if
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
                  !      Send dead seeds to the litter pool (this time, only those that    !
                  ! were not harvested).                                                   !
                  !------------------------------------------------------------------------!
                  bseeds_mort_litter = cpatch%bseeds(ico) * cpatch%nplant(ico)             &
                                     * seedling_mortality(ipft)
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



                  !------ Update crop yield. ----------------------------------------------!
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
                  call update_cohort_derived_props(cpatch,ico,cpoly%lsl(isi))
                  !------------------------------------------------------------------------!


                  !----- Update annual average carbon balances for mortality. -------------!
                  cpatch%cb          (prev_month,ico) = cpatch%cb          (13,ico)
                  cpatch%cb_lightmax (prev_month,ico) = cpatch%cb_lightmax (13,ico)
                  cpatch%cb_moistmax (prev_month,ico) = cpatch%cb_moistmax (13,ico)
                  cpatch%cb_mlmax    (prev_month,ico) = cpatch%cb_mlmax    (13,ico)
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
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  ! MLO. We now update the heat capacity and the vegetation internal       !
                  !      energy.  Since no energy or water balance is done here, we simply !
                  !      update the energy in order to keep the same temperature and water !
                  !      as before.  Internal energy is an extensive variable, we just     !
                  !      account for the difference in the heat capacity to update it.     !
                  !------------------------------------------------------------------------!
                  old_leaf_hcap = cpatch%leaf_hcap(ico)
                  old_wood_hcap = cpatch%wood_hcap(ico)
                  call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdeada(ico)                  &
                                    ,cpatch%bsapwooda(ico),cpatch%bbarka(ico)              &
                                    ,cpatch%nplant(ico),cpatch%pft(ico)                    &
                                    ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
                  call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap)
                  call is_resolvable(csite,ipa,ico)
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

               end do cohortloop
               !---------------------------------------------------------------------------!

               !----- Age the patch if this is not agriculture. ---------------------------!
               if (csite%dist_type(ipa) /= 1) csite%age(ipa) = csite%age(ipa) + 1.0/12.0
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
                                         ,bdeadb,bstorage,maxh,f_bseeds,f_growth,f_bstorage)
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
      f_bstorage = 1.0 - f_growth - f_bseeds
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
end module structural_growth
!==========================================================================================!
!==========================================================================================!
