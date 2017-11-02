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
   subroutine dbdead_dt(cgrid, month)
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
                                     , h_edge                      ! ! intent(in)
      use decomp_coms         , only : f_labile                    ! ! intent(in)
      use disturb_coms        , only : cl_fseeds_harvest           ! ! intent(in)
      use ed_max_dims         , only : n_pft                       & ! intent(in)
                                     , n_dbh                       ! ! intent(in)
      use ed_misc_coms        , only : ibigleaf                    & ! intent(in)
                                     , current_time                ! ! intent(in)
      use ed_therm_lib        , only : calc_veg_hcap               & ! function
                                     , update_veg_energy_cweh      ! ! function
      use ed_misc_coms        , only : igrass                      ! ! intent(in)
      use physiology_coms     , only : ddmort_const                & ! intent(in)
                                     , iddmort_scheme              & ! intent(in)
                                     , cbr_scheme                  ! ! intent(in)
      use fuse_fiss_utils     , only : sort_cohorts                ! ! subroutine
      use stable_cohorts      , only : is_resolvable               ! ! subroutine
      use update_derived_utils, only : update_cohort_derived_props & ! subroutine
                                     , update_vital_rates          ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      integer          , intent(in) :: month
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
      integer                       :: imonth
      real                          :: balive_in
      real                          :: bdead_in
      real                          :: bleaf_in
      real                          :: hite_in
      real                          :: dbh_in
      real                          :: nplant_in
      real                          :: bstorage_in
      real                          :: agb_in
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
      real                          :: f_bdead
      real                          :: balive_mort_litter
      real                          :: bstorage_mort_litter
      real                          :: struct_litter
      real                          :: maxh !< maximum patch height
      real                          :: mort_litter
      real                          :: seed_litter
      real                          :: net_seed_N_uptake
      real                          :: net_stem_N_uptake
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      logical          , parameter  :: printout  = .false.
      character(len=17), parameter  :: fracfile  = 'struct_growth.txt'
      !----- Locally saved variables. -----------------------------------------------------!
      logical          , save       :: first_time = .true.
      !------------------------------------------------------------------------------------!


      !----- First time, and the user wants to print the output.  Make a header. ----------!
      if (first_time) then

         !----- Make the header. ----------------------------------------------------------!
         if (printout) then
            open (unit=66,file=fracfile,status='replace',action='write')
            write (unit=66,fmt='(20(a,1x))')                                               &
              ,      '  YEAR',     '  MONTH',      '   DAY',      '   PFT',      '   ICO'  &
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
      prev_month = 1 + modulo(month-2,12)
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


                  !----- Remember inputs in order to calculate increments later on. -------!
                  balive_in   = cpatch%balive  (ico)
                  bdead_in    = cpatch%bdead   (ico)
                  bleaf_in    = cpatch%bleaf   (ico)
                  hite_in     = cpatch%hite    (ico)
                  dbh_in      = cpatch%dbh     (ico)
                  nplant_in   = cpatch%nplant  (ico)
                  bstorage_in = cpatch%bstorage(ico)
                  agb_in      = cpatch%agb     (ico)
                  ba_in       = cpatch%basarea (ico)
                  bag_in      = sum(cpoly%basal_area_growth(ipft,:,isi))
                  bam_in      = sum(cpoly%basal_area_mort(ipft,:,isi))
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


                  !----- Calculate litter owing to mortality. -----------------------------!
                  balive_mort_litter   = - cpatch%balive(ico)   * cpatch%monthly_dndt(ico)
                  bstorage_mort_litter = - cpatch%bstorage(ico) * cpatch%monthly_dndt(ico)
                  struct_litter        = - cpatch%bdead(ico)    * cpatch%monthly_dndt(ico)
                  mort_litter          = balive_mort_litter + bstorage_mort_litter         &
                                       + struct_litter
                  !------------------------------------------------------------------------!



                  !----- Reset monthly_dndt. ----------------------------------------------!
                  cpatch%monthly_dndt  (ico) = 0.0
                  cpatch%monthly_dlnndt(ico) = 0.0
                  !------------------------------------------------------------------------!


                  !----- Determine how to distribute what is in bstorage. -----------------!
                  call plant_structural_allocation(cpatch%pft(ico),cpatch%hite(ico)        &
                                                  ,cpatch%dbh(ico),cgrid%lat(ipy)          &
                                                  ,cpatch%phenology_status(ico)            &
                                                  ,bdead_in, bstorage_in                   &
                                                  ,f_bseeds,f_bdead, maxh)
                  !------------------------------------------------------------------------!


                  !----- Grow plants; bdead gets fraction f_bdead of bstorage. ------------!
                  cpatch%bdead(ico) = cpatch%bdead(ico) + f_bdead * cpatch%bstorage(ico)
                  !------------------------------------------------------------------------!


                  if (ibigleaf == 0 ) then
                     !------ NPP allocation to wood and coarse roots in KgC /m2 -----------!
                     cpatch%today_NPPwood(ico) = agf_bs(ipft)*f_bdead*cpatch%bstorage(ico) &
                                                * cpatch%nplant(ico)
                     cpatch%today_NPPcroot(ico) = (1. - agf_bs(ipft)) * f_bdead            &
                                                * cpatch%bstorage(ico) * cpatch%nplant(ico)
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Rebalance the plant nitrogen uptake considering the actual alloc- !
                  ! ation to structural growth.  This is necessary because c2n_stem does   !
                  ! not necessarily equal c2n_storage.                                     !
                  !------------------------------------------------------------------------!
                  net_stem_N_uptake = (cpatch%bdead(ico) - bdead_in) * cpatch%nplant(ico)  &
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
                  seed_litter        = cpatch%bseeds(ico) * cpatch%nplant(ico)             &
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
                  cpatch%bstorage(ico) = cpatch%bstorage(ico) * (1.0 - f_bdead - f_bseeds)
                  !------------------------------------------------------------------------!



                  !----- Finalize litter inputs. ------------------------------------------!
                  csite%fsc_in(ipa) = csite%fsc_in(ipa)                                    &
                                    + f_labile(ipft) * balive_mort_litter                  &
                                    + bstorage_mort_litter + seed_litter
                  csite%fsn_in(ipa) = csite%fsn_in(ipa)                                    &
                                    + f_labile(ipft) * balive_mort_litter / c2n_leaf(ipft) &
                                    + bstorage_mort_litter/ c2n_storage                    &
                                    + seed_litter / c2n_recruit(ipft)
                  csite%ssc_in(ipa) = csite%ssc_in(ipa) + struct_litter                    &
                                    + (1.0 - f_labile(ipft)) * balive_mort_litter
                  csite%ssl_in(ipa) = csite%ssl_in(ipa)                                    &
                                    + ( (1.0 - f_labile(ipft)) * balive_mort_litter        &
                                      + struct_litter ) * l2n_stem                         &
                                    / c2n_stem(cpatch%pft(ico))
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
                  ! MLO. We now update the heat capacity and the vegetation internal       !
                  !      energy.  Since no energy or water balance is done here, we simply !
                  !      update the energy in order to keep the same temperature and water !
                  !      as before.  Internal energy is an extensive variable, we just     !
                  !      account for the difference in the heat capacity to update it.     !
                  !------------------------------------------------------------------------!
                  old_leaf_hcap = cpatch%leaf_hcap(ico)
                  old_wood_hcap = cpatch%wood_hcap(ico)
                  call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdead(ico)                   &
                                    ,cpatch%bsapwooda(ico),cpatch%bbark(ico)               &
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
   end subroutine dbdead_dt
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will compute the seed allocation and carbon balance stuff, but it !
   ! won't apply to cohorts.                                                               !
   ! IMPORTANT: Do not change the order of operations below unless you know what you are   !
   !            doing.  Changing the order can affect the C/N budgets.                     !
   !---------------------------------------------------------------------------------------!
   subroutine dbdead_dt_eq_0(cgrid, month)
      use ed_state_vars       , only : edtype                      & ! structure
                                     , polygontype                 & ! structure
                                     , sitetype                    & ! structure
                                     , patchtype                   ! ! structure
      use pft_coms            , only : seedling_mortality          & ! intent(in)
                                     , c2n_storage                 & ! intent(in)
                                     , c2n_recruit                 & ! intent(in)
                                     , c2n_stem                    & ! intent(in)
                                     , agf_bs                      & ! intent(in)
                                     , cbr_severe_stress           & ! intent(in)
                                     , is_grass                    ! ! intent(in)
      use ed_max_dims         , only : n_pft                       & ! intent(in)
                                     , n_dbh                       ! ! intent(in)
      use ed_therm_lib        , only : calc_veg_hcap               & ! function
                                     , update_veg_energy_cweh      ! ! function
      use ed_misc_coms        , only : igrass                      ! ! intent(in)
      use physiology_coms     , only : ddmort_const                & ! intent(in)
                                     , iddmort_scheme              & ! intent(in)
                                     , cbr_scheme                  ! ! intent(in)
      use stable_cohorts      , only : is_resolvable               ! ! intent(in)
      use update_derived_utils, only : update_cohort_derived_props & ! subroutine
                                     , update_vital_rates          ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      integer          , intent(in) :: month
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
      integer                       :: imonth
      real                          :: cb_act
      real                          :: cb_lightmax
      real                          :: cb_moistmax
      real                          :: cb_mlmax
      real                          :: cbr_light
      real                          :: cbr_moist
      real                          :: cbr_ml
      real                          :: balive_in
      real                          :: bleaf_in
      real                          :: broot_in
      real                          :: bsapwooda_in
      real                          :: bsapwoodb_in
      real                          :: bbark_in
      real                          :: bdead_in
      real                          :: hite_in
      real                          :: dbh_in
      real                          :: nplant_in
      real                          :: bstorage_in
      real                          :: agb_in
      real                          :: ba_in
      real                          :: phenstatus_in
      real                          :: lai_in
      real                          :: wai_in
      real                          :: cai_in
      integer                       :: krdepth_in
      real                          :: f_bseeds
      real                          :: f_bdead
      real                          :: balive_mort_litter
      real                          :: bstorage_mort_litter
      real                          :: struct_litter
      real                          :: mort_litter
      real                          :: seed_litter
      real                          :: net_seed_N_uptake
      real                          :: net_stem_N_uptake
      real                          :: old_leaf_hcap
      real                          :: old_wood_hcap
      !------------------------------------------------------------------------------------!

      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         !----- Initialization. -----------------------------------------------------------!
         cpoly%basal_area(:,:,:) = 0.0
         cpoly%agb(:,:,:)        = 0.0

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               cohortloop: do ico = 1,cpatch%ncohorts
                  !----- Assigning an alias for PFT type. ---------------------------------!
                  ipft    = cpatch%pft(ico)


                  !------------------------------------------------------------------------!
                  !      Remember inputs in order to calculate increments and revert back  !
                  ! to these values later on.                                              !
                  !------------------------------------------------------------------------!
                  balive_in     = cpatch%balive          (ico)
                  bdead_in      = cpatch%bdead           (ico)
                  bleaf_in      = cpatch%bleaf           (ico)
                  broot_in      = cpatch%broot           (ico)
                  bsapwooda_in  = cpatch%bsapwooda       (ico)
                  bsapwoodb_in  = cpatch%bsapwoodb       (ico)
                  bbark_in      = cpatch%bbark           (ico)
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
                  !------------------------------------------------------------------------!
                  !    Flush time derivatives of population .                              !
                  !------------------------------------------------------------------------!
                  cpatch%monthly_dndt  (ico) = 0.0
                  cpatch%monthly_dlnndt(ico) = 0.0
                  !------------------------------------------------------------------------!


                  !----- Calculate litter owing to mortality. -----------------------------!
                  balive_mort_litter   = - cpatch%balive(ico)   * cpatch%monthly_dndt(ico)
                  bstorage_mort_litter = - cpatch%bstorage(ico) * cpatch%monthly_dndt(ico)
                  struct_litter        = - cpatch%bdead(ico)    * cpatch%monthly_dndt(ico)
                  mort_litter          = balive_mort_litter + bstorage_mort_litter         &
                                       + struct_litter

                  !----- Determine how to distribute what is in bstorage. -----------------!
                  call plant_structural_allocation(cpatch%pft(ico),cpatch%hite(ico)        &
                                                  ,cpatch%dbh(ico),cgrid%lat(ipy)          &
                                                  ,cpatch%phenology_status(ico)            &
                                                  ,bdead_in, bstorage_in                   &
                                                  ,f_bseeds,f_bdead, cpatch%hite(1))
                  !------------------------------------------------------------------------!



                  !----- Grow plants; bdead gets fraction f_bdead of bstorage. ------------!
                  cpatch%bdead(ico) = cpatch%bdead(ico) + f_bdead * cpatch%bstorage(ico)


                  !------ NPP allocation to wood and coarse roots in KgC /m2 --------------!
                  cpatch%today_NPPwood(ico) = agf_bs(ipft) * f_bdead * cpatch%bstorage(ico)&
                                             * cpatch%nplant(ico)
                  cpatch%today_NPPcroot(ico) = (1. - agf_bs(ipft)) * f_bdead               &
                                             * cpatch%bstorage(ico) * cpatch%nplant(ico)

                  !------------------------------------------------------------------------!
                  !      Rebalance the plant nitrogen uptake considering the actual alloc- !
                  ! ation to structural growth.  This is necessary because c2n_stem does   !
                  ! not necessarily equal c2n_storage.                                     !
                  !------------------------------------------------------------------------!
                  net_stem_N_uptake = (cpatch%bdead(ico) - bdead_in) * cpatch%nplant(ico)  &
                                    * ( 1.0 / c2n_stem(cpatch%pft(ico)) - 1.0 / c2n_storage)

                  !------------------------------------------------------------------------!
                  !      Calculate total seed production and seed litter.  The seed pool   !
                  ! gets a fraction f_bseeds of bstorage.                                  !
                  !------------------------------------------------------------------------!
                  cpatch%bseeds(ico) = f_bseeds * cpatch%bstorage(ico)

                  cpatch%today_NPPseeds(ico) = f_bseeds * cpatch%bstorage(ico)             &
                                             * cpatch%nplant(ico)

                  !------------------------------------------------------------------------!
                  ! ALS. If agriculture: set seedling_mortality very low or zero           !
                  !      to keep all of the seeds for harvest later in the season          !
                  !------------------------------------------------------------------------!
                  seed_litter        = cpatch%bseeds(ico) * cpatch%nplant(ico)             &
                                     * seedling_mortality(ipft)


                  !------------------------------------------------------------------------!
                  !      Rebalance the plant nitrogen uptake considering the actual alloc- !
                  ! ation to seeds.  This is necessary because c2n_recruit does not have   !
                  ! to be equal to c2n_storage.                                            !
                  !------------------------------------------------------------------------!
                  net_seed_N_uptake = cpatch%bseeds(ico) * cpatch%nplant(ico)              &
                                    * (1.0 / c2n_recruit(ipft) - 1.0 / c2n_storage)
                  !------------------------------------------------------------------------!



                  !----- Decrement the storage pool. --------------------------------------!
                  cpatch%bstorage(ico) = cpatch%bstorage(ico) * (1.0 - f_bdead - f_bseeds)
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
                  if (month == 1) then
                     prev_month = 12
                  else
                     prev_month = month - 1
                  end if
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
                     !------ Storage is accounted. ----------------------------------------!
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

                  !----- Compute the relative carbon balance. -----------------------------!
                  if (is_grass(ipft).and. igrass==1) then
                     !----- Grass loop, use past month's carbon balance only. -------------!
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

                  !----- Soil moisture-related carbon balance. ----------------------------!
                  if (cb_moistmax > 0.0) then
                     cbr_moist = min(1.0, cb_act / cb_moistmax )
                  else
                     cbr_moist = cbr_severe_stress(ipft)
                  end if

                  !----- Soil moisture+light related carbon balance. ----------------------!
                  if (cb_mlmax > 0.0) then
                     cbr_ml    = min(1.0, cb_act / cb_mlmax )
                  else
                     cbr_ml    = cbr_severe_stress(ipft)
                  end if

                  !------------------------------------------------------------------------!
                  !  calculate CBR according to the specified CBR_SCHEME                   !
                  !------------------------------------------------------------------------!
                  select case (cbr_scheme)
                  case (0)
                    !----- CBR from absolute max CB ---------------------------------------!
                    cpatch%cbr_bar(ico) = max(cbr_ml, cbr_severe_stress(ipft))

                  case (1)
                    !----- CBR from combination of light & moist CBR ----------------------!
                    !----- Relative carbon balance: a combination of the two factors. -----!
                    if ( cbr_light <= cbr_severe_stress(ipft) .and.                        &
                      cbr_moist <= cbr_severe_stress(ipft)       ) then
                      cpatch%cbr_bar(ico) = cbr_severe_stress(ipft)
                    else
                      cpatch%cbr_bar(ico) = cbr_severe_stress(ipft)                        &
                              + ( cbr_light - cbr_severe_stress(ipft) )                    &
                              * ( cbr_moist - cbr_severe_stress(ipft) )                    &
                              / (        ddmort_const  * cbr_moist                         &
                                + (1.0 - ddmort_const) * cbr_light                         &
                                - cbr_severe_stress(ipft) )
                    end if

                  case (2)
                    !----- CBR from most limiting CBR -------------------------------------!
                    cpatch%cbr_bar(ico) = max(min(cbr_moist, cbr_light)                    &
                                             , cbr_severe_stress(ipft))

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
                  !     Revert back to previous values:                                    !
                  !------------------------------------------------------------------------!
                  cpatch%balive          (ico) = balive_in
                  cpatch%bdead           (ico) = bdead_in
                  cpatch%bleaf           (ico) = bleaf_in
                  cpatch%broot           (ico) = broot_in
                  cpatch%bsapwooda       (ico) = bsapwooda_in
                  cpatch%bsapwoodb       (ico) = bsapwoodb_in
                  cpatch%bbark           (ico) = bbark_in
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
                  call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdead(ico)                   &
                                    ,cpatch%bsapwooda(ico),cpatch%bbark(ico)               &
                                    ,cpatch%nplant(ico),cpatch%pft(ico)                    &
                                    ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico) )
                  call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap)
                  call is_resolvable(csite,ipa,ico)
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!


      return
   end subroutine dbdead_dt_eq_0
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will decide the partition of storage biomass into seeds and dead  !
   ! (structural) biomass.                                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine plant_structural_allocation(ipft,hite,dbh,lat,phen_status,bdead,bstorage     &
                                         ,f_bseeds,f_bdead,maxh)
      use pft_coms      , only : phenology    & ! intent(in)
                               , repro_min_h  & ! intent(in)
                               , r_fract      & ! intent(in)
                               , st_fract     & ! intent(in)
                               , dbh_crit     & ! intent(in)
                               , hgt_max      & ! intent(in)
                               , is_grass     & ! intent(in)
                               , is_liana     ! ! intent(in)
      use ed_misc_coms  , only : current_time & ! intent(in)
                               , igrass       & ! intent(in)
                               , ibigleaf     ! ! intent(in)
      use consts_coms   , only : r_tol_trunc  ! ! intent(in)
      use allometry     , only : dbh2bd       & ! intent(in)
                               , h2dbh        ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      integer, intent(in)  :: ipft
      real   , intent(in)  :: hite
      real   , intent(in)  :: dbh
      real   , intent(in)  :: lat
      real   , intent(in)  :: bdead     !> Current dead biomass
      real   , intent(in)  :: bstorage  !> Current storage pool
      real   , intent(in)  :: maxh      !> Height of the tallest cohort in the patch
      integer, intent(in)  :: phen_status
      real   , intent(out) :: f_bseeds
      real   , intent(out) :: f_bdead
      !----- Local variables --------------------------------------------------------------!
      real                           :: bd_target !> Target Bd to reach maxh height
      real                           :: delta_bd  !> Target Bd - actual Bd
      logical                        :: late_spring
      logical          , parameter   :: printout  = .false.
      character(len=13), parameter   :: fracfile  = 'storalloc.txt'
      !----- Locally saved variables. -----------------------------------------------------!
      logical          , save        :: first_time = .true.
      !------------------------------------------------------------------------------------!


      !----- First time, and the user wants to print the output.  Make a header. ----------!
      if (first_time) then

         !----- Make the header. ----------------------------------------------------------!
         if (printout) then
            open (unit=66,file=fracfile,status='replace',action='write')
            write (unit=66,fmt='(15(a,1x))')                                               &
              ,'        YEAR','       MONTH','         DAY','         PFT','   PHENOLOGY'  &
              ,' PHEN_STATUS',' LATE_SPRING','       GRASS','      HEIGHT','   REPRO_HGT'  &
              ,'         DBH','    DBH_CRIT','   F_STORAGE','     F_SEEDS','     F_BDEAD'
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


      select case (ibigleaf)
      case (0)
         !---------------------------------------------------------------------------------!
         !      Size and age structure.  Calculate fraction of bstorage going to bdead and !
         ! reproduction.  First we must make sure that the plant should do something here. !
         ! A plant should not allocate anything to reproduction or growth if it is not the !
         ! right time of year (for cold deciduous plants), or if the plants are actively   !
         ! dropping leaves or off allometry.                                               !
         !---------------------------------------------------------------------------------!
         if ((phenology(ipft) /= 2   .or.  late_spring) .and. phen_status == 0)    then
            !------------------------------------------------------------------------------!
            !      This is where allocation to seeds is occuring.  It will need to be      !
            ! modified but I'm leaving it for later --- GRASSES!  Want to add a functional !
            ! form to constrain this throughout the season - also consider moving this to  !
            ! growth_balive since it isn't actually structural growth                      !
            !------------------------------------------------------------------------------!
            if (is_grass(ipft) .and. igrass == 1) then
               !----- New grasses. --------------------------------------------------------!
               if ( hite >= ( (1.0-r_tol_trunc) * hgt_max(ipft)) ) then 
                  !------------------------------------------------------------------------!
                  !   Grasses have reached the maximum height, stop growing in size and    !
                  ! send everything to reproduction.                                       !
                  !------------------------------------------------------------------------!
                  f_bseeds = 1.0 - st_fract(ipft)
                  !------------------------------------------------------------------------!
               elseif (hite < ( (1.0-r_tol_trunc) * repro_min_h(ipft) ) ) then
                  !----- The plant is too short, invest as much as it can in growth. ------!
                  f_bseeds = 0.0
                  !------------------------------------------------------------------------!
               else ! repro_min_h < hite< hgt_max
                  !----- Medium-sized grass, use prescribed reproduction rate. ------------!
                  f_bseeds = r_fract(ipft)
                  !------------------------------------------------------------------------!
               end if
               f_bdead  = 0.0
               !---------------------------------------------------------------------------!

            elseif (is_liana(ipft)) then
               !---------------------------------------------------------------------------!
               !    Lianas: we must check height relative to the rest of the local plant   !
               ! community.                                                                !
               !---------------------------------------------------------------------------!
               if (hite >= ( (1.0-r_tol_trunc) *  maxh)) then
                  f_bseeds = merge(0.0, r_fract(ipft)                                      &
                                  ,hite < ( (1.0-r_tol_trunc) * repro_min_h(ipft)) )
                  f_bdead  = 0.0
               else
                  bd_target = dbh2bd(h2dbh(maxh,ipft),ipft)
                  delta_bd = bd_target - bdead
                  !------------------------------------------------------------------------!
                  !    If bstorage is 0 or lianas have already reached their bd_target     !
                  ! don't grow otherwise invest what is needed (or everything in case it's !
                  ! not enough) to reach bd_target.                                        !
                  !    For seeds first check if the cohort has reached repro_min_h. In     !
                  ! case so, then check that it has enough left to invest r_fract in       !
                  ! reproduction.  In case so, invest r_fract in reproduction and the rest !
                  ! will stay in storage.  Otherwise, invest everything in reproduction.   !
                  ! Finally in case the liana hasn't reached repro_min_h, leave everything !
                  ! in storage.    !                                                       !
                  !------------------------------------------------------------------------!
                  f_bdead   = merge(0.0                                                    &
                                   ,min(delta_bd / bstorage, 1.0)                          &
                                   ,bstorage * delta_bd <= 0.0)
                  f_bseeds  = merge( 0.0, merge( r_fract(ipft), 1.0 - f_bdead              &
                                               , r_fract(ipft) < 1.0 - f_bdead)            &
                                   , hite <= repro_min_h(ipft))
               end if
               !---------------------------------------------------------------------------!


            elseif (hite < ((1.0-r_tol_trunc) * repro_min_h(ipft))) then
               !----- The tree is too short, invest as much as it can in growth. ----------!
               f_bseeds = 0.0
               f_bdead  = 1.0 - st_fract(ipft) - f_bseeds 
               !---------------------------------------------------------------------------!
            else
               !----- Medium-sized tree, use prescribed reproduction rate. ----------------!
               f_bseeds = r_fract(ipft)
               f_bdead  = 1.0 - st_fract(ipft) - f_bseeds 
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         else  
            !----- Plant should not allocate carbon to seeds or grow new biomass. ---------!
            f_bdead  = 0.0
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
         if ((phenology(ipft) /= 2   .or.  late_spring) .and. phen_status == 0)    then
            !------------------------------------------------------------------------------!
            ! A plant should only grow if it is the right time of year (for cold deciduous !
            ! plants), or if the plants are not actively dropping leaves or off allometry. !
            !------------------------------------------------------------------------------!
            f_bseeds = 1.0 - st_fract(ipft)
            f_bdead  = 0.0
         else
            f_bdead  = 0.0
            f_bseeds = 0.0
         end if
      end select
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      if (printout) then
         open (unit=66,file=fracfile,status='old',position='append',action='write')
         write (unit=66,fmt='(6(i12,1x),2(11x,l1,1x),7(f12.4,1x))')                        &
               current_time%year,current_time%month,current_time%date,ipft,phenology(ipft) &
              ,phen_status,late_spring,is_grass(ipft),hite,repro_min_h(ipft),dbh           &
              ,dbh_crit(ipft),st_fract(ipft),f_bseeds,f_bdead
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
