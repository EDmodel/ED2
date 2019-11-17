!==========================================================================================!
!==========================================================================================!
!     This module contains the main reproduction sub-routines.                             !
!------------------------------------------------------------------------------------------!
module reproduction
   contains

   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine will drive the reproduction, based on its carbon availability    !
   ! and PFT-specific reproduction properties.  No reproduction will happen if the user    !
   ! didn't want it, in which case the seedling biomass will go to the litter pools.       !
   !---------------------------------------------------------------------------------------!
   subroutine reproduction_driver(cgrid,month,veget_dyn_on)
      use ed_state_vars       , only : edtype                      & ! structure
                                     , polygontype                 & ! structure
                                     , sitetype                    & ! structure
                                     , patchtype                   & ! structure
                                     , allocate_patchtype          & ! subroutine
                                     , copy_patchtype              & ! subroutine
                                     , deallocate_patchtype        ! ! subroutine
      use met_driver_coms     , only : met_driv_state              ! ! structure
      use pft_coms            , only : recruittype                 & ! structure
                                     , zero_recruit                & ! subroutine
                                     , copy_recruit                & ! subroutine
                                     , seedling_mortality          & ! intent(in)
                                     , min_recruit_size            & ! intent(in)
                                     , one_plant_c                 & ! intent(in)
                                     , c2n_recruit                 & ! intent(in)
                                     , include_pft                 & ! intent(in)
                                     , include_pft_ag              & ! intent(in)
                                     , include_pft_fp              & ! intent(in)
                                     , q                           & ! intent(in)
                                     , qsw                         & ! intent(in)
                                     , qbark                       & ! intent(in)
                                     , agf_bs                      & ! intent(in)
                                     , hgt_min                     & ! intent(in)
                                     , plant_min_temp              ! ! intent(in)
      use ed_max_dims         , only : n_pft                       ! ! intent(in)
      use fuse_fiss_utils     , only : sort_cohorts                & ! subroutine
                                     , terminate_cohorts           & ! subroutine
                                     , old_fuse_cohorts            & ! subroutine
                                     , new_fuse_cohorts            & ! subroutine
                                     , split_cohorts               & ! subroutine
                                     , rescale_patches             ! ! subroutine
      use phenology_coms      , only : repro_scheme                ! ! intent(in)
      use mem_polygons        , only : maxcohort                   ! ! intent(in)
      use consts_coms         , only : pio4                        & ! intent(in)
                                     , t3ple                       ! ! intent(in)
      use ed_therm_lib        , only : calc_veg_hcap               & ! function
                                     , update_veg_energy_cweh      ! ! function
      use allometry           , only : size2bl                     & ! function
                                     , size2bd                     & ! function
                                     , h2dbh                       & ! function
                                     , size2bt                     & ! function
                                     , size2xb                     & ! function
                                     , ed_biomass                  & ! function
                                     , area_indices                & ! subroutine
                                     , size2krdepth                ! ! function
      use grid_coms           , only : nzg                         ! ! intent(in)
      use ed_misc_coms        , only : ibigleaf                    & ! intent(in)
                                     , current_time                ! ! intent(in)
      use phenology_aux       , only : pheninit_balive_bstorage    ! ! intent(in)
      use stable_cohorts      , only : is_resolvable               ! ! function
      use update_derived_utils, only : update_patch_derived_props  & ! sub-routine
                                     , update_site_derived_props   & ! sub-routine
                                     , update_cohort_plastic_trait ! ! sub-routine
      use fusion_fission_coms , only : ifusion                     ! ! intent(in)
      use ed_type_init        , only : init_ed_cohort_vars         ! ! sub-routine
      use plant_hydro         , only : rwc2tw                      & ! sub-routine
                                     , twi2twe                     ! ! sub-routine
      use physiology_coms     , only : trait_plasticity_scheme     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)        , target     :: cgrid
      integer             , intent(in) :: month
      logical             , intent(in) :: veget_dyn_on
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)   , pointer          :: cpoly
      type(sitetype)      , pointer          :: csite
      type(met_driv_state), pointer          :: cmet
      type(patchtype)     , pointer          :: cpatch
      type(patchtype)     , pointer          :: temppatch
      type(recruittype)   , dimension(n_pft) :: recruit
      type(recruittype)                      :: rectest
      integer                                :: ipy
      integer                                :: isi
      integer                                :: ipa
      integer                                :: ico
      integer                                :: ipft
      integer                                :: inew
      integer                                :: imon
      integer                                :: ncohorts_new
      logical                                :: late_spring
      logical                                :: allow_pft
      logical                                :: make_recruit
      real                                   :: elim_nplant
      real                                   :: elim_lai
      real                                   :: nplant_inc
      real                                   :: bleaf_plant
      real                                   :: bdead_plant
      real                                   :: broot_plant
      real                                   :: bsapwood_plant
      real                                   :: bbark_plant
      real                                   :: balive_plant
      real                                   :: rec_bdead
      real                                   :: rec_biomass
      logical             , parameter        :: printout  = .false.
      character(len=17)   , parameter        :: fracfile  = 'repro_details.txt'
      !----- Saved variables --------------------------------------------------------------!
      logical             , save             :: first_time = .true.
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    If this is the first time, check whether the user wants reproduction.  If not,  !
      ! kill all potential recruits and send their biomass to the litter pool.             !
      !------------------------------------------------------------------------------------!
      if (first_time) then
         !---- Halt reproduction by killing all seedlings. --------------------------------!
         if (repro_scheme == 0 .or. (.not. veget_dyn_on)) then
            seedling_mortality(1:n_pft) = 1.0
         end if
         !---------------------------------------------------------------------------------!

         !----- Make the header. ----------------------------------------------------------!
         if (printout .and. veget_dyn_on) then
            open (unit=66,file=fracfile,status='replace',action='write')
            write (unit=66,fmt='(11(a,1x))')                                               &
                     '  YEAR',     '  MONTH',      '   DAY',      '   IPA',      '   ICO'  &
              ,      '   PFT','  REC_NPLANT',' REC_BIOMASS','       REPRO','MIN_REC_SIZE'  &
              ,'MAKE_RECRUIT'
            close (unit=66,status='keep')
         end if
         !---------------------------------------------------------------------------------!

         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Make sure repro is set to zero in case vegetation dynamics is off.             !
      !------------------------------------------------------------------------------------!
      if (.not. veget_dyn_on) then
         offpolyloop: do ipy = 1,cgrid%npolygons
            cpoly => cgrid%polygon(ipy)

            offsiteloop: do isi = 1,cpoly%nsites
               csite => cpoly%site(isi)

               !---------------------------------------------------------------------------!
               !    Zero all reproduction stuff...                                         !
               !---------------------------------------------------------------------------!
               offpatchloop: do ipa = 1,csite%npatches
                  call zero_recruit(n_pft,recruit)
                  csite%repro(:,ipa)  = 0.0
               end do offpatchloop
               !---------------------------------------------------------------------------!

               !----- Reset minimum monthly temperature. ----------------------------------!
               cpoly%min_monthly_temp(isi) = huge(1.)
               !---------------------------------------------------------------------------!
            end do offsiteloop
            !------------------------------------------------------------------------------!
         end do offpolyloop
         return
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Decide which vegetation structure to use.                                      !
      !------------------------------------------------------------------------------------!
      select case (ibigleaf)
      case (0)
         !---------------------------------------------------------------------------------!
         !                        Size- and age_structure reproduction                     !
         !---------------------------------------------------------------------------------!


         !----- The big loops start here. -------------------------------------------------!
         polyloop: do ipy = 1,cgrid%npolygons

            !------------------------------------------------------------------------------!
            !     Check whether this is late spring/early summer.  This is needed for      !
            ! temperate broadleaf deciduous trees.  Late spring means June in the Northern !
            ! Hemisphere, or December in the Southern Hemisphere.                          !
            !------------------------------------------------------------------------------!
            late_spring = (cgrid%lat(ipy) >= 0.0 .and. month == 6) .or.                    &
                          (cgrid%lat(ipy) < 0.0 .and. month == 12)

            cpoly => cgrid%polygon(ipy)
            siteloop_sort: do isi = 1,cpoly%nsites
               csite => cpoly%site(isi)

               !---------------------------------------------------------------------------!
               !      Cohorts may have grown differently, so we need to sort them by size. !
               !---------------------------------------------------------------------------!
               patchloop_sort: do ipa = 1,csite%npatches
                  cpatch => csite%patch(ipa)
                  call sort_cohorts(cpatch)
               end do patchloop_sort
            end do siteloop_sort

            !------- Update the repro arrays. ---------------------------------------------!
            call seed_dispersal(cpoly,late_spring)
            !------------------------------------------------------------------------------!

            siteloop: do isi = 1,cpoly%nsites
               csite => cpoly%site(isi)
               cmet  => cpoly%met(isi)

               !---------------------------------------------------------------------------!
               !    For the recruitment to happen, four requirements must be met:          !
               !    1.  PFT is included in this simulation;                                !
               !    2.  It is not too cold (min_monthly_temp > plant_min_temp - 5)         !
               !    3.  We are dealing with EITHER a non-agriculture patch OR              !
               !        a PFT that could exist in an agricultural patch.                   !
               !    4.  There must be sufficient carbon to form the recruits.              !
               !---------------------------------------------------------------------------!
               patchloop: do ipa = 1,csite%npatches
                  inew = 0
                  call zero_recruit(n_pft,recruit)
                  cpatch => csite%patch(ipa)

                  !---- This time we loop over PFTs, not cohorts. -------------------------!
                  pftloop: do ipft = 1, n_pft

                     !---------------------------------------------------------------------!
                     !     Check whether to include this PFT or not.  The decision depends !
                     ! on the following decisions.                                         !
                     ! 1.  The PFT is included in this simulation.   In case of            !
                     !     agriculture or forest plantation, it must also be allowed in    !
                     !     such patches.                                                   !
                     ! 2.  The temperature is not limiting reproduction.                   !
                     ! 3.  The user wants reproduction to occur.                           !
                     !---------------------------------------------------------------------!
                     select case (csite%dist_type(ipa))
                     case (1)
                        !----- Agriculture (cropland or pasture). -------------------------!
                        allow_pft =                                                        &
                           include_pft_ag(ipft)                                      .and. &
                           cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and. &
                           repro_scheme                /= 0
                        !------------------------------------------------------------------!

                     case (2)
                        !----- Forest plantation. -----------------------------------------!
                        allow_pft =                                                        &
                           include_pft_fp(ipft)                                      .and. &
                           cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and. &
                           repro_scheme                /= 0
                        !------------------------------------------------------------------!

                     case default
                        !----- Primary or secondary vegetation. ---------------------------!
                        allow_pft =                                                        &
                           include_pft(ipft)                                         .and. &
                           cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and. &
                           repro_scheme                /= 0
                        !------------------------------------------------------------------!
                     end select
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !     If this PFT is allowed, check for recruitment.                  !
                     !---------------------------------------------------------------------!
                     if (allow_pft) then
                        !------------------------------------------------------------------!
                        !       Initialise recruit, and assume that the new cohorts are in !
                        ! thermal equilibrium with the canopy air space.                   !
                        !------------------------------------------------------------------!
                        rectest%pft          = ipft
                        rectest%leaf_temp    = csite%can_temp (ipa)
                        rectest%wood_temp    = csite%can_temp (ipa)
                        rectest%leaf_temp_pv = csite%can_temp (ipa)
                        rectest%wood_temp_pv = csite%can_temp (ipa)
                        rectest%leaf_vpdef   = csite%can_vpdef(ipa)
                        !------------------------------------------------------------------!

                        !------------------------------------------------------------------!
                        !    Recruits start at minimum height and dbh and bleaf are        !
                        ! calculated from that.                                            !
                        !------------------------------------------------------------------!
                        rectest%hite      = hgt_min(ipft)
                        rectest%dbh       = h2dbh(rectest%hite, ipft)
                        rectest%krdepth   = size2krdepth(rectest%hite,rectest%dbh          &
                                                        ,rectest%pft,cpoly%lsl(isi))
                        rec_bdead         = size2bd(rectest%dbh,rectest%hite,ipft)
                        rectest%bdeada    =         agf_bs(ipft)   * rec_bdead
                        rectest%bdeadb    = ( 1.0 - agf_bs(ipft) ) * rec_bdead

                        call pheninit_balive_bstorage(nzg,rectest%pft,rectest%krdepth      &
                                                     ,rectest%hite,rectest%dbh             &
                                                     ,csite%soil_water(:,ipa)              &
                                                     ,cpoly%ntext_soil(:,isi)              &
                                                     ,rectest%paw_avg,rectest%elongf       &
                                                     ,rectest%phenology_status             &
                                                     ,rectest%bleaf,rectest%broot          &
                                                     ,rectest%bsapwooda,rectest%bsapwoodb  &
                                                     ,rectest%bbarka,rectest%bbarkb        &
                                                     ,rectest%bstorage,rectest%cb          &
                                                     ,rectest%cb_lightmax                  &
                                                     ,rectest%cb_moistmax                  &
                                                     ,rectest%cb_mlmax,rectest%cbr_bar)

                        !------------------------------------------------------------------!
                        !     Find balive: we cannot use the allometry built-in function   !
                        ! because rectest is not a patchtype structure.                    !
                        !------------------------------------------------------------------!
                        rectest%balive = rectest%bleaf     + rectest%broot                 &
                                       + rectest%bsapwooda + rectest%bsapwoodb             &
                                       + rectest%bbarka    + rectest%bbarkb
                        !------------------------------------------------------------------!

                        !------------------------------------------------------------------!
                        !     Find the expected population from the reproduction stocks.   !
                        ! MLO Note: seed rain is now included in csite%repro (Check sub-   !
                        ! routine seed_dispersal).                                         !
                        !------------------------------------------------------------------!
                        rectest%nplant    = csite%repro(ipft,ipa)                          &
                                          / ( rectest%balive                               &
                                            + rectest%bdeada + rectest%bdeadb              &
                                            + rectest%bstorage)
                        !------------------------------------------------------------------!



                        !----- Find the total biomass for this potential new cohort. ------!
                        rec_biomass = rectest%nplant * ( rectest%balive                    &
                                                       + rectest%bdeada + rectest%bdeadb   &
                                                       + rectest%bstorage                )
                        !------------------------------------------------------------------!

                        !----- Decide whether it is possible to build a recruit. ----------!
                        make_recruit = rec_biomass >= min_recruit_size(ipft)
                        !------------------------------------------------------------------!


                        !------------------------------------------------------------------!
                        !     Print additional information for debugging.                  !
                        !------------------------------------------------------------------!
                        if (printout) then
                           open (unit=66,file=fracfile,status='old',position='append'      &
                                ,action='write')
                           write(unit=66,fmt='(6(1x,i6),4(1x,f12.6),1(12x,l1))')           &
                              current_time%year,current_time%month,current_time%date       &
                              ,ipa,ico,ipft,rectest%nplant,rec_biomass                     &
                              ,csite%repro(ipft,ipa),min_recruit_size(ipft),make_recruit
                           close (unit=66,status='keep')
                        end if
                        !------------------------------------------------------------------!


                        !------------------------------------------------------------------!
                        !      Create new recruit in case there is enough biomass.         !
                        !------------------------------------------------------------------!
                        if (make_recruit) then

                           !----- Add new recruit. ----------------------------------------!
                           inew = inew + 1
                           call copy_recruit(rectest,recruit(inew))
                           !---------------------------------------------------------------!


                           !----- Reset the carbon available for reproduction. ------------!
                           csite%repro(ipft,ipa) = 0.0
                           !---------------------------------------------------------------!
                        end if
                        !------------------------------------------------------------------!
                     else
                        !------------------------------------------------------------------!
                        !     This PFT shouldn't exist... at least not on this patch.      !
                        ! Send the seed litter to the soil pools for decomposition.        !
                        !                                                                  !
                        !  ALS: dont send all seeds to litter!  Keep it for harvesting?    !
                        !  MLO: byield now captures seeds that are harvested.              !
                        !-------------------------------------------- ---------------------!
                        csite%fast_grnd_N(ipa) = csite%fast_grnd_N(ipa)                    &
                                               + csite%repro(ipft,ipa) / c2n_recruit(ipft)
                        csite%fast_grnd_C(ipa) = csite%fast_grnd_C(ipa)                    &
                                               + csite%repro(ipft,ipa)
                        csite%repro(ipft,ipa)  = 0.0
                        !------------------------------------------------------------------!
                     end if
                     !---------------------------------------------------------------------!
                  end do pftloop
                  !------------------------------------------------------------------------!


                  !----- Update the number of cohorts with the recently created. ----------!
                  ncohorts_new = cpatch%ncohorts + inew
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     The number of recruits is now known. If there is any recruit, then !
                  ! we allocate the temporary patch vector with the current number plus    !
                  ! the number of recruits.                                                !
                  !------------------------------------------------------------------------!
                  if (ncohorts_new > cpatch%ncohorts) then
                     nullify(temppatch)
                     allocate(temppatch)
                     call allocate_patchtype(temppatch,cpatch%ncohorts)

                     !----- Fill the temp space with the current patches. -----------------!
                     call copy_patchtype(cpatch,temppatch,1                                &
                                        ,cpatch%ncohorts,1,cpatch%ncohorts)
                     !---------------------------------------------------------------------!

                     !----- Deallocate the current patch. ---------------------------------!
                     call deallocate_patchtype(cpatch)
                     !---------------------------------------------------------------------!

                     !----- Reallocate the current site. ----------------------------------!
                     call allocate_patchtype(cpatch,ncohorts_new)
                     !---------------------------------------------------------------------!

                     !----- Transfer the temp values back in. -----------------------------!
                     call copy_patchtype(temppatch,cpatch,1,temppatch%ncohorts             &
                                        ,1,temppatch%ncohorts)
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !     Add recruits.                                                   !
                     !---------------------------------------------------------------------!
                     inew = 0
                     recloop: do ico = temppatch%ncohorts+1,ncohorts_new
                        !------------------------------------------------------------------!
                        !     Add the recruits, copying the information from the recruit-  !
                        ! ment table, and derive other variables or assume standard        !
                        ! initial values.                                                  !
                        !------------------------------------------------------------------!
                        inew = inew + 1

                        !----- Copy from recruitment table (I). ---------------------------!
                        cpatch%pft(ico)       = recruit(inew)%pft
                        cpatch%hite(ico)      = recruit(inew)%hite
                        cpatch%dbh(ico)       = recruit(inew)%dbh
                        !------------------------------------------------------------------!

                        !----- Carry out standard initialization. -------------------------!
                        call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi),nzg             &
                                                ,cpoly%ntext_soil(:,isi))
                        !------------------------------------------------------------------!


                        !----- Copy from recruitment table (II). --------------------------!
                        cpatch%nplant             (ico) = recruit(inew)%nplant 
                        cpatch%bdeada             (ico) = recruit(inew)%bdeada
                        cpatch%bdeadb             (ico) = recruit(inew)%bdeadb
                        cpatch%paw_avg            (ico) = recruit(inew)%paw_avg
                        cpatch%elongf             (ico) = recruit(inew)%elongf
                        cpatch%phenology_status   (ico) = recruit(inew)%phenology_status
                        cpatch%bleaf              (ico) = recruit(inew)%bleaf
                        cpatch%broot              (ico) = recruit(inew)%broot
                        cpatch%bsapwooda          (ico) = recruit(inew)%bsapwooda
                        cpatch%bsapwoodb          (ico) = recruit(inew)%bsapwoodb
                        cpatch%bbarka             (ico) = recruit(inew)%bbarka
                        cpatch%bbarkb             (ico) = recruit(inew)%bbarkb
                        cpatch%balive             (ico) = recruit(inew)%balive
                        cpatch%bstorage           (ico) = recruit(inew)%bstorage
                        cpatch%leaf_temp          (ico) = recruit(inew)%leaf_temp
                        cpatch%wood_temp          (ico) = recruit(inew)%wood_temp
                        cpatch%leaf_temp_pv       (ico) = recruit(inew)%leaf_temp_pv
                        cpatch%wood_temp_pv       (ico) = recruit(inew)%wood_temp_pv
                        cpatch%leaf_vpdef         (ico) = recruit(inew)%leaf_vpdef
                        do imon=1,13
                           cpatch%cb         (imon,ico) = recruit(inew)%cb         (imon)
                           cpatch%cb_lightmax(imon,ico) = recruit(inew)%cb_lightmax(imon)
                           cpatch%cb_moistmax(imon,ico) = recruit(inew)%cb_moistmax(imon)
                           cpatch%cb_mlmax   (imon,ico) = recruit(inew)%cb_mlmax   (imon)
                        end do
                        cpatch%cbr_bar            (ico) = recruit(inew)%cbr_bar
                        !------------------------------------------------------------------!


                        !----- Initialise the next variables with zeroes... ---------------!
                        cpatch%leaf_water(ico) = 0.0
                        cpatch%wood_water(ico) = 0.0
                        !------------------------------------------------------------------!


                        !------------------------------------------------------------------!
                        !      Because internal water is likely not be zero, we must       !
                        ! ensure that leaf and wood liquid water fractions are consistent  !
                        ! with temperature.                                                !
                        !------------------------------------------------------------------!
                        if ( cpatch%leaf_temp(ico) == t3ple) then
                           cpatch%leaf_fliq (ico) = 0.5
                        elseif ( cpatch%leaf_temp(ico) > t3ple) then
                           cpatch%leaf_fliq (ico) = 1.0
                        else
                           cpatch%leaf_fliq (ico) = 0.0
                        end if
                        if ( cpatch%wood_temp(ico) == t3ple) then
                           cpatch%wood_fliq (ico) = 0.5
                        elseif ( cpatch%wood_temp(ico) > t3ple) then
                           cpatch%wood_fliq (ico) = 1.0
                        else
                           cpatch%wood_fliq (ico) = 0.0
                        end if
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !    Compute initial AGB and Basal Area.  Their derivatives will   !
                        ! be zero.                                                         !
                        !------------------------------------------------------------------!
                        cpatch%agb      (ico) = ed_biomass(cpatch, ico)
                        cpatch%btimber  (ico) = size2bt   ( cpatch%dbh       (ico)         &
                                                          , cpatch%hite      (ico)         &
                                                          , cpatch%bdeada    (ico)         &
                                                          , cpatch%bsapwooda (ico)         &
                                                          , cpatch%bbarka    (ico)         &
                                                          , cpatch%pft       (ico) )
                        cpatch%thbark   (ico) = size2xb   ( cpatch%dbh       (ico)         &
                                                          , cpatch%hite      (ico)         &
                                                          , cpatch%bbarka    (ico)         &
                                                          , cpatch%bbarkb    (ico)         &
                                                          , cpatch%pft       (ico) )
                        cpatch%basarea  (ico) = pio4 * cpatch%dbh(ico)  * cpatch%dbh(ico)
                        cpatch%dagb_dt  (ico) = 0.0
                        cpatch%dlnagb_dt(ico) = 0.0
                        cpatch%dba_dt   (ico) = 0.0
                        cpatch%dlnba_dt (ico) = 0.0
                        cpatch%ddbh_dt  (ico) = 0.0
                        cpatch%dlndbh_dt(ico) = 0.0
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Set new_recruit_flag to 1 indicates that this cohort is      !
                        ! included when we tally agb_recruit, basal_area_recruit.          !
                        !------------------------------------------------------------------!
                        cpatch%new_recruit_flag(ico) = 1
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Update plastic traits (SLA, Vm0).  This must be done before  !
                        ! calculating LAI.                                                 !
                        !------------------------------------------------------------------!
                        select case (trait_plasticity_scheme)
                        case (0) 
                           !----- Trait plasticity is disabled, do nothing. ---------------!
                           continue
                           !---------------------------------------------------------------!
                        case (-2,2) ! Update trait every month
                           !----- Allow recruits to start adapted to their environment. ---!
                           call update_cohort_plastic_trait(cpatch,ico)
                           !---------------------------------------------------------------!
                        end select
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !    Obtain derived properties.                                    !
                        !------------------------------------------------------------------!
                        !----- Find LAI, WAI, and CAI. ------------------------------------!
                        call area_indices(cpatch, ico)
                        !----- Find heat capacity and vegetation internal energy. ---------!
                        call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdeada(ico)            &
                                          ,cpatch%bsapwooda(ico),cpatch%bbarka(ico)        &
                                          ,cpatch%nplant(ico),cpatch%pft(ico)              &
                                          ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))
                        !----- Find total internal water content. -------------------------!
                        call rwc2tw(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico)              &
                                   ,cpatch%bleaf(ico),cpatch%bsapwooda(ico)                &
                                   ,cpatch%bsapwoodb(ico),cpatch%bdeada(ico)               &
                                   ,cpatch%bdeadb(ico),cpatch%broot(ico),cpatch%dbh(ico)   &
                                   ,cpatch%pft(ico),cpatch%leaf_water_int(ico)             &
                                   ,cpatch%wood_water_int(ico))
                        call twi2twe(cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico) &
                                    ,cpatch%nplant(ico),cpatch%leaf_water_im2(ico)         &
                                    ,cpatch%wood_water_im2(ico))
                        !------------------------------------------------------------------!
                        !    Set internal energy to 0., then update it so the effect of    !
                        ! recruits to the carbon balance is accounted for.                 !
                        !------------------------------------------------------------------!
                        cpatch%leaf_energy(ico) = 0.0
                        cpatch%wood_energy(ico) = 0.0
                        call update_veg_energy_cweh(csite,ipa,ico,0.,0.,0.,0.,.false.)
                        !----- Update flags for the biophysical integrator. ---------------!
                        call is_resolvable(csite,ipa,ico)
                        !------------------------------------------------------------------!

                        !----- Update number of cohorts in this site. ---------------------!
                        csite%cohort_count(ipa) = csite%cohort_count(ipa) + 1
                        !------------------------------------------------------------------!
                     end do recloop
                     !---------------------------------------------------------------------!



                     !---- Remove the temporary patch. ------------------------------------!
                     call deallocate_patchtype(temppatch)
                     deallocate(temppatch)
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!
               end do patchloop
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !   Now that recruitment has occured, fuse, terminate, split, and re-sort.  !
               !---------------------------------------------------------------------------!
               update_patch_loop: do ipa = 1,csite%npatches
                  cpatch => csite%patch(ipa)

                  !----- Update the cohort distribution. ----------------------------------!
                  if (cpatch%ncohorts > 0 .and. maxcohort >= 0) then
                     select case (ifusion)
                     case (0)
                        call old_fuse_cohorts(csite,ipa,cpoly%lsl(isi),.false.)
                     case (1)
                        call new_fuse_cohorts(csite,ipa,cpoly%lsl(isi),.false.)
                     end select
                     call terminate_cohorts(csite,ipa,cmet,elim_nplant,elim_lai)
                     call split_cohorts(csite,ipa,cpoly%green_leaf_factor(:,isi))
                  end if
                  !------------------------------------------------------------------------!


                  !----- Sort cohorts by height. ------------------------------------------!
                  call sort_cohorts(cpatch)
                  !------------------------------------------------------------------------!

                  !----- Update the number of cohorts (this is redundant...). -------------!
                  csite%cohort_count(ipa) = cpatch%ncohorts
                  !------------------------------------------------------------------------!


                  !----- Since cohorts may have changed, update patch properties... -------!
                  call update_patch_derived_props(csite,ipa,.true.)
                  !------------------------------------------------------------------------!
               end do update_patch_loop
               !---------------------------------------------------------------------------!

               !----- Since patch properties may have changed, update site properties... --!
               call update_site_derived_props(cpoly,0,isi)
               !---------------------------------------------------------------------------!


               !----- Reset minimum monthly temperature. ----------------------------------!
               cpoly%min_monthly_temp(isi) = huge(1.)
               !---------------------------------------------------------------------------!
            end do siteloop
            !------------------------------------------------------------------------------!
         end do polyloop
         !---------------------------------------------------------------------------------!



      case (1)
         !---------------------------------------------------------------------------------!
         !                                    'big leaf' ED                                !
         !  Growth and reproduction are done together as there is no vertical structure in !
         ! big leaf ED so the cohorts cannot grow vertically.  Therefore daily NPP is      !
         ! is accumulated and monthly the nplant of each cohort is increased (1 cohort per !
         ! patch and 1 patch per pft and disturbance type).
         !---------------------------------------------------------------------------------!


         !----- The big loops start here. -------------------------------------------------!
         polyloop_big: do ipy = 1,cgrid%npolygons
            cpoly => cgrid%polygon(ipy)
            late_spring = (cgrid%lat(ipy) >= 0.0 .and. month == 6) .or.                    &
                          (cgrid%lat(ipy) < 0.0 .and. month == 12)

            !------- Update the repro arrays. ---------------------------------------------!
            call seed_dispersal(cpoly,late_spring)
            !------------------------------------------------------------------------------!

            siteloop_big: do isi = 1,cpoly%nsites
               csite => cpoly%site(isi)

               patchloop_big: do ipa = 1,csite%npatches
                  cpatch => csite%patch(ipa)

                  !------------------------------------------------------------------------!
                  !    There should only be ONE cohort... if there are more, crash         !
                  !------------------------------------------------------------------------!
                  if (cpatch%ncohorts > 1) then
                    write (unit=*,fmt='(a,1x,es12.5)') ' + PATCH   : ',ipa
                    write (unit=*,fmt='(a,1x,es12.5)') ' + NCOHORTS: ',cpatch%ncohorts
                    call fatal_error('NCOHORTS must be 1 for big-leaf runs'                &
                                     ,'reproduction' ,'reproduction.f90')
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      "Loop" over cohorts.  Reproduction does not create new patches,   !
                  ! instead it will add population to the existing cohort.                 !
                  !------------------------------------------------------------------------!
                  cohortloop_big: do ico = 1, cpatch%ncohorts


                     !------ Current PFT. -------------------------------------------------!
                     ipft = cpatch%pft(ico)
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     Check whether to include this PFT or not.  The decision depends !
                     ! on the following decisions.                                         !
                     ! 1.  In case of agriculture or forest plantation, it must also be    !
                     ! allowed in such patches.                                            !
                     ! 2.  The temperature is not limiting reproduction.                   !
                     ! 3.  The user wants reproduction to occur.                           !
                     !---------------------------------------------------------------------!
                     select case (csite%dist_type(ipa))
                     case (1)
                        !----- Agriculture (cropland or pasture). -------------------------!
                        allow_pft =                                                        &
                           include_pft_ag(ipft)                                      .and. &
                           cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and. &
                           repro_scheme                /= 0
                        !------------------------------------------------------------------!

                     case (2)
                        !----- Forest plantation. -----------------------------------------!
                        allow_pft =                                                        &
                           include_pft_fp(ipft)                                      .and. &
                           cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and. &
                           repro_scheme                /= 0
                        !------------------------------------------------------------------!

                     case default
                        !----- Primary or secondary vegetation. ---------------------------!
                        allow_pft =                                                        &
                           cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and. &
                           repro_scheme                /= 0
                        !------------------------------------------------------------------!
                     end select
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !     Update cohort properties in case reproduction is allowed.       !
                     !---------------------------------------------------------------------!
                     if (allow_pft) then

                        !------------------------------------------------------------------!
                        !    Plants don't have size distribution, so use the standard      !
                        ! value of one plant to find a population increase that is         !
                        ! consistent with the expected average biomass.                    !
                        !------------------------------------------------------------------!
                        nplant_inc         = csite%repro(ipft,ipa) / one_plant_c(ipft)
                        cpatch%nplant(ico) = cpatch%nplant(ico) + nplant_inc
                        !------------------------------------------------------------------!


                        !----- Reset the carbon available for reproduction. ---------------!
                        csite%repro(ipft,ipa) = 0.0
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !    Will only reproduce/grow if on-allometry so dont' have to     !
                        ! worry about elongation factor.                                   !
                        !------------------------------------------------------------------!
                        bleaf_plant     = size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft) 
                        broot_plant     = bleaf_plant * q    (ipft)
                        bsapwood_plant  = bleaf_plant * qsw  (ipft) * cpatch%hite(ico)
                        bbark_plant     = bleaf_plant * qbark(ipft) * cpatch%hite(ico)
                        balive_plant    = bleaf_plant + broot_plant + bsapwood_plant       &
                                        + bbark_plant
                        bdead_plant     = size2bd(cpatch%dbh(ico),cpatch%hite(ico),ipft)
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !      Update the productivity terms.                              !
                        !------------------------------------------------------------------!
                        cpatch%today_nppleaf(ico)   =  nplant_inc * bleaf_plant
                        cpatch%today_nppfroot(ico)  =  nplant_inc * broot_plant
                        cpatch%today_nppsapwood(ico)=  nplant_inc * bsapwood_plant
                        cpatch%today_nppbark(ico)   =  nplant_inc * bbark_plant
                        cpatch%today_nppwood(ico)   = agf_bs(ipft) * nplant_inc            &
                                                    * bdead_plant
                        cpatch%today_nppcroot(ico)  = (1. - agf_bs(ipft)) * nplant_inc     &
                                                    * bdead_plant
                        !------------------------------------------------------------------!


                        !------------------------------------------------------------------!
                        !     Update the derived properties since the population may have  !
                        ! changed.                                                         !
                        !------------------------------------------------------------------!
                        !----- Find LAI, WAI, and CAI. ------------------------------------!
                        call area_indices(cpatch, ico)
                        !----- Find heat capacity and vegetation internal energy. ---------!
                        call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdeada(ico)            &
                                          ,cpatch%bsapwooda(ico),cpatch%bbarka(ico)        &
                                          ,cpatch%nplant(ico),cpatch%pft(ico)              &
                                          ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))
                        !----- Find total internal water content. -------------------------!
                        call rwc2tw(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico)              &
                                   ,cpatch%bleaf(ico),cpatch%bsapwooda(ico)                &
                                   ,cpatch%bsapwoodb(ico),cpatch%bdeada(ico)               &
                                   ,cpatch%bdeadb(ico),cpatch%broot(ico),cpatch%dbh(ico)   &
                                   ,cpatch%pft(ico),cpatch%leaf_water_int(ico)             &
                                   ,cpatch%wood_water_int(ico))
                        call twi2twe(cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico) &
                                    ,cpatch%nplant(ico),cpatch%leaf_water_im2(ico)         &
                                    ,cpatch%wood_water_im2(ico))
                        !------------------------------------------------------------------!
                        !    Set internal energy to 0., then update it so the effect of    !
                        ! recruits to the carbon balance is accounted for.                 !
                        !------------------------------------------------------------------!
                        cpatch%leaf_energy(ico) = 0.0
                        cpatch%wood_energy(ico) = 0.0
                        call update_veg_energy_cweh(csite,ipa,ico,0.,0.,0.,0.,.false.)
                        !----- Update flags for the biophysical integrator. ---------------!
                        call is_resolvable(csite,ipa,ico)
                        !------------------------------------------------------------------!
                     else
                        !------------------------------------------------------------------!
                        !     This PFT shouldn't exist... at least not on this patch.      !
                        ! Send the seed litter to the soil pools for decomposition.        !
                        !                                                                  !
                        !  ALS: dont send all seeds to litter!  Keep it for harvesting?    !
                        !  MLO: not here because these are the PFTs that are not allowed   !
                        !       on this patch (for example, trees on agriculture patch).   !
                        !       Perhaps when we convert bseeds to repro, we should ask     !
                        !       how much of the seed pool should be lost to harvest.       !
                        !------------------------------------------------------------------!
                        csite%fast_grnd_N(ipa) = csite%fast_grnd_N(ipa)                    &
                                               + csite%repro(ipft,ipa) / c2n_recruit(ipft)
                        csite%fast_grnd_C(ipa) = csite%fast_grnd_C(ipa)                    &
                                               + csite%repro(ipft,ipa)
                        csite%repro(ipft,ipa)  = 0.0
                        !------------------------------------------------------------------!
                     end if
                     !---------------------------------------------------------------------!
                  end do cohortloop_big
                  !------------------------------------------------------------------------!
               end do patchloop_big
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !     Update derived properties.                                            !
               !---------------------------------------------------------------------------!
               update_patch_loop_big: do ipa = 1,csite%npatches
                  cpatch => csite%patch(ipa)

                  !----- Update the number of cohorts (this is redundant...). -------------!
                  csite%cohort_count(ipa) = cpatch%ncohorts
                  !------------------------------------------------------------------------!

                  !----- Since cohorts may have changed, update patch properties... -------!
                  call update_patch_derived_props(csite,ipa,.true.)
                  !------------------------------------------------------------------------!
               end do update_patch_loop_big
               !---------------------------------------------------------------------------!

               !----- Since patch properties may have changed, update site properties... --!
               call update_site_derived_props(cpoly,0,isi)
               !---------------------------------------------------------------------------!

               !----- Reset minimum monthly temperature. ----------------------------------!
               cpoly%min_monthly_temp(isi) = huge(1.)
               !---------------------------------------------------------------------------!
            end do siteloop_big
            !------------------------------------------------------------------------------!
         end do polyloop_big
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end subroutine reproduction_driver
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine will update the repro structure.  The way the dispersal is going !
   ! to happen depends on whether the user allows dispersal between sites or within sites  !
   ! only.  Broadleaf deciduous PFTs do their dispersal only in late spring.  This is be-  !
   ! cause of their unique storage respiration.  Other PFTs disperse seeds at any time of  !
   ! the year.                                                                             !
   !---------------------------------------------------------------------------------------!
   subroutine seed_dispersal(cpoly,late_spring)
      use ed_max_dims        , only : n_pft                 ! ! intent(in)
      use phenology_coms     , only : repro_scheme          ! ! intent(in)
      use ed_state_vars      , only : polygontype           & ! structure
                                    , sitetype              & ! structure
                                    , patchtype             ! ! structure
      use pft_coms           , only : recruittype           & ! structure
                                    , nonlocal_dispersal    & ! intent(in)
                                    , seedling_mortality    & ! intent(in)
                                    , phenology             & ! intent(in)
                                    , plant_min_temp        & ! intent(in)
                                    , seed_rain             & ! intent(in)
                                    , one_plant_c           & ! intent(in)
                                    , include_pft           & ! intent(in)
                                    , include_pft_ag        & ! intent(in)
                                    , include_pft_fp        ! ! intent(in)
      use ed_misc_coms       , only : ibigleaf              & ! intent(in)
                                    , frqsumi               ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(polygontype), target     :: cpoly       ! Current polygon             [     ---]
      logical          , intent(in) :: late_spring ! Is it late spring           [     T|F]
      !----- Local variables. -------------------------------------------------------------!
      type(sitetype)   , pointer    :: csite        ! Current site               [     ---]
      type(sitetype)   , pointer    :: donsite      ! Donor site                 [     ---]
      type(sitetype)   , pointer    :: recsite      ! Receptor site              [     ---]
      type(patchtype)  , pointer    :: donpatch     ! Donor patch                [     ---]
      logical                       :: allow_pft    ! Flag: is this PFT allowed? [     ---]
      integer                       :: isi          ! Site counter               [     ---]
      integer                       :: ipa          ! Patch counter              [     ---]
      integer                       :: ipft         ! PFT counter                [     ---]
      integer                       :: recsi        ! Receptor site counter      [     ---]
      integer                       :: donsi        ! Donor site counter         [     ---]
      integer                       :: recpa        ! Receptor patch counter     [     ---]
      integer                       :: donpa        ! Donor patch counter        [     ---]
      integer                       :: donco        ! Donor cohort counter       [     ---]
      integer                       :: donpft       ! Donor PFT                  [     ---]
      real                          :: bseedling    ! Surviving seedling biomass [  kgC/m2]
      real                          :: bseed_stays  ! Seedl. biomass that stays  [  kgC/m2]
      real                          :: bseed_maygo  ! Seedl. biomass that may go [  kgC/m2]
      real                          :: bseed_xpatch ! Seedl. X-patch exchange    [  kgC/m2]
      real                          :: seed_rain_c  ! Seed rain carbon input     [  kgC/m2]
      real, dimension(n_pft)        :: pft_seedrain ! Sum of all ext. seed rain  [  kgC/m2]
      real, dimension(n_pft)        :: pft_repro0   ! Seed bank before dispersal [  kgC/m2]
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    The following variables are used to ensure reproduction is conserving carbon.   !
      !------------------------------------------------------------------------------------!
      pft_seedrain(:) = 0.0
      pft_repro0  (:) = 0.0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Here we decide how to disperse seeds based on the reproduction scheme.         !
      !------------------------------------------------------------------------------------!
      select case (repro_scheme)
      case (0)
         !------ No reproduction, quit. ---------------------------------------------------!
         return
         !---------------------------------------------------------------------------------!
      case default
         !---------------------------------------------------------------------------------!
         !     Seed rain should apply to all patches and sites, include them to the seed   !
         ! bank.                                                                           !
         !---------------------------------------------------------------------------------!
         siteloop0: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)


            !------------------------------------------------------------------------------!
            !      Loop over the patches and PFTs to include seed_rain.                    !
            !------------------------------------------------------------------------------!
            seedloop: do ipa = 1,csite%npatches
               pftloop: do ipft=1,n_pft
                  !------------------------------------------------------------------------!
                  !     Save original seed bank, for budget check.                         !
                  !------------------------------------------------------------------------!
                  pft_repro0(ipft) = pft_repro0(ipft)                                      &
                                   + csite%repro(ipft,ipa) * csite%area(ipa)               &
                                   * cpoly%area(isi)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Check whether to include this PFT or not.  The decision depends on !
                  ! the following decisions.                                               !
                  ! 1.  The PFT is included in this simulation.   In case of agriculture   !
                  !     or forest plantation, it must also be allowed in such patches.     !
                  ! 2.  The temperature is not limiting reproduction.                      !
                  ! 3.  The user wants reproduction to occur.                              !
                  !------------------------------------------------------------------------!
                  select case (csite%dist_type(ipa))
                  case (1)
                     !----- Agriculture (cropland or pasture). ----------------------------!
                     allow_pft =                                                           &
                        include_pft_ag(ipft)                                      .and.    &
                        cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and.    &
                        repro_scheme                /= 0
                     !---------------------------------------------------------------------!

                  case (2)
                     !----- Forest plantation. --------------------------------------------!
                     allow_pft =                                                           &
                        include_pft_fp(ipft)                                      .and.    &
                        cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and.    &
                        repro_scheme                /= 0
                     !---------------------------------------------------------------------!

                  case default
                     !----- Primary or secondary vegetation. ------------------------------!
                     allow_pft =                                                           &
                        include_pft(ipft)                                         .and.    &
                        cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and.    &
                        repro_scheme                /= 0
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    In case this PFT is allowed, update the seed pool.  Also update the !
                  ! seed rain budget checks.                                               !
                  !------------------------------------------------------------------------!
                  if (allow_pft) then
                     seed_rain_c                 = seed_rain(ipft) * one_plant_c(ipft)
                     csite%repro      (ipft,ipa) = csite%repro(ipft,ipa) + seed_rain_c
                     csite%cbudget_seedrain(ipa) = csite%cbudget_seedrain(ipa)             &
                                                 + seed_rain_c * frqsumi
                     pft_seedrain(ipft)          = pft_seedrain(ipft)                      &
                                                 + seed_rain_c                             &
                                                 * csite%area(ipa) * cpoly%area(isi)
                  end if
                  !------------------------------------------------------------------------!
               end do pftloop
               !---------------------------------------------------------------------------!
            end do seedloop
            !------------------------------------------------------------------------------!
         end do siteloop0
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      In this select block we decide the seed dispersal depending the cross-site    !
      !------------------------------------------------------------------------------------!
      select case (repro_scheme)
      case (0)
         !------ This should never happen, but just in case, quit. ------------------------!
         return
         !---------------------------------------------------------------------------------!
      case (1)
         !---------------------------------------------------------------------------------!
         !     Seeds are dispersed amongst patches that belong to the same site, but they  !
         ! cannot go outside their native site.                                            !
         !---------------------------------------------------------------------------------!
         siteloop1: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)
            !------------------------------------------------------------------------------!
            !      Loop over the donor cohorts.                                            !
            !------------------------------------------------------------------------------!
            donpaloop1: do donpa = 1,csite%npatches
               donpatch => csite%patch(donpa)

               doncoloop1: do donco = 1, donpatch%ncohorts

                  !----- Define an alias for PFT. -----------------------------------------!
                  donpft = donpatch%pft(donco)
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !    Find the biomass of survivor seedlings.  Units: kgC/m2              !
                  !------------------------------------------------------------------------!
                  if (phenology(donpft) /= 2 .or. late_spring) then
                     bseedling   = donpatch%nplant(donco) * donpatch%bseeds(donco)         &
                                 * (1.0 - seedling_mortality(donpft))
                     select case (ibigleaf)
                     case (0)
                        bseed_stays = bseedling * (1.0 - nonlocal_dispersal(donpft))
                        bseed_maygo = bseedling * nonlocal_dispersal(donpft)
                     case (1)
                        !---- if bigleaf cannot disperse seedlings to other patches -------!
                        bseed_stays = bseedling
                        bseed_maygo = 0.
                     end select

                  else
                     !----- Not a good time for reproduction.  No seedlings. --------------!
                     bseedling     = 0.
                     bseed_stays   = 0.
                     bseed_maygo   = 0.
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Subtract the seeds that leave the patch from the seed rain.  Some  !
                  ! of them will land again in this patch, and we correct for this further !
                  ! down.                                                                  !
                  !------------------------------------------------------------------------!
                  csite%cbudget_seedrain(donpa) = csite%cbudget_seedrain(donpa)            &
                                                - bseed_maygo * frqsumi
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     This is the seed biomass that crosses patches, in units of         !
                  ! kgC/m2_recp_patch.                                                     !
                  !                                                                        !
                  !  Variable names for rationale:                                         !
                  !  ------------------------------                                        !
                  ! DPY - Seeds that may leave donor patch     [           kgC/m2_polygon] !
                  ! DPA - Seeds that may leave donor patch     [       kgC/m2_donor_patch] !
                  ! AD  - Area of the donor patch, in          [m2_donor_patch/m2_polygon] !
                  ! RPY - Seeds that land in the recp. patch   [           kgC/m2_polygon] !
                  ! RPA - Seeds that land in the recp. patch   [        kgC/m2_recp_patch] !
                  ! AR  - Area of the receptor patch           [ m2_recp_patch/m2_polygon] !
                  !                                                                        !
                  !  Rationale:                                                            !
                  ! -------------                                                          !
                  !                                                                        !
                  ! (1) DPY = DPA * AD + 0 * (1-AD)   (General patch->polygon definition)  !
                  ! (2) RPY = RPA * AR + 0 * (1-AR)   (General patch->polygon definition)  !
                  ! (3) RPY = DPY * AR                (Probabilty is prop. to area)        !
                  ! (4) RPY = DPA * AD * AR           (1->3)                               !
                  ! (5) RPA = DPA * AD                (4->2, regardless of the patch)      !
                  !------------------------------------------------------------------------!
                  bseed_xpatch = bseed_maygo * csite%area(donpa)
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !   Spread the seedlings across all patches in this site.                !
                  !------------------------------------------------------------------------!
                  recpaloop1: do recpa = 1,csite%npatches

                     !---------------------------------------------------------------------!
                     !     Add the non-local dispersal evenly across all patches,          !
                     ! including the donor patch.  We must scale the biomass by the        !
                     ! combined area of this patch and site so the total carbon is         !
                     ! preserved.                                                          !
                     !---------------------------------------------------------------------!
                     csite%repro(donpft,recpa) = csite%repro(donpft,recpa) + bseed_xpatch
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    Donor and receptor patches are different.  Account for cross-    !
                     ! patch seed exchange in the seed rain pool.  Some of them are        !
                     ! from the patch itself, but they should be added because we          !
                     ! subtracted all the non-local dispersal outside the receptor site    !
                     ! loop.                                                               !
                     !---------------------------------------------------------------------!
                     csite%cbudget_seedrain(recpa) = csite%cbudget_seedrain(recpa)         &
                                                   + bseed_xpatch * frqsumi
                     !---------------------------------------------------------------------!




                     !---------------------------------------------------------------------!
                     !      Include the local dispersal if this is the donor patch.        !
                     !---------------------------------------------------------------------!
                     if (recpa == donpa) then
                        csite%repro(donpft,recpa) = csite%repro(donpft,recpa) + bseed_stays
                     end if
                     !---------------------------------------------------------------------!
                  end do recpaloop1
                  !------------------------------------------------------------------------!
               end do doncoloop1
               !---------------------------------------------------------------------------!
            end do donpaloop1
            !------------------------------------------------------------------------------!
         end do siteloop1
         !---------------------------------------------------------------------------------!

      case default

         !---------------------------------------------------------------------------------!
         !     Seeds are dispersed amongst patches that belong to the same polygon.  They  !
         ! are allowed to go from one site to the other, but they cannot go outside their  !
         ! native polygon.                                                                 !
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !      Loop over the donor cohorts.                                               !
         !---------------------------------------------------------------------------------!
         donsiloop2: do donsi = 1,cpoly%nsites
            donsite => cpoly%site(donsi)

            donpaloop2: do donpa = 1,donsite%npatches
               donpatch => donsite%patch(donpa)

               doncoloop2: do donco = 1, donpatch%ncohorts

                  !----- Define an alias for PFT. -----------------------------------------!
                  donpft = donpatch%pft(donco)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Find the biomass of survivor seedlings.  Units: kgC/m2              !
                  !------------------------------------------------------------------------!
                  if (phenology(donpft) /= 2 .or. late_spring) then
                     bseedling   = donpatch%nplant(donco) * donpatch%bseeds(donco)         &
                                 * (1.0 - seedling_mortality(donpft))

                     select case (ibigleaf)
                     case (0)
                        bseed_stays = bseedling * (1.0 - nonlocal_dispersal(donpft))
                        bseed_maygo = bseedling * nonlocal_dispersal(donpft)
                     case (1)
                        !---- if bigleaf cannot disperse seedlings to other patches -------!
                        bseed_stays = bseedling
                        bseed_maygo = 0.
                     end select
                  else
                     !----- Not a good time for reproduction.  No seedlings. --------------!
                     bseedling   = 0.
                     bseed_stays = 0.
                     bseed_maygo = 0.
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Subtract the seeds that leave the patch from the seed rain.  Some  !
                  ! of them will land again in this patch, and we correct for this further !
                  ! down.                                                                  !
                  !------------------------------------------------------------------------!
                  csite%cbudget_seedrain(donpa) = csite%cbudget_seedrain(donpa)            &
                                                - bseed_maygo * frqsumi
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     This is the seed biomass that crosses patches, in units of         !
                  ! kgC/m2_recp_patch.                                                     !
                  !                                                                        !
                  !  Variable names for rationale:                                         !
                  !  ------------------------------                                        !
                  ! DPY - Seeds that may leave donor patch     [           kgC/m2_polygon] !
                  ! DPA - Seeds that may leave donor patch     [       kgC/m2_donor_patch] !
                  ! AD  - Area of the donor patch, in          [m2_donor_patch/m2_polygon] !
                  ! RPY - Seeds that land in the recp. patch   [           kgC/m2_polygon] !
                  ! RPA - Seeds that land in the recp. patch   [        kgC/m2_recp_patch] !
                  ! AR  - Area of the receptor patch           [ m2_recp_patch/m2_polygon] !
                  !                                                                        !
                  !  Rationale:                                                            !
                  ! -------------                                                          !
                  !                                                                        !
                  ! (1) DPY = DPA * AD + 0 * (1-AD)   (General patch->polygon definition)  !
                  ! (2) RPY = RPA * AR + 0 * (1-AR)   (General patch->polygon definition)  !
                  ! (3) RPY = DPY * AR                (Probabilty is prop. to area)        !
                  ! (4) RPY = DPA * AD * AR           (1->3)                               !
                  ! (5) RPA = DPA * AD                (4->2, regardless of the patch)      !
                  !------------------------------------------------------------------------!
                  bseed_xpatch = bseed_maygo * csite%area(donpa) * cpoly%area(donsi)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !   Spread the seedlings across all patches in this polygon.             !
                  !------------------------------------------------------------------------!
                  recsiloop2: do recsi = 1,cpoly%nsites
                     recsite => cpoly%site(recsi)
                     recpaloop2: do recpa = 1,recsite%npatches

                        !------------------------------------------------------------------!
                        !     Add the non-local dispersal evenly across all patches,       !
                        ! including the donor patch.  We must scale the biomass by the     !
                        ! combined area of this patch and site so the total carbon is      !
                        ! preserved.                                                       !
                        !------------------------------------------------------------------!
                        recsite%repro(donpft,recpa) = recsite%repro(donpft,recpa)          &
                                                    + bseed_xpatch
                        !------------------------------------------------------------------!


                        !------------------------------------------------------------------!
                        !    Donor and receptor patches are different.  Account for cross- !
                        ! patch seed exchange in the seed rain pool.  Some of them are     !
                        ! from the patch itself, but they should be added because we       !
                        ! subtracted all the non-local dispersal outside the receptor site !
                        ! loop.                                                            !
                        !------------------------------------------------------------------!
                        csite%cbudget_seedrain(recpa) = csite%cbudget_seedrain(recpa)      &
                                                      + bseed_xpatch * frqsumi
                        !------------------------------------------------------------------!




                        !------------------------------------------------------------------!
                        !      Include the local dispersal if this is the donor patch.     !
                        !------------------------------------------------------------------!
                        if (recpa == donpa .and. recsi == donsi) then
                           recsite%repro(donpft,recpa) = recsite%repro(donpft,recpa)       &
                                                       + bseed_stays
                        end if
                        !------------------------------------------------------------------!

                     end do recpaloop2
                     !---------------------------------------------------------------------!
                  end do recsiloop2
                  !------------------------------------------------------------------------!
               end do doncoloop2
               !---------------------------------------------------------------------------!
            end do donpaloop2
            !------------------------------------------------------------------------------!
         end do donsiloop2
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Before we leave the seeds grow, we check that the carbon is being conserved.  !
      !------------------------------------------------------------------------------------!
      call check_budget_bseeds(cpoly,pft_repro0,pft_seedrain)
      !------------------------------------------------------------------------------------!

      return
   end subroutine seed_dispersal
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This sub-routine checks that seed biomass and seed bank match.                   !
   !---------------------------------------------------------------------------------------!
   subroutine check_budget_bseeds(cpoly,pft_repro0,pft_seedrain)
      use ed_state_vars      , only : polygontype           & ! structure
                                    , sitetype              & ! structure
                                    , patchtype             ! ! structure
      use ed_max_dims        , only : n_pft                 ! ! intent(in)
      use pft_coms           , only : seedling_mortality    & ! intent(in)
                                    , min_recruit_size      ! ! intent(in)
      use consts_coms        , only : r_tol_trunc           ! ! intent(in)
      use ed_misc_coms       , only : current_time          ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(polygontype)        , target     :: cpoly        ! Current polygon     [    ---]
      real   , dimension(n_pft), intent(in) :: pft_repro0   ! Previous seed bank  [ kgC/m2]
      real   , dimension(n_pft), intent(in) :: pft_seedrain ! Total seed rain     [ kgC/m2]
      !----- Local variables. -------------------------------------------------------------!
      type(sitetype)           , pointer    :: csite        ! Current site        [    ---]
      type(patchtype)          , pointer    :: cpatch       ! Current patch       [    ---]
      integer                               :: isi          ! Site counter        [    ---]
      integer                               :: ipa          ! Patch counter       [    ---]
      integer                               :: ico          ! Cohort counter      [    ---]
      integer                               :: ipft         ! PFT counter         [    ---]
      real   , dimension(n_pft)             :: toler_bseeds ! Scale for tolerance [ kgC/m2]
      real   , dimension(n_pft)             :: pft_bseeds   ! Total viable bseeds [ kgC/m2]
      real   , dimension(n_pft)             :: pft_repro    ! Total seed bank     [ kgC/m2]
      real   , dimension(n_pft)             :: resid_bseeds ! Residual seed stock [ kgC/m2]
      logical, dimension(n_pft)             :: ok_seedbank  ! Consistent seedbank [ kgC/m2]
      !----- Local constants. -------------------------------------------------------------!
      character(len=11)        , parameter  :: fmth='(a,7(1x,a))'
      character(len=23)        , parameter  :: fmtp='(i5,6(1x,es12.5),1x,l1)'
      character(len=27)        , parameter  :: fmtt='(a,i4.4,2(1x,i2.2),1x,f6.0)'
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Initialise seed checking pools.                                                 !
      !------------------------------------------------------------------------------------!
      pft_repro (:) = 0.0
      pft_bseeds(:) = 0.0
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Loop through all the sites, patches and cohorts, and add the total seed        !
      ! biomass.  
      !------------------------------------------------------------------------------------!
      siteloop: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         !----- Patch loop. ---------------------------------------------------------------!
         patchloop: do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)

            !----- First loop: PFT loop, in which we add seeds in the seed bank. ----------!
            pftreproloop: do ipft = 1, n_pft
               pft_repro(ipft) = pft_repro(ipft)                                           &
                               + csite%repro(ipft,ipa) * csite%area(ipa) * cpoly%area(isi)
            end do pftreproloop
            !------------------------------------------------------------------------------!

            !----- Second loop: Add seeds that went to the seed bank. ---------------------!
            cohortloop: do ico = 1,cpatch%ncohorts
               ipft = cpatch%pft(ico)
               pft_bseeds(ipft) = pft_bseeds(ipft)                                         &
                                + cpatch%nplant(ico) * (1.0 - seedling_mortality(ipft))    &
                                * cpatch%bseeds(ico) * csite%area(ipa) * cpoly%area(isi)
            end do cohortloop
            !------------------------------------------------------------------------------!
         end do patchloop
         !---------------------------------------------------------------------------------!
      end do siteloop
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     In this loop, we check that the total seed bank is consistent for each PFT.    !
      !------------------------------------------------------------------------------------!
      pftcheckloop: do ipft=1,n_pft
         resid_bseeds(ipft) = pft_repro (ipft) - pft_bseeds(ipft) - pft_seedrain(ipft)     &
                            - pft_repro0(ipft)
         toler_bseeds(ipft) = r_tol_trunc * max(pft_repro(ipft),min_recruit_size(ipft))
         ok_seedbank (ipft) = abs(resid_bseeds(ipft)) < toler_bseeds(ipft)
      end do pftcheckloop
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     In case anything is inconsitent, we print the information on screen and stop   !
      ! the simulation.                                                                    !
      !------------------------------------------------------------------------------------!
      if (.not. all(ok_seedbank(:))) then
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|====================================================|'
         write(unit=*,fmt='(a)')  '|          !!!   Bseeds budget failed   !!!          |'
         write(unit=*,fmt='(a)')  '|----------------------------------------------------|'
         write(unit=*,fmt=fmtt )  ' TIME                : ',current_time%year              &
                                                           ,current_time%month             &
                                                           ,current_time%date              &
                                                           ,current_time%time
         write(unit=*,fmt='(a)')  '|----------------------------------------------------|'
         write(unit=*,fmt=fmth )  ' IPFT','REPRO_BEFORE','   REPRO_NOW','      BSEEDS'     &
                                         ,'   SEED_RAIN','       RESID','       TOLER'     &
                                         ,'  ACCEPTABLE'
         do ipft=1,n_pft
            write(unit=*,fmt=fmtp ) ipft,pft_repro0(ipft),pft_repro(ipft),pft_bseeds(ipft) &
                                        ,pft_seedrain(ipft),resid_bseeds(ipft)             &
                                        ,toler_bseeds(ipft),ok_seedbank(ipft)
         end do
         write(unit=*,fmt='(a)')  '|----------------------------------------------------|'
         write(unit=*,fmt='(a)')  ' '
         call fatal_error('Budget check has failed, see message above.'                    &
                         ,'check_budget_bseeds','reproduction.f90')
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine check_budget_bseeds
   !=======================================================================================!
   !=======================================================================================!
end module reproduction
!==========================================================================================!
!==========================================================================================!
