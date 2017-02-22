!==========================================================================================!
!==========================================================================================!
!      This subroutine will drive the reproduction, based on its carbon availability and   !
! PFT-specific reproduction properties.  No reproduction will happen if the user didn't    !
! want it, in which case the seedling biomass will go to the litter pools.                 !
!------------------------------------------------------------------------------------------!
subroutine reproduction(cgrid, month)
   use ed_state_vars      , only : edtype                   & ! structure
                                 , polygontype              & ! structure
                                 , sitetype                 & ! structure
                                 , patchtype                & ! structure
                                 , allocate_patchtype       & ! subroutine
                                 , copy_patchtype           & ! subroutine
                                 , deallocate_patchtype     ! ! subroutine
   use pft_coms           , only : recruittype              & ! structure
                                 , zero_recruit             & ! subroutine
                                 , copy_recruit             & ! subroutine
                                 , seedling_mortality       & ! intent(in)
                                 , c2n_stem                 & ! intent(in)
                                 , l2n_stem                 & ! intent(in)
                                 , min_recruit_size         & ! intent(in)
                                 , one_plant_c              & ! intent(in)
                                 , c2n_recruit              & ! intent(in)
                                 , seed_rain                & ! intent(in)
                                 , include_pft              & ! intent(in)
                                 , include_pft_ag           & ! intent(in)
                                 , include_pft_fp           & ! intent(in)
                                 , qsw                      & ! intent(in)
                                 , q                        & ! intent(in)
                                 , agf_bs                   & ! intent(in)
                                 , sla                      & ! intent(in)
                                 , hgt_min                  & ! intent(in)
                                 , plant_min_temp           ! ! intent(in)
   use decomp_coms        , only : f_labile                 ! ! intent(in)
   use ed_max_dims        , only : n_pft                    ! ! intent(in)
   use fuse_fiss_utils    , only : sort_cohorts             & ! subroutine
                                 , terminate_cohorts        & ! subroutine
                                 , fuse_cohorts             & ! subroutine
                                 , split_cohorts            & ! subroutine
                                 , rescale_patches          ! ! subroutine
   use phenology_coms     , only : repro_scheme             ! ! intent(in)
   use mem_polygons       , only : maxcohort                ! ! intent(in)
   use consts_coms        , only : pio4                     ! ! intent(in)
   use ed_therm_lib       , only : calc_veg_hcap            ! ! function
   use allometry          , only : dbh2bd                   & ! function
                                 , size2bl                  & ! function
                                 , h2dbh                    & ! function
                                 , ed_biomass               & ! function
                                 , area_indices             & ! subroutine
                                 , dbh2krdepth              ! ! function
   use grid_coms          , only : nzg                      ! ! intent(in)
   use ed_misc_coms       , only : ibigleaf                 ! ! intent(in)
   use phenology_aux      , only : pheninit_balive_bstorage ! ! intent(in)
   use budget_utils       , only : update_budget            ! ! sub-routine
   use therm_lib          , only : cmtl2uext                ! ! function
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)     , target     :: cgrid
   integer          , intent(in) :: month
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype), pointer          :: cpoly
   type(sitetype)   , pointer          :: csite
   type(patchtype)  , pointer          :: cpatch
   type(patchtype)  , pointer          :: temppatch
   type(recruittype), dimension(n_pft) :: recruit
   type(recruittype)                   :: rectest
   integer                             :: ipy
   integer                             :: isi
   integer                             :: ipa
   integer                             :: ico
   integer                             :: ipft
   integer                             :: inew
   integer                             :: ncohorts_new
   logical                             :: late_spring
   logical                             :: allow_pft
   real                                :: elim_nplant
   real                                :: elim_lai
   real                                :: nplant_inc
   real                                :: bleaf_plant
   real                                :: bdead_plant
   real                                :: broot_plant
   real                                :: bsapwood_plant
   real                                :: balive_plant
   real                                :: rec_biomass
   !----- Saved variables -----------------------------------------------------------------!
   logical          , save             :: first_time = .true.
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    If this is the first time, check whether the user wants reproduction.  If not,     !
   ! kill all potential recruits and send their biomass to the litter pool.                !
   !---------------------------------------------------------------------------------------!
   if (first_time) then
      if (repro_scheme == 0) seedling_mortality(1:n_pft) = 1.0
      first_time = .false.
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Decide which vegetation structure to use.                                         !
   !---------------------------------------------------------------------------------------!
   select case (ibigleaf)
   case (0)
      !------------------------------------------------------------------------------------!
      !                        Size- and age_structure reproduction                        !
      !------------------------------------------------------------------------------------!


      !----- The big loops start here. -------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
      
         !---------------------------------------------------------------------------------!
         !     Check whether this is late spring/early summer.  This is needed for         !
         ! temperate broadleaf deciduous trees.  Late spring means June in the Northern    !
         ! Hemisphere, or December in the Southern Hemisphere.                             !
         !---------------------------------------------------------------------------------!
         late_spring = (cgrid%lat(ipy) >= 0.0 .and. month == 6) .or.                       &
                       (cgrid%lat(ipy) < 0.0 .and. month == 12)

         cpoly => cgrid%polygon(ipy)
         siteloop_sort: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !------------------------------------------------------------------------------!
            !      Cohorts may have grown differently, so we need to sort them by size.    !
            !------------------------------------------------------------------------------!
            patchloop_sort: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)
               call sort_cohorts(cpatch)
            end do patchloop_sort
         end do siteloop_sort

         !------- Update the repro arrays. ------------------------------------------------!
         call seed_dispersal(cpoly,late_spring)
         !---------------------------------------------------------------------------------!

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !------------------------------------------------------------------------------!
            !    For the recruitment to happen, four requirements must be met:             !
            !    1.  PFT is included in this simulation;                                   !
            !    2.  It is not too cold (min_monthly_temp > plant_min_temp - 5)            !
            !    3.  We are dealing with EITHER a non-agriculture patch OR                 !
            !        a PFT that could exist in an agricultural patch.                      !
            !    4.  There must be sufficient carbon to form the recruits.                 !
            !------------------------------------------------------------------------------!
            patchloop: do ipa = 1,csite%npatches
               inew = 0
               call zero_recruit(n_pft,recruit)
               cpatch => csite%patch(ipa)

               !---- This time we loop over PFTs, not cohorts. ----------------------------!
               pftloop: do ipft = 1, n_pft

                  !------------------------------------------------------------------------!
                  !     Check whether to include this PFT or not.  The decision depends    !
                  ! on the following decisions.                                            !
                  ! 1.  The PFT is included in this simulation.   In case of agriculture   !
                  !     or forest plantation, it must also be allowed in such patches.     !
                  ! 2.  The temperature is not limiting reproduction.                      !
                  ! 3.  The user wants reproduction to occur.                              !
                  !------------------------------------------------------------------------!
                  select case (csite%dist_type(ipa))
                  case (1)
                     !----- Agriculture (cropland or pasture). ----------------------------!
                     allow_pft =                                                           &
                           include_pft_ag(ipft)                                      .and. &
                           cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and. &
                           repro_scheme                /= 0
                     !---------------------------------------------------------------------!

                  case (2)
                     !----- Forest plantation. --------------------------------------------!
                     allow_pft =                                                           &
                           include_pft_fp(ipft)                                      .and. &
                           cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and. &
                           repro_scheme                /= 0
                     !---------------------------------------------------------------------!

                  case default
                     !----- Primary or secondary vegetation. ------------------------------!
                     allow_pft =                                                           &
                           include_pft(ipft)                                         .and. &
                           cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and. &
                           repro_scheme                /= 0
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     If this PFT is allowed, check for recruitment.                     !
                  !------------------------------------------------------------------------!
                  if (allow_pft) then
                     !---------------------------------------------------------------------!
                     !       Initialise recruit, and assume that the new cohorts are in    !
                     ! thermal equilibrium with the canopy air space.                      !
                     !---------------------------------------------------------------------!
                     rectest%pft          = ipft
                     rectest%leaf_temp    = csite%can_temp (ipa)
                     rectest%wood_temp    = csite%can_temp (ipa)
                     rectest%leaf_temp_pv = csite%can_temp (ipa)
                     rectest%wood_temp_pv = csite%can_temp (ipa)
                     rectest%leaf_vpdef   = csite%can_vpdef(ipa)
                     !---------------------------------------------------------------------!

                     !---------------------------------------------------------------------!
                     !    Recruits start at minimum height and dbh and bleaf are           !
                     ! calculated from that.                                               !
                     !---------------------------------------------------------------------!
                     rectest%hite      = hgt_min(ipft)
                     rectest%dbh       = h2dbh(rectest%hite, ipft)
                     rectest%krdepth   = dbh2krdepth(rectest%hite,rectest%dbh              &
                                                    ,rectest%pft,cpoly%lsl(isi))
                     rectest%bdead     = dbh2bd(rectest%dbh, ipft)

                     call pheninit_balive_bstorage(nzg,rectest%pft,rectest%krdepth         &
                                                  ,rectest%hite,rectest%dbh                &
                                                  ,csite%soil_water(:,ipa)                 &
                                                  ,cpoly%ntext_soil(:,isi)                 &
                                                  ,cpoly%green_leaf_factor(:,isi)          &
                                                  ,rectest%paw_avg,rectest%elongf          &
                                                  ,rectest%phenology_status                &
                                                  ,rectest%bleaf,rectest%broot             &
                                                  ,rectest%bsapwooda,rectest%bsapwoodb     &
                                                  ,rectest%balive,rectest%bstorage)

                     !---------------------------------------------------------------------!
                     !     Find the expected population from the reproduction stocks, and  !
                     ! also include the "seed_rain" term, which represents sources of      !
                     ! seed that are external to the polygon.                              !
                     !---------------------------------------------------------------------!
                     rectest%nplant    = csite%repro(ipft,ipa)                             &
                                       / ( rectest%balive + rectest%bdead                  &
                                         + rectest%bstorage)                               &
                                       + seed_rain(ipft)
                     !---------------------------------------------------------------------!



                     !----- Find the total biomass for this potential new cohort. ---------!
                     rec_biomass = rectest%nplant * ( rectest%balive + rectest%bdead       &
                                                    + rectest%bstorage )
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !      Create new recruit in case there is enough biomass.  In case   !
                     ! repro_scheme is 3, plants do not reproduce when they are not on     !
                     ! allometry.                                                          !
                     !---------------------------------------------------------------------!
                     select case (repro_scheme)
                     case (3)
                        if ( rec_biomass              >  min_recruit_size(ipft) .and.      &
                             rectest%phenology_status == 0                      ) then

                           !----- Add new recruit. ----------------------------------------!
                           inew = inew + 1
                           call copy_recruit(rectest,recruit(inew))
                           !---------------------------------------------------------------!

                           !----- Reset the carbon available for reproduction. ------------!
                           csite%repro(ipft,ipa) = 0.0                          
                           !---------------------------------------------------------------!
                        end if
                        !------------------------------------------------------------------!

                     case default
                        if (rec_biomass > min_recruit_size(ipft)) then

                           !----- Add new recruit. ----------------------------------------!
                           inew = inew + 1
                           call copy_recruit(rectest,recruit(inew))
                           !---------------------------------------------------------------!


                           !----- Reset the carbon available for reproduction. ------------!
                           csite%repro(ipft,ipa) = 0.0                          
                           !---------------------------------------------------------------!
                        end if
                     end select
                     !---------------------------------------------------------------------!
                  else
                     !---------------------------------------------------------------------!
                     !     This PFT shouldn't exist... at least not on this patch.   Send  !
                     ! the seed litter to the soil pools for decomposition.                !
                     !                                                                     !
                     !     ALS: dont send all seeds to litter!  Keep it for harvesting?    !
                     !     MLO: not here because these are the PFTs that are not allowed   !
                     !          on this patch (for example, trees on agriculture patch).   !
                     !          Perhaps when we convert bseeds to repro, we should ask     !
                     !          how much of the seed pool should be lost to harvest.       !
                     !---------------------------------------------------------------------!
                     csite%fast_soil_N(ipa) = csite%fast_soil_N(ipa)                       &
                                            + csite%repro(ipft,ipa) / c2n_recruit(ipft)
                     csite%fast_soil_C(ipa) = csite%fast_soil_C(ipa)                       &
                                            + csite%repro(ipft,ipa)
                     csite%repro(ipft,ipa)  = 0.0
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!
               end do pftloop
               !---------------------------------------------------------------------------!


               !----- Update the number of cohorts with the recently created. -------------!
               ncohorts_new = cpatch%ncohorts + inew
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     The number of recruits is now known. If there is any recruit, then we !
               ! allocate the temporary patch vector with the current number plus the      !
               ! number of recruits.                                                       !
               !---------------------------------------------------------------------------!
               if (ncohorts_new > cpatch%ncohorts) then
                  nullify(temppatch)
                  allocate(temppatch)
                  call allocate_patchtype(temppatch,cpatch%ncohorts)

                  !----- Fill the temp space with the current patches. --------------------!
                  call copy_patchtype(cpatch,temppatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
                  !------------------------------------------------------------------------!

                  !----- Deallocate the current patch. ------------------------------------!
                  call deallocate_patchtype(cpatch)
                  !------------------------------------------------------------------------!

                  !----- Reallocate the current site. -------------------------------------!
                  call allocate_patchtype(cpatch,ncohorts_new)
                  !------------------------------------------------------------------------!

                  !----- Transfer the temp values back in. --------------------------------!
                  call copy_patchtype(temppatch,cpatch,1,temppatch%ncohorts                &
                                     ,1,temppatch%ncohorts)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Add recruits.                                                      !
                  !------------------------------------------------------------------------!
                  inew = 0
                  recloop: do ico = temppatch%ncohorts+1,ncohorts_new
                     !---------------------------------------------------------------------!
                     !     Add the recruits, copying the information from the recruitment  !
                     ! table, and derive other variables or assume standard initial        !
                     ! values.                                                             !
                     !---------------------------------------------------------------------!
                     inew = inew + 1

                     !----- Copy from recruitment table (I). ------------------------------!
                     cpatch%pft(ico)       = recruit(inew)%pft
                     cpatch%hite(ico)      = recruit(inew)%hite
                     cpatch%dbh(ico)       = recruit(inew)%dbh
                     !---------------------------------------------------------------------!

                     !----- Carry out standard initialization. ----------------------------!
                     call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi))
                     !---------------------------------------------------------------------!


                     !----- Copy from recruitment table (II). -----------------------------!
                     cpatch%nplant          (ico) = recruit(inew)%nplant 
                     cpatch%bdead           (ico) = recruit(inew)%bdead 
                     cpatch%paw_avg         (ico) = recruit(inew)%paw_avg
                     cpatch%elongf          (ico) = recruit(inew)%elongf
                     cpatch%phenology_status(ico) = recruit(inew)%phenology_status
                     cpatch%bleaf           (ico) = recruit(inew)%bleaf
                     cpatch%broot           (ico) = recruit(inew)%broot
                     cpatch%bsapwooda       (ico) = recruit(inew)%bsapwooda
                     cpatch%bsapwoodb       (ico) = recruit(inew)%bsapwoodb
                     cpatch%balive          (ico) = recruit(inew)%balive
                     cpatch%bstorage        (ico) = recruit(inew)%bstorage
                     cpatch%leaf_temp       (ico) = recruit(inew)%leaf_temp
                     cpatch%wood_temp       (ico) = recruit(inew)%wood_temp
                     cpatch%leaf_temp_pv    (ico) = recruit(inew)%leaf_temp_pv
                     cpatch%wood_temp_pv    (ico) = recruit(inew)%wood_temp_pv
                     cpatch%leaf_vpdef      (ico) = recruit(inew)%leaf_vpdef
                     !---------------------------------------------------------------------!


                     !----- Initialise the next variables with zeroes... ------------------!
                     cpatch%leaf_water(ico) = 0.0
                     cpatch%leaf_fliq (ico) = 0.0
                     cpatch%wood_water(ico) = 0.0
                     cpatch%wood_fliq (ico) = 0.0
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    Compute initial AGB and Basal Area.  Their derivatives will be   !
                     ! zero.                                                               !
                     !---------------------------------------------------------------------!
                     cpatch%agb      (ico) = ed_biomass( cpatch%bdead     (ico)            &
                                                       , cpatch%bleaf     (ico)            &
                                                       , cpatch%bsapwooda (ico)            &
                                                       , cpatch%pft       (ico) )
                     cpatch%basarea  (ico) = pio4 * cpatch%dbh(ico)  * cpatch%dbh(ico)
                     cpatch%dagb_dt  (ico) = 0.0
                     cpatch%dlnagb_dt(ico) = 0.0
                     cpatch%dba_dt   (ico) = 0.0
                     cpatch%dlnba_dt (ico) = 0.0
                     cpatch%ddbh_dt  (ico) = 0.0
                     cpatch%dlndbh_dt(ico) = 0.0
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !     Set new_recruit_flag to 1 indicates that this cohort is         !
                     ! included when we tally agb_recruit, basal_area_recruit.             !
                     !---------------------------------------------------------------------!
                     cpatch%new_recruit_flag(ico) = 1
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    Obtain derived properties.                                       !
                     !---------------------------------------------------------------------!
                     !----- Find LAI, WAI, and CAI. ---------------------------------------!
                     call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico)                &
                                      ,cpatch%bdead(ico),cpatch%balive(ico)                &
                                      ,cpatch%dbh(ico),cpatch%hite(ico),cpatch%pft(ico)    &
                                      ,cpatch%sla(ico),cpatch%lai(ico),cpatch%wai(ico)     &
                                      ,cpatch%crown_area(ico),cpatch%bsapwooda(ico))
                     !----- Find heat capacity and vegetation internal energy. ------------!
                     call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdead(ico)                &
                                       ,cpatch%bsapwooda(ico),cpatch%nplant(ico)           &
                                       ,cpatch%pft(ico)                                    &
                                       ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))

                     cpatch%leaf_energy(ico) = cmtl2uext(cpatch%leaf_hcap (ico)            &
                                                        ,cpatch%leaf_water(ico)            &
                                                        ,cpatch%leaf_temp (ico)            &
                                                        ,cpatch%leaf_fliq (ico))
                     cpatch%wood_energy(ico) = cmtl2uext(cpatch%wood_hcap (ico)            &
                                                        ,cpatch%wood_water(ico)            &
                                                        ,cpatch%wood_temp (ico)            &
                                                        ,cpatch%wood_fliq (ico))

                     call is_resolvable(csite,ipa,ico)
                     !---------------------------------------------------------------------!

                     !----- Update number of cohorts in this site. ------------------------!
                     csite%cohort_count(ipa) = csite%cohort_count(ipa) + 1
                     !---------------------------------------------------------------------!
                  end do recloop
                  !------------------------------------------------------------------------!



                  !---- Remove the temporary patch. ---------------------------------------!
                  call deallocate_patchtype(temppatch)
                  deallocate(temppatch)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !   Now that recruitment has occured, terminate, fuse, split, and re-sort.     !
            !------------------------------------------------------------------------------!
            update_patch_loop: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Update the cohort distribution. -------------------------------------!
               if(cpatch%ncohorts > 0 .and. maxcohort >= 0) then
                  call terminate_cohorts(csite,ipa,elim_nplant,elim_lai)
                  call fuse_cohorts(csite,ipa, cpoly%green_leaf_factor(:,isi)              &
                                   ,cpoly%lsl(isi),.false.)
                  call split_cohorts(cpatch, cpoly%green_leaf_factor(:,isi),cpoly%lsl(isi))
               end if
               !---------------------------------------------------------------------------!


               !----- Sort cohorts by height. ---------------------------------------------!
               call sort_cohorts(cpatch)
               !---------------------------------------------------------------------------!

               !----- Update the number of cohorts (this is redundant...). ----------------!
               csite%cohort_count(ipa) = cpatch%ncohorts
               !---------------------------------------------------------------------------!


               !----- Since cohorts may have changed, update patch properties... ----------!
               call update_patch_derived_props(csite,ipa)
               call update_budget(csite,cpoly%lsl(isi),ipa,ipa)
               !---------------------------------------------------------------------------!
            end do update_patch_loop
            !------------------------------------------------------------------------------!

            !----- Since patch properties may have changed, update site properties... -----!
            call update_site_derived_props(cpoly,0,isi)
            !------------------------------------------------------------------------------!


            !----- Reset minimum monthly temperature. -------------------------------------!
            cpoly%min_monthly_temp(isi) = huge(1.)
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!



   case (1)
      !------------------------------------------------------------------------------------!
      !                                    'big leaf' ED                                   !
      !  Growth and reproduction are done together as there is no vertical structure in    !
      ! big leaf ED so the cohorts cannot grow vertically.  Therefore daily NPP is         !
      ! is accumulated and monthly the nplant of each cohort is increased (1 cohort per    !
      ! patch and 1 patch per pft and disturbance type).
      !------------------------------------------------------------------------------------!

      
      !----- The big loops start here. ----------------------------------------------------!
      polyloop_big: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         late_spring = (cgrid%lat(ipy) >= 0.0 .and. month == 6) .or.                       &
                       (cgrid%lat(ipy) < 0.0 .and. month == 12)

         !------- Update the repro arrays. ------------------------------------------------!
         call seed_dispersal(cpoly,late_spring)
         !---------------------------------------------------------------------------------!

         siteloop_big: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            patchloop_big: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !---------------------------------------------------------------------------!
               !    There should only be ONE cohort... if there are more, crash            !
               !---------------------------------------------------------------------------!
               if (cpatch%ncohorts > 1) then
                 write (unit=*,fmt='(a,1x,es12.5)') ' + PATCH   : ',ipa
                 write (unit=*,fmt='(a,1x,es12.5)') ' + NCOHORTS: ',cpatch%ncohorts
                 call fatal_error('NCOHORTS can never be greater than 1 for big-leaf runs' &
                                  ,'reproduction' ,'reproduction.f90')
               end if 
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      "Loop" over cohorts.  Reproduction does not create new patches,      !
               ! instead it will add population to the existing cohort.                    !
               !---------------------------------------------------------------------------!
               cohortloop_big: do ico = 1, cpatch%ncohorts  


                  !------ Current PFT. ----------------------------------------------------!
                  ipft = cpatch%pft(ico)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Check whether to include this PFT or not.  The decision depends    !
                  ! on the following decisions.                                            !
                  ! 1.  In case of agriculture or forest plantation, it must also be       !
                  ! allowed in such patches.                                               !
                  ! 2.  The temperature is not limiting reproduction.                      !
                  ! 3.  The user wants reproduction to occur.                              !
                  !------------------------------------------------------------------------!
                  select case (csite%dist_type(ipa))
                  case (1)
                     !----- Agriculture (cropland or pasture). ----------------------------!
                     allow_pft =                                                           &
                           include_pft_ag(ipft)                                      .and. &
                           cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and. &
                           repro_scheme                /= 0
                     !---------------------------------------------------------------------!

                  case (2)
                     !----- Forest plantation. --------------------------------------------!
                     allow_pft =                                                           &
                           include_pft_fp(ipft)                                      .and. &
                           cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and. &
                           repro_scheme                /= 0
                     !---------------------------------------------------------------------!

                  case default
                     !----- Primary or secondary vegetation. ------------------------------!
                     allow_pft =                                                           &
                           cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and. &
                           repro_scheme                /= 0
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Update cohort properties in case reproduction is allowed.          !
                  !------------------------------------------------------------------------!
                  if (allow_pft) then

                     !---------------------------------------------------------------------!
                     !    Plants don't have size distribution, so use the standard value   !
                     ! of one plant to find a population increase that is consistent with  !
                     ! the expected average biomass.                                       !
                     !---------------------------------------------------------------------!
                     nplant_inc         = csite%repro(ipft,ipa) / one_plant_c(ipft)
                     cpatch%nplant(ico) = cpatch%nplant(ico) + nplant_inc
                     !---------------------------------------------------------------------!


                     !----- Reset the carbon available for reproduction. ------------------!
                     csite%repro(ipft,ipa) = 0.0
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !    Will only reproduce/grow if on-allometry so dont' have to worry  !
                     ! about elongation factor.                                            !
                     !---------------------------------------------------------------------!
                     bleaf_plant     = size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft) 
                     broot_plant     = bleaf_plant * q(ipft)
                     bsapwood_plant  = bleaf_plant * qsw(ipft) * cpatch%hite(ico)
                     balive_plant    = bleaf_plant + broot_plant + bsapwood_plant
                     bdead_plant     = dbh2bd(cpatch%dbh(ico),ipft)
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !      Update the productivity terms.                                 !
                     !---------------------------------------------------------------------!
                     cpatch%today_nppleaf(ico)   =  nplant_inc * bleaf_plant
                     cpatch%today_nppfroot(ico)  =  nplant_inc * broot_plant
                     cpatch%today_nppsapwood(ico)=  nplant_inc * bsapwood_plant
                     cpatch%today_nppwood(ico)   = agf_bs(ipft) * nplant_inc * bdead_plant
                     cpatch%today_nppcroot(ico)  = (1. - agf_bs(ipft)) * nplant_inc        &
                                                 * bdead_plant
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     Update the derived properties since the population may have     !
                     ! changed.                                                            !
                     !---------------------------------------------------------------------!
                     !----- Find LAI, WAI, and CAI. ---------------------------------------!
                     call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico)                &
                                      ,cpatch%bdead(ico),cpatch%balive(ico)                &
                                      ,cpatch%dbh(ico),cpatch%hite(ico)                    &
                                      ,cpatch%pft(ico),cpatch%sla(ico)                     &
                                      ,cpatch%lai(ico),cpatch%wai(ico)                     &
                                      ,cpatch%crown_area(ico),cpatch%bsapwooda(ico))
                     !----- Find heat capacity and vegetation internal energy. ------------!
                     call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdead(ico)                &
                                       ,cpatch%bsapwooda(ico),cpatch%nplant(ico)           &
                                       ,cpatch%pft(ico)                                    &
                                       ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))
                     cpatch%leaf_energy(ico) = cmtl2uext(cpatch%leaf_hcap (ico)            &
                                                        ,cpatch%leaf_water(ico)            &
                                                        ,cpatch%leaf_temp (ico)            &
                                                        ,cpatch%leaf_fliq (ico))
                     cpatch%wood_energy(ico) = cmtl2uext(cpatch%wood_hcap (ico)            &
                                                        ,cpatch%wood_water(ico)            &
                                                        ,cpatch%wood_temp (ico)            &
                                                        ,cpatch%wood_fliq (ico))
                     !----- Update flags for the biophysical integrator. ------------------!
                     call is_resolvable(csite,ipa,ico)
                     !---------------------------------------------------------------------!
                  else
                     !---------------------------------------------------------------------!
                     !     This PFT shouldn't exist... at least not on this patch.   Send  !
                     ! the seed litter to the soil pools for decomposition.                !
                     !                                                                     !
                     !     ALS: dont send all seeds to litter!  Keep it for harvesting?    !
                     !     MLO: not here because these are the PFTs that are not allowed   !
                     !          on this patch (for example, trees on agriculture patch).   !
                     !          Perhaps when we convert bseeds to repro, we should ask     !
                     !          how much of the seed pool should be lost to harvest.       !
                     !---------------------------------------------------------------------!
                     csite%fast_soil_N(ipa) = csite%fast_soil_N(ipa)                       &
                                            + csite%repro(ipft,ipa) / c2n_recruit(ipft)
                     csite%fast_soil_C(ipa) = csite%fast_soil_C(ipa)                       &
                                            + csite%repro(ipft,ipa)
                     csite%repro(ipft,ipa)  = 0.0
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!
               end do cohortloop_big
               !---------------------------------------------------------------------------!
            end do patchloop_big
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Update derived properties.                                               !
            !------------------------------------------------------------------------------!
            update_patch_loop_big: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Update the number of cohorts (this is redundant...). ----------------!
               csite%cohort_count(ipa) = cpatch%ncohorts
               !---------------------------------------------------------------------------!

               !----- Since cohorts may have changed, update patch properties... ----------!
               call update_patch_derived_props(csite,ipa)
               call update_budget(csite,cpoly%lsl(isi),ipa,ipa)
               !---------------------------------------------------------------------------!
            end do update_patch_loop_big
            !------------------------------------------------------------------------------!

            !----- Since patch properties may have changed, update site properties... -----!
            call update_site_derived_props(cpoly,0,isi)
            !------------------------------------------------------------------------------!

            !----- Reset minimum monthly temperature. -------------------------------------!
            cpoly%min_monthly_temp(isi) = huge(1.)
            !------------------------------------------------------------------------------!
         end do siteloop_big
         !---------------------------------------------------------------------------------!
      end do polyloop_big
      !------------------------------------------------------------------------------------!
   end select
   !---------------------------------------------------------------------------------------!

   return
end subroutine reproduction
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!      This subroutine will bypass all reproduction.                                       !
!------------------------------------------------------------------------------------------!
subroutine reproduction_eq_0(cgrid, month)
   use ed_state_vars      , only : edtype                & ! structure
                                 , polygontype           & ! structure
                                 , sitetype              & ! structure
                                 , patchtype             & ! structure
                                 , allocate_patchtype    & ! subroutine
                                 , copy_patchtype        & ! subroutine
                                 , deallocate_patchtype  ! ! subroutine
   use pft_coms           , only : recruittype           & ! structure
                                 , zero_recruit          & ! subroutine
                                 , copy_recruit          & ! subroutine
                                 , seedling_mortality    & ! intent(in)
                                 , c2n_stem              & ! intent(in)
                                 , l2n_stem              & ! intent(in)
                                 , min_recruit_size      & ! intent(in)
                                 , c2n_recruit           & ! intent(in)
                                 , seed_rain             & ! intent(in)
                                 , include_pft           & ! intent(in)
                                 , include_pft_ag        & ! intent(in)
                                 , qsw                   & ! intent(in)
                                 , q                     & ! intent(in)
                                 , sla                   & ! intent(in)
                                 , hgt_min               & ! intent(in)
                                 , plant_min_temp        ! ! intent(in)
   use decomp_coms        , only : f_labile              ! ! intent(in)
   use ed_max_dims        , only : n_pft                 ! ! intent(in)
   use fuse_fiss_utils    , only : sort_cohorts          & ! subroutine
                                 , terminate_cohorts     & ! subroutine
                                 , fuse_cohorts          & ! subroutine
                                 , split_cohorts         ! ! subroutine
   use phenology_coms     , only : repro_scheme          ! ! intent(in)
   use mem_polygons       , only : maxcohort             ! ! intent(in)
   use consts_coms        , only : pio4                  ! ! intent(in)
   use ed_therm_lib       , only : calc_veg_hcap         ! ! function
   use grid_coms          , only : nzg                   ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)     , target     :: cgrid
   integer          , intent(in) :: month
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype), pointer          :: cpoly
   type(sitetype)   , pointer          :: csite
   type(patchtype)  , pointer          :: cpatch
   type(patchtype)  , pointer          :: temppatch
   type(recruittype), dimension(n_pft) :: recruit
   type(recruittype)                   :: rectest
   integer                             :: ipy
   integer                             :: isi
   integer                             :: ipa
   integer                             :: ico
   !----- Saved variables -----------------------------------------------------------------!
   logical          , save             :: first_time = .true.
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    If this is the first time, check whether the user wants reproduction.  If not,     !
   ! kill all potential recruits and send their biomass to the litter pool.                !
   !---------------------------------------------------------------------------------------!
   seedling_mortality(1:n_pft) = 1.0
   !---------------------------------------------------------------------------------------!


   !----- The big loops start here. -------------------------------------------------------!
   polyloop: do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      siteloop: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         !---------------------------------------------------------------------------------!
         !    Zero all reproduction stuff...                                               !
         !---------------------------------------------------------------------------------!
         patchloop: do ipa = 1,csite%npatches
            call zero_recruit(n_pft,recruit)
            csite%repro(:,ipa)  = 0.0
         end do patchloop
         !---------------------------------------------------------------------------------!

         !----- Reset minimum monthly temperature. ----------------------------------------!
         cpoly%min_monthly_temp(isi) = huge(1.)
         !---------------------------------------------------------------------------------!
      end do siteloop
      !------------------------------------------------------------------------------------!
   end do polyloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine reproduction_eq_0
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine will update the repro structure.  The way the dispersal is going to !
! happen depends on whether the user allows dispersal between sites or within sites only.  !
! Broadleaf deciduous PFTs do their dispersal only in late spring.  This is because of     !
! their unique storage respiration.  Other PFTs disperse seeds at any time of the year.    !
!------------------------------------------------------------------------------------------!
subroutine seed_dispersal(cpoly,late_spring)
   use phenology_coms     , only : repro_scheme          ! ! intent(in)
   use ed_state_vars      , only : polygontype           & ! structure
                                 , sitetype              & ! structure
                                 , patchtype             ! ! structure
   use pft_coms           , only : recruittype           & ! structure
                                 , nonlocal_dispersal    & ! intent(in)
                                 , seedling_mortality    & ! intent(in)
                                 , phenology             ! ! intent(in)
   use ed_misc_coms       , only : ibigleaf              ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(polygontype), target     :: cpoly       ! Current polygon               [      ---]
   logical          , intent(in) :: late_spring ! Is it late spring             [      T|F]
   !----- Local variables. ----------------------------------------------------------------!
   type(sitetype)   , pointer    :: csite       ! Current site                  [      ---]
   type(sitetype)   , pointer    :: donsite     ! Donor site                    [      ---]
   type(sitetype)   , pointer    :: recsite     ! Receptor site                 [      ---]
   type(patchtype)  , pointer    :: donpatch    ! Donor patch                   [      ---]
   integer                       :: isi         ! Site counter                  [      ---]
   integer                       :: recsi       ! Receptor site counter         [      ---]
   integer                       :: donsi       ! Donor site counter            [      ---]
   integer                       :: recpa       ! Receptor patch counter        [      ---]
   integer                       :: donpa       ! Donor patch counter           [      ---]
   integer                       :: donco       ! Donor cohort counter          [      ---]
   integer                       :: donpft      ! Donor PFT                     [      ---]
   real                          :: bseedling   ! Surviving seedling biomass    [   kgC/m2]
   real                          :: bseed_stays ! Seedling biomass that stays   [   kgC/m2]
   real                          :: bseed_maygo ! Seedling biomass that may go  [   kgC/m2]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Here we decide how to disperse seeds based on the reproduction scheme.            !
   !---------------------------------------------------------------------------------------!
   select case (repro_scheme)
   case (0)
       !------ No reproduction, quit. -----------------------------------------------------!
       return
   case (1)
      !------------------------------------------------------------------------------------!
      !     Seeds are dispersed amongst patches that belong to the same site, but they     !
      ! cannot go outside their native site.                                               !
      !------------------------------------------------------------------------------------!
      siteloop1: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         !---------------------------------------------------------------------------------!
         !      Loop over the donor cohorts.                                               !
         !---------------------------------------------------------------------------------!
         donpaloop1: do donpa = 1,csite%npatches
            donpatch => csite%patch(donpa)

            doncoloop1: do donco = 1, donpatch%ncohorts

               !----- Define an alias for PFT. --------------------------------------------!
               donpft = donpatch%pft(donco)

               !---------------------------------------------------------------------------!
               !    Find the biomass of survivor seedlings.  Units: kgC/m2                 !
               !---------------------------------------------------------------------------!
               if (phenology(donpft) /= 2 .or. late_spring) then
                  bseedling   = donpatch%nplant(donco) * donpatch%bseeds(donco)            &
                              * (1.0 - seedling_mortality(donpft))
                  select case (ibigleaf)
                  case (0)
                     bseed_stays = bseedling * (1.0 - nonlocal_dispersal(donpft))
                     bseed_maygo = bseedling * nonlocal_dispersal(donpft)
                  case (1)
                     !---- if bigleaf cannot disperse seedlings to other patches ----------!
                     bseed_stays = bseedling
                     bseed_maygo = 0.
                  end select

               else
                  !----- Not a good time for reproduction.  No seedlings. -----------------!
                  bseedling     = 0.
                  bseed_stays   = 0.
                  bseed_maygo   = 0.
               end if
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !   Spread the seedlings across all patches in this site.                   !
               !---------------------------------------------------------------------------!
               recpaloop1: do recpa = 1,csite%npatches

                  !------------------------------------------------------------------------!
                  !     Add the non-local dispersal evenly across all patches, including   !
                  ! the donor patch.  We must scale the biomass by the combined area of    !
                  ! this patch and site so the total carbon is preserved.                  !
                  !------------------------------------------------------------------------!
                  csite%repro(donpft,recpa) = csite%repro(donpft,recpa)                    &
                                            + bseed_maygo * csite%area(donpa)
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !      Include the local dispersal if this is the donor patch.           !
                  !------------------------------------------------------------------------!
                  if (recpa == donpa) then
                     csite%repro(donpft,recpa) = csite%repro(donpft,recpa) + bseed_stays
                  end if
                  !------------------------------------------------------------------------!

               end do recpaloop1
               !---------------------------------------------------------------------------!
            end do doncoloop1
            !------------------------------------------------------------------------------!
         end do donpaloop1
         !---------------------------------------------------------------------------------!
      end do siteloop1 
      !------------------------------------------------------------------------------------!

   case (2,3)

      !------------------------------------------------------------------------------------!
      !     Seeds are dispersed amongst patches that belong to the same polygon.  They are !
      ! allowed to go from one site to the other, but they cannot go outside their native  !
      ! polygon.                                                                           !
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Loop over the donor cohorts.                                                  !
      !------------------------------------------------------------------------------------!
      donsiloop2: do donsi = 1,cpoly%nsites
         donsite => cpoly%site(donsi)

         donpaloop2: do donpa = 1,donsite%npatches
            donpatch => donsite%patch(donpa)

            doncoloop2: do donco = 1, donpatch%ncohorts

               !----- Define an alias for PFT. --------------------------------------------!
               donpft = donpatch%pft(donco)

               !---------------------------------------------------------------------------!
               !    Find the biomass of survivor seedlings.  Units: kgC/m2                 !
               !---------------------------------------------------------------------------!
               if (phenology(donpft) /= 2 .or. late_spring) then
                  bseedling   = donpatch%nplant(donco) * donpatch%bseeds(donco)            &
                              * (1.0 - seedling_mortality(donpft))

                  select case (ibigleaf)
                  case (0)
                     bseed_stays = bseedling * (1.0 - nonlocal_dispersal(donpft))
                     bseed_maygo = bseedling * nonlocal_dispersal(donpft)
                  case (1)
                     !---- if bigleaf cannot disperse seedlings to other patches ----------!
                     bseed_stays = bseedling
                     bseed_maygo = 0.
                  end select
               else
                  !----- Not a good time for reproduction.  No seedlings. -----------------!
                  bseedling   = 0.
                  bseed_stays = 0.
                  bseed_maygo = 0.
               end if
               !---------------------------------------------------------------------------!
               
               !---------------------------------------------------------------------------!
               !   Spread the seedlings across all patches in this polygon.                !
               !---------------------------------------------------------------------------!
               recsiloop2: do recsi = 1,cpoly%nsites
                  recsite => cpoly%site(recsi)
                  recpaloop2: do recpa = 1,recsite%npatches

                     !---------------------------------------------------------------------!
                     !     Add the non-local dispersal evenly across all patches,          !
                     ! including the donor patch.  We must scale the biomass by the        !
                     ! combined area of this patch and site so the total carbon is         !
                     ! preserved.                                                          !
                     !---------------------------------------------------------------------!
                     recsite%repro(donpft,recpa) = recsite%repro(donpft,recpa)             &
                                                 + bseed_maygo * recsite%area(donpa)       &
                                                 * cpoly%area(donsi)
                     !---------------------------------------------------------------------!

                     !---------------------------------------------------------------------!
                     !      Include the local dispersal if this is the donor patch.        !
                     !---------------------------------------------------------------------!
                     if (recpa == donpa .and. recsi == donsi) then
                        recsite%repro(donpft,recpa) = recsite%repro(donpft,recpa)          &
                                                    + bseed_stays
                     end if
                     !---------------------------------------------------------------------!
                     
                  end do recpaloop2
                  !------------------------------------------------------------------------!
               end do recsiloop2
               !---------------------------------------------------------------------------!
            end do doncoloop2
            !------------------------------------------------------------------------------!
         end do donpaloop2
         !---------------------------------------------------------------------------------!
      end do donsiloop2
      !------------------------------------------------------------------------------------!


   end select
   return
end subroutine seed_dispersal
!==========================================================================================!
!==========================================================================================!
