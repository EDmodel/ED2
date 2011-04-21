!==========================================================================================!
!==========================================================================================!
!      This subroutine will drive the reproduction, based on its carbon availability and   !
! PFT-specific reproduction properties.  No reproduction will happen if the user didn't    !
! want it, in which case the seedling biomass will go to the litter pools.                 !
!------------------------------------------------------------------------------------------!
subroutine reproduction(cgrid, month)
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
                                 , nonlocal_dispersal    & ! intent(in)
                                 , seedling_mortality    & ! intent(in)
                                 , phenology             & ! intent(in)
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
   use ed_therm_lib       , only : calc_hcapveg          ! ! function
   use allometry          , only : dbh2bd                & ! function
                                 , dbh2bl                & ! function
                                 , h2dbh                 & ! function
                                 , ed_biomass            & ! function
                                 , area_indices          ! ! subroutine
   use grid_coms          , only : nzg                   ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)     , target     :: cgrid
   integer          , intent(in) :: month
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype), pointer          :: cpoly
   type(sitetype)   , pointer          :: csite
   type(patchtype)  , pointer          :: cpatch, dpatch, temppatch
   type(recruittype), dimension(n_pft) :: recruit
   type(recruittype)                   :: rectest
   integer                             :: ipy, isi, ipa, ico, ipft
   integer                             :: recp, donp
   integer                             :: inew, ncohorts_new
   logical                             :: late_spring
   real                                :: elim_nplant
   real                                :: elim_lai
   !----- Saved variables -----------------------------------------------------------------!
   logical          , save             :: first_time=.true.
   !----- External functions. -------------------------------------------------------------!
   logical          , external         :: is_resolvable  ! The cohort can be resolved.
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    If this is the first time, check whether the user wants reproduction.  If not,     !
   ! kill all potential recruits and send their biomass to the litter pool.                !
   !---------------------------------------------------------------------------------------!
   if (first_time) then
      if (repro_scheme == 0) seedling_mortality(1:n_pft) = 1.0
      first_time = .false.
   end if


   !----- The big loops start here. -------------------------------------------------------!
   polyloop: do ipy = 1,cgrid%npolygons
   
      !------------------------------------------------------------------------------------!
      !     Check whether this is late spring/early summer.  This is needed for temperate  !
      ! broadleaf deciduous trees.  Late spring means June in the Northern Hemisphere, or  !
      ! December in the Southern Hemisphere.                                               !
      !------------------------------------------------------------------------------------!
      late_spring = (cgrid%lat(ipy) >= 0.0 .and. month == 6) .or.                          &
                    (cgrid%lat(ipy) < 0.0 .and. month == 12)

      cpoly => cgrid%polygon(ipy)
      siteloop: do isi = 1,cpoly%nsites
         csite => cpoly%site(isi)

         !---------------------------------------------------------------------------------!
         !      Cohorts may have grown differently, so we need to sort them by size.       !
         !---------------------------------------------------------------------------------!
         sortloop: do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)
            call sort_cohorts(cpatch)
         end do sortloop

         !---------------------------------------------------------------------------------!
         !   We now perform the seed dispersal.  Some considerations:                      !
         !   1.  For non-local dispersal, seeds are dispersed throughout the site, not the !
         !       polygon;                                                                  !
         !   2.  Note that broadleaf deciduous PFTs do their dispersal only in late        !
         !       spring.  This is because of their unique storage respiration.  Other PFTs !
         !       do dispersal any time of the year.                                        !
         !---------------------------------------------------------------------------------!

         recploop: do recp = 1,csite%npatches     !recp = ptarget
            
            donploop: do donp = 1,csite%npatches  !donp = psource

               
               !----- Non-local, gridcell-wide dispersal. ---------------------------------!
               dpatch => csite%patch(donp)
               doncloop: do ico=1,dpatch%ncohorts

                  !----- Defining an alias for PFT. ---------------------------------------!
                  ipft = dpatch%pft(ico)

                  if(phenology(ipft) /= 2 .or. late_spring )then
                     !---------------------------------------------------------------------!
                     !    This means that we are not dealing with broad leaf deciduous,    !
                     ! or, if it is, then were are in the reproduction season.  Fill the   !
                     ! reproduction table.                                                 !
                     !---------------------------------------------------------------------!
                     csite%repro(ipft,recp) = csite%repro(ipft,recp)                       &
                                            + nonlocal_dispersal(ipft)*dpatch%nplant(ico)  &
                                            * (1.0 - seedling_mortality(ipft))             &
                                            * dpatch%bseeds(ico) * csite%area(recp)
                  end if
               end do doncloop
               
               !----- Local dispersal (seeds stay in this patch). -------------------------!
               if (recp == donp) then

                  cpatch => csite%patch(donp)
                  selfcloop: do ico=1,cpatch%ncohorts

                     !----- Defining an alias for PFT. ------------------------------------!
                     ipft = cpatch%pft(ico)

                     if (phenology(ipft) /= 2 .or. late_spring) then
                        csite%repro(ipft,recp) = csite%repro(ipft,recp)                    &
                             + cpatch%nplant(ico) * (1.0 - nonlocal_dispersal(ipft))       &
                             * ( 1.0 - seedling_mortality(ipft) ) * cpatch%bseeds(ico)
                     end if
                  end do selfcloop
               end if
            end do donploop
         end do recploop

         !---------------------------------------------------------------------------------!
         !    For the recruitment to happen, four requirements must be met:                !
         !    1.  PFT is included in this simulation;                                      !
         !    2.  It is not too cold (min_monthly_temp > plant_min_temp - 5)               !
         !    3.  We are dealing with EITHER a non-agriculture patch OR                    !
         !        a PFT that could exist in an agricultural patch.                         !
         !    4.  There must be sufficient carbon to form the recruits.                    !
         !---------------------------------------------------------------------------------!
         patchloop: do ipa = 1,csite%npatches
            inew = 0
            call zero_recruit(n_pft,recruit)
            cpatch => csite%patch(ipa)

            !---- This time we loop over PFTs, not cohorts. -------------------------------!
            pftloop: do ipft = 1, n_pft

               !---------------------------------------------------------------------------!
               !    Check to make sure we are including the PFT and that it is not too     !
               ! cold.                                                                     !
               !---------------------------------------------------------------------------!
               if(include_pft(ipft)           == 1                          .and.          &
                  cpoly%min_monthly_temp(isi) >= plant_min_temp(ipft) - 5.0 .and.          &
                  repro_scheme == 1 ) then

                  !------------------------------------------------------------------------!
                  !     Make sure that this is not agriculture or that it is okay for this !
                  ! PFT to be in an agriculture patch.                                     !
                  !------------------------------------------------------------------------!
                  if(csite%dist_type(ipa) /= 1 .or. include_pft_ag(ipft) == 1) then

                     !---------------------------------------------------------------------!
                     !    We assign the recruit in the temporary recruitment structure.    !
                     !---------------------------------------------------------------------!
                     rectest%pft      = ipft
                     rectest%veg_temp = csite%can_temp(ipa)
                     rectest%hite     = hgt_min(ipft)
                     rectest%dbh      = h2dbh(rectest%hite, ipft)
                     rectest%bdead    = dbh2bd(rectest%dbh, rectest%hite, ipft)
                     rectest%bleaf    = dbh2bl(rectest%dbh, ipft)
                     rectest%balive   = rectest%bleaf                                      &
                                      * (1.0 + q(ipft) + qsw(ipft) * rectest%hite)
                     rectest%nplant   = csite%repro(ipft,ipa)                              &
                                      / (rectest%balive + rectest%bdead)

                     if(include_pft(ipft) == 1) then
                        rectest%nplant = rectest%nplant + seed_rain(ipft)
                     end if

                     ! If there is enough carbon, form the recruits.
                     if ( rectest%nplant * (rectest%balive + rectest%bdead) >              &
                          min_recruit_size(ipft)) then
                        inew = inew + 1
                        call copy_recruit(rectest,recruit(inew))

                        !----- Reset the carbon available for reproduction. ---------------!
                        csite%repro(ipft,ipa) = 0.0

                     end if
                  else
                     !---------------------------------------------------------------------!
                     !     If we have reached this branch, we are in an agricultural       !
                     ! patch.  Send the seed litter to the soil pools for decomposition.   !
                     !---------------------------------------------------------------------!
                     csite%fast_soil_N(ipa) = csite%fast_soil_N(ipa)                       &
                                            + csite%repro(ipft,ipa) / c2n_recruit(ipft)
                     csite%fast_soil_C(ipa) = csite%fast_soil_C(ipa)                       &
                                            + csite%repro(ipft,ipa)
                     csite%repro(ipft,ipa)  = 0.0
                  end if
               end if
            end do pftloop

            !----- Update the number of cohorts with the recently created. ----------------!
            ncohorts_new = cpatch%ncohorts + inew
            
            !------------------------------------------------------------------------------!
            !     The number of recruits is now known. If there is any recruit, then we    !
            ! allocate the temporary patch vector with the current number plus the number  !
            ! of recruits.                                                                 !
            !------------------------------------------------------------------------------!
            if (ncohorts_new > cpatch%ncohorts) then
               nullify(temppatch)
               allocate(temppatch)
               call allocate_patchtype(temppatch,cpatch%ncohorts)

               !----- Fill the temp space with the current patches. -----------------------!
               call copy_patchtype(cpatch,temppatch,1,cpatch%ncohorts,1,cpatch%ncohorts)

               !----- Deallocate the current patch. ---------------------------------------!
               call deallocate_patchtype(cpatch)

               !----- Reallocate the current site. ----------------------------------------!
               call allocate_patchtype(cpatch,ncohorts_new)

               !----- Transfer the temp values back in. -----------------------------------!
               call copy_patchtype(temppatch,cpatch,1,temppatch%ncohorts                   &
                                  ,1,temppatch%ncohorts)

               inew = 0
               recloop: do ico = temppatch%ncohorts+1,ncohorts_new
                  !------------------------------------------------------------------------!
                  !     Add the recruits, copying the information from the recruitment     !
                  ! table, and derive other variables or assume standard initial values.   !
                  !------------------------------------------------------------------------!
                  inew = inew + 1

                  !----- Copy from recruitment table (I). ---------------------------------!
                  cpatch%pft(ico)       = recruit(inew)%pft
                  cpatch%hite(ico)      = recruit(inew)%hite
                  cpatch%dbh(ico)       = recruit(inew)%dbh
                  !------------------------------------------------------------------------!

                  !----- Carry out standard initialization. -------------------------------!
                  call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi))
                  !------------------------------------------------------------------------!


                  !----- Copy from recruitment table (II). --------------------------------!
                  cpatch%bdead(ico)     = recruit(inew)%bdead
                  cpatch%nplant(ico)    = recruit(inew)%nplant
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Even though we brought leaf biomass and biomass of the active      !
                  ! tissues, we will make them consistent with the initial amount of water !
                  ! available.  This is done inside pheninit_alive_storage.                !
                  !------------------------------------------------------------------------!
                  call pheninit_balive_bstorage(nzg,csite,ipa,ico,cpoly%ntext_soil(:,isi))
                  !------------------------------------------------------------------------!


                  !----- Assign temperature after init_ed_cohort_vars... ------------------!
                  cpatch%veg_temp(ico)  = recruit(inew)%veg_temp

                  !----- Initialise the next variables with zeroes... ---------------------!
                  cpatch%veg_water(ico) = 0.0
                  cpatch%veg_fliq(ico)  = 0.0
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !    Computing initial AGB and Basal Area. Their derivatives will be     !
                  ! zero.                                                                  !
                  !------------------------------------------------------------------------!
                  cpatch%agb(ico)     = ed_biomass(cpatch%bdead(ico),cpatch%balive(ico)    &
                                                  ,cpatch%bleaf(ico),cpatch%pft(ico)       &
                                                  ,cpatch%hite(ico),cpatch%bstorage(ico)   &
                                                  ,cpatch%bsapwood(ico))
                  cpatch%basarea(ico) = pio4 * cpatch%dbh(ico)  * cpatch%dbh(ico)
                  cpatch%dagb_dt(ico) = 0.0
                  cpatch%dba_dt(ico)  = 0.0
                  cpatch%ddbh_dt(ico) = 0.0
                  !------------------------------------------------------------------------!
                  !     Setting new_recruit_flag to 1 indicates that this cohort is        !
                  ! included when we tally agb_recruit, basal_area_recruit.                !
                  !------------------------------------------------------------------------!
                  cpatch%new_recruit_flag(ico) = 1
                  
                  !------------------------------------------------------------------------!
                  !    Obtain derived properties.                                          !
                  !------------------------------------------------------------------------!
                  !----- Find LAI, WPA, WAI. ----------------------------------------------!
                  call area_indices(cpatch%nplant(ico),cpatch%bleaf(ico),cpatch%bdead(ico) &
                                   ,cpatch%balive(ico),cpatch%dbh(ico), cpatch%hite(ico)   &
                                   ,cpatch%pft(ico),cpatch%sla(ico), cpatch%lai(ico)       &
                                   ,cpatch%wpa(ico),cpatch%wai(ico)                        &
                                   ,cpatch%crown_area(ico),cpatch%bsapwood(ico))
                  !----- Find heat capacity and vegetation internal energy. ---------------!
                  cpatch%hcapveg(ico) = calc_hcapveg(cpatch%bleaf(ico),cpatch%bdead(ico)   &
                                                    ,cpatch%balive(ico),cpatch%nplant(ico) &
                                                    ,cpatch%hite(ico),cpatch%pft(ico)      &
                                                    ,cpatch%phenology_status(ico)          &
                                                    ,cpatch%bsapwood(ico))
                  cpatch%veg_energy(ico) = cpatch%hcapveg(ico) * cpatch%veg_temp(ico)
                  cpatch%resolvable(ico) = is_resolvable(csite,ipa,ico                     &
                                                        ,cpoly%green_leaf_factor(:,isi))

                  !----- Update number of cohorts in this site. ---------------------------!
                  csite%cohort_count(ipa) = csite%cohort_count(ipa) + 1
               end do recloop

               !---- Remove the temporary patch. ------------------------------------------!
               call deallocate_patchtype(temppatch)
               deallocate(temppatch)
            end if
         end do patchloop
         
         
         !---------------------------------------------------------------------------------!
         !   Now that recruitment has occured, terminate, fuse, split, and re-sort.        !
         !---------------------------------------------------------------------------------!
         update_patch_loop: do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)

            if(cpatch%ncohorts > 0 .and. maxcohort >= 0) then
               call terminate_cohorts(csite,ipa,elim_nplant,elim_lai)
               call fuse_cohorts(csite,ipa, cpoly%green_leaf_factor(:,isi),cpoly%lsl(isi))                         
               call split_cohorts(cpatch, cpoly%green_leaf_factor(:,isi),cpoly%lsl(isi))
            end if

            !----- Sort the cohorts by height. --------------------------------------------!
            call sort_cohorts(cpatch)

            !----- Update the number of cohorts (this is redundant...). -------------------!
            csite%cohort_count(ipa) = cpatch%ncohorts

            !----- Since cohorts may have changed, update patch properties... -------------!
            call update_patch_derived_props(csite,cpoly%lsl(isi),cpoly%met(isi)%prss,ipa)
            call update_budget(csite,cpoly%lsl(isi),ipa,ipa)
         end do update_patch_loop

         !----- Since patch properties may have changed, update site properties... --------!
         call update_site_derived_props(cpoly,0,isi)
         
         !----- Reset minimum monthly temperature. ----------------------------------------!
         cpoly%min_monthly_temp(isi) = huge(1.)
      end do siteloop
   end do polyloop
   return
end subroutine reproduction
!==========================================================================================!
!==========================================================================================!
