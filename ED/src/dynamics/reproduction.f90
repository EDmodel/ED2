subroutine reproduction_ar(cgrid, month)

  use ed_state_vars,only:edtype,polygontype,sitetype,patchtype,&
       allocate_patchtype,copy_patchtype,deallocate_patchtype
    
  use misc_coms, only: dtlsm
  use pft_coms, only: nonlocal_dispersal, seedling_mortality, phenology,  &
       sla, c2n_stem, l2n_stem, c2n_recruit, seed_rain, include_pft,  &
       include_pft_ag, qsw, q, hgt_min, plant_min_temp
  use decomp_coms, only: f_labile
  use fusion_fission_coms, only: min_recruit_size
  use max_dims, only: n_pft
  use fuse_fiss_utils_ar, only: sort_cohorts_ar, terminate_cohorts_ar,fuse_cohorts_ar &
                            ,split_cohorts_ar
  use phenology_coms,only:repro_scheme
  
  use consts_coms, only: t3ple
  use canopy_air_coms, only: hcapveg_ref,heathite_min
  use mem_sites, only: maxcohort
  use ed_therm_lib,only : calc_hcapveg
  use allometry, only: dbh2bd, dbh2bl, h2dbh
  
  implicit none

  type(edtype),target       :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch,temppatch
  integer :: ipy,isi,ipa,ico,ipa1,ipa2

  integer, intent(in) :: month

  real :: nplant
  real :: balive
  real :: bleaf
  real :: bdead
  real :: dbh
  real :: hite
  integer :: pft

  integer :: inew,ncohorts_new
  real,dimension(n_pft,9) :: recruit_array
  
  logical, save :: first_time=.true.

  if(repro_scheme .eq. 0) seedling_mortality(1:n_pft) = 1.0

  allocate(temppatch)

  do ipy = 1,cgrid%npolygons
     
     cpoly => cgrid%polygon(ipy)
     
     do isi = 1,cpoly%nsites
        
        csite => cpoly%site(isi)
        
        !!initialize
        csite%repro = 0.0
   
        !--------------------------
        ! First, sort the cohorts after the growth
        !----------------------------
        do ipa = 1,csite%npatches
           cpatch => csite%patch(ipa)
           call sort_cohorts_ar(cpatch)
        enddo

        !--------------------------------------------------------------------
        ! SEED DISPERSAL:  
        !    (1) For non-local dispersal, seeds are dispersed throughout the 
        !        site, not the polygon.
        !    (2) Note that broadleaf deciduous PFTs do their dispersal only 
        !        in June (northern hemisphere) or December (southern
        !        hemisphere).  This is because of their unique storage 
        !        respiration.  All other PFTs do dispersal every month.
        !--------------------------------------------------------------------

        do ipa1 = 1,csite%npatches     !ipa1 = ptarget
           
           do ipa2 = 1,csite%npatches  !ipa2 = psource

              
              ! Non-local, gridcell-wide dispersal
              cpatch => csite%patch(ipa2)
              do ico=1,cpatch%ncohorts
                 
                 if(phenology(cpatch%pft(ico)) /= 2   .or.  &
                      (cgrid%lat(ipy) >= 0.0 .and. month == 6) .or.   &
                      (cgrid%lat(ipy) < 0.0 .and. month == 12) )then
                    ! NOT broad leaf deciduous OR Jun in north OR Dec in south
                    csite%repro(cpatch%pft(ico),ipa1) = csite%repro(cpatch%pft(ico),ipa1) +   &
                         nonlocal_dispersal(cpatch%pft(ico)) * cpatch%nplant(ico) *   &
                         ( 1.0 - seedling_mortality(cpatch%pft(ico)) ) &
                         * cpatch%bseeds(ico) * csite%area(ipa1)

!                    print*,ipa1,ico,cpatch%pft(ico),csite%repro(cpatch%pft(ico),ipa1),"PARTS", &
!                         nonlocal_dispersal(cpatch%pft(ico)),cpatch%nplant(ico),   &
!                         ( 1.0 - seedling_mortality(cpatch%pft(ico)) ), &
!                         cpatch%bseeds(ico),csite%area(ipa2)


                 endif

              enddo
              
              ! Local dispersal (seeds stay in this patch)
              if (ipa1.eq.ipa2) then

                 cpatch => csite%patch(ipa2)
                 do ico=1,cpatch%ncohorts
                    if(phenology(cpatch%pft(ico)) /= 2 .or.  &
                         (cgrid%lat(ipy) >= 0.0 .and. month == 6) .or.   &
                         (cgrid%lat(ipy) < 0.0 .and. month == 12) )then

                       ! NOT broadleaf decid OR Jun in north OR Dec in south
                       csite%repro(cpatch%pft(ico),ipa1) = csite%repro(cpatch%pft(ico),ipa1) +   &
                            cpatch%nplant(ico) * (1.0 - nonlocal_dispersal(cpatch%pft(ico)))   &
                            * ( 1.0 - seedling_mortality(cpatch%pft(ico)) ) * cpatch%bseeds(ico)
                    endif
                 enddo

              endif
           enddo
        enddo

        !--------------------------------------------------------------------
        ! RECRUITMENT: Form new recruits if:
        !   a) PFT is included in this simulation
        !   b) it is not winter (min_monthly_temp > plant_min_temp - 5)
        !   c) we are dealing with EITHER a non-agriculture patch OR 
        !      a PFT that could exist in an agricultural patch.
        !   d) finally, there must be sufficient carbon to form the recruits.
        !--------------------------------------------------------------------

        do ipa = 1,csite%npatches     !ipa1 = ptarget
           
           inew = 0
           cpatch => csite%patch(ipa)
           recruit_array = 0.0
           
           do pft = 1, n_pft
           
              ! Check to make sure we are including the PFT and that 
              ! it is not too cold.

              if(include_pft(pft) == 1 .and.   &
                   cpoly%min_monthly_temp(isi) >= plant_min_temp(pft) - 5.0 .and. &
                   repro_scheme > 0)then
                 
                 ! Make sure that this is not agriculture or that it is 
                 ! OK for this PFT to be in an agriculture patch.
                 if(csite%dist_type(ipa) /= 1 .or. include_pft_ag(pft) == 1)then

                    ! Generate specs for this PFT
                    hite = hgt_min(pft)
                    dbh = h2dbh(hite, pft)
                    bdead = dbh2bd(dbh, hite, pft)
                    bleaf = dbh2bl(dbh, pft)
                    balive = bleaf * (1.0 + q(pft) + qsw(pft) * hite)

                    nplant = csite%repro(pft,ipa) / (balive + bdead)

                    if(include_pft(pft) == 1) nplant = nplant + seed_rain(pft)

                    ! If there is enough carbon, form the recruits.
                    if( (nplant * (balive + bdead)) > min_recruit_size)then
                       
                       inew = inew + 1
                       
                       recruit_array(inew,1) = real(pft)
                       recruit_array(inew,2) = nplant
                       recruit_array(inew,3) = hite
                       recruit_array(inew,4) = dbh
                       recruit_array(inew,5) = bdead
                       recruit_array(inew,6) = bleaf
                       recruit_array(inew,7) = balive
                       recruit_array(inew,8) = bleaf * nplant * sla(pft)
                       recruit_array(inew,9) = csite%can_temp(ipa)
                       
                    endif
                 else
                    
                    ! If we have reached this branch, we are in an 
                    ! agricultural patch.  Send the seed litter to the soil
                    ! pools for decomposition.
                    csite%fast_soil_N(ipa) = csite%fast_soil_N(ipa)   &
                         + csite%repro(pft,ipa) / c2n_recruit(pft)
                    csite%fast_soil_C(ipa) = csite%fast_soil_C(ipa) +   &
                         csite%repro(pft,ipa)
                    csite%repro(pft,ipa) = 0.0
                    
                 endif
              endif
              
              ! Reset the carbon available for reproduction.
              
              csite%repro(pft,ipa) = 0.0
              
           enddo

           ncohorts_new = cpatch%ncohorts + inew
           
           ! The number of recruits is now known. Allocate the temporary
           ! patch vector with the current number plus the number of recruits
           
           call allocate_patchtype(temppatch,cpatch%ncohorts)

           ! Fill the temp space with the current patches
           
           call copy_patchtype(cpatch,temppatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
           
           ! Deallocate the current patch
           
           if (cpatch%ncohorts > 0) call deallocate_patchtype(cpatch)

           ! Reallocate the current site

           call allocate_patchtype(cpatch,ncohorts_new)
           
           ! Transfer the temp values back in
           
           call copy_patchtype(temppatch,cpatch,1,temppatch%ncohorts,1,temppatch%ncohorts)
           
           inew = 0
           do ico = temppatch%ncohorts+1,ncohorts_new
              inew = inew + 1
              cpatch%pft(ico)       = int(recruit_array(inew,1))
              cpatch%nplant (ico)   = recruit_array(inew,2)
              cpatch%hite(ico)      = recruit_array(inew,3)
              cpatch%dbh(ico)       = recruit_array(inew,4)
              cpatch%bdead(ico)     = recruit_array(inew,5)
              cpatch%bleaf(ico)     = recruit_array(inew,6)
              cpatch%phenology_status(ico) = 0
              cpatch%balive(ico)    = recruit_array(inew,7)
              cpatch%lai(ico)       = recruit_array(inew,8)
              cpatch%bstorage(ico)  = 0.0
              cpatch%veg_temp(ico)  = recruit_array(inew,9)
              cpatch%veg_water(ico) = 0.0

              !----- Because we assigned no water, the internal energy 
              !      is simply hcapveg*T
              
              cpatch%hcapveg(ico) = calc_hcapveg(cpatch%bleaf(ico),cpatch%nplant(ico) &
                                                ,cpatch%lai(ico),cpatch%pft(ico)      &
                                                ,cpatch%phenology_status(ico))
              cpatch%veg_energy(ico) = cpatch%hcapveg(ico) * cpatch%veg_temp(ico)
              
              ! Setting new_recruit_flag to 1 indicates that 
              ! this cohort is included when we tally agb_recruit,
              ! basal_area_recruit.
              cpatch%new_recruit_flag(ico) = 1

              ! Carry out standard initialization.
              call init_ed_cohort_vars_array(cpatch,ico,cpoly%lsl(isi))
              
              csite%cohort_count(ipa) = csite%cohort_count(ipa) + 1
                            
           enddo

           ! Remove the temporary patch
           
           if (temppatch%ncohorts > 0) call deallocate_patchtype(temppatch)
           
        enddo
        
        
        !----------------------------
        ! Re-sort, terminate, fuse, and split.
        !----------------------------
        
        do ipa = 1,csite%npatches
           cpatch => csite%patch(ipa)

           if(cpatch%ncohorts>0 .and. maxcohort >= 0) then

              call terminate_cohorts_ar(csite,ipa)
              call fuse_cohorts_ar(csite,ipa, cpoly%green_leaf_factor(:,isi), cpoly%lsl(isi))                         
              call split_cohorts_ar(cpatch, cpoly%green_leaf_factor(:,isi), cpoly%lsl(isi))
              
           endif

           csite%cohort_count(ipa) = cpatch%ncohorts !THIS IS REDUNDANT

           call update_patch_derived_props_ar(csite, cpoly%lsl(isi), cpoly%met(isi)%rhos, ipa)

        enddo

        ! Update site properties.
        call update_site_derived_props_ar(cpoly, 0,isi)
        
        cpoly%min_monthly_temp(isi) = huge(1.)
        
     enddo
     
  enddo


  deallocate(temppatch)
  
  first_time = .false.

  return
end subroutine reproduction_ar
