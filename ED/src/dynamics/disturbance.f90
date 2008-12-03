module disturbance_utils_ar

  use ed_state_vars,only: &
       allocate_patchtype,   &
       copy_patchtype,       &
       deallocate_patchtype, &
       allocate_sitetype,    &
       deallocate_sitetype,  &
       copy_sitetype_mask

  use fuse_fiss_utils_ar,only: &
       fuse_cohorts_ar,      &
       terminate_cohorts_ar, &
       split_cohorts_ar

  contains

  subroutine apply_disturbances_ar(cgrid)

    use ed_state_vars,only: edtype,polygontype,sitetype,patchtype
    use misc_coms, only: current_time
    use disturb_coms, only: treefall_age_threshold,  &
         min_new_patch_area
    use max_dims, only: n_dist_types, n_pft, n_dbh

    implicit none

    type(edtype),target       :: cgrid
    type(polygontype),pointer :: cpoly
    type(sitetype),pointer    :: csite
    type(sitetype),pointer    :: tsite
    type(patchtype),pointer   :: qpatch 
    integer :: ipy,isi,ipa
    
    integer :: q
    real :: area

    real :: dA
    real :: area_fac
    real, dimension(n_pft, n_dbh) :: initial_agb
    real, dimension(n_pft, n_dbh) :: initial_basal_area

    integer,allocatable :: disturb_mask(:)
    integer :: onsp

    allocate(tsite)  ! The disturbance site
    
    do ipy = 1,cgrid%npolygons
       
     cpoly => cgrid%polygon(ipy)

     
     do isi = 1,cpoly%nsites
        
        csite => cpoly%site(isi)

        
        ! Store AGB, basal area profiles in memory.
        call update_site_derived_props_ar(cpoly, 1,isi)
        initial_agb(1:n_pft, 1:n_dbh) = cpoly%agb(1:n_pft, 1:n_dbh,isi)
        initial_basal_area(1:n_pft, 1:n_dbh) = cpoly%basal_area(1:n_pft, 1:n_dbh,isi)
        
        ! First take care of harvesting: secondary -> secondary and 
        ! primary -> secondary.
        call apply_forestry_ar(cpoly,isi, current_time%year, cpoly%met(isi)%rhos)
        
        ! Update the cut output variables
        call update_site_derived_props_ar(cpoly, 1,isi)
        
        cpoly%agb_cut(1:n_pft, 1:n_dbh,isi) = cpoly%agb_cut(1:n_pft, 1:n_dbh,isi) +  &
             initial_agb(1:n_pft, 1:n_dbh) - cpoly%agb(1:n_pft, 1:n_dbh,isi)
        
        cpoly%basal_area_cut(1:n_pft, 1:n_dbh,isi) =   &
             cpoly%basal_area_cut(1:n_pft, 1:n_dbh,isi) + &
             initial_basal_area(1:n_pft, 1:n_dbh) - cpoly%basal_area(1:n_pft, 1:n_dbh,isi)
        
        
        onsp = csite%npatches    ! Original Number (of) Site Patches


        ! Create a temporary site with vectors containing all current patches
        ! as well as n_dist_types patches.
        ! Create the newly disturbed patches in here, and depending on
        ! How many are created, repopulate the existing site's patch vectors
        
        call allocate_sitetype(tsite,onsp)

        allocate(disturb_mask(onsp + n_dist_types))
        disturb_mask = 0
        disturb_mask(1:onsp) = 1

        ! Transfer the origial patch values into the front end of the temp's space
        call copy_sitetype_mask(csite,tsite,disturb_mask(1:onsp), &
             sum(disturb_mask),sum(disturb_mask))

        ! Reallocate and transfer them back
        
        call deallocate_sitetype(csite)
        
        call allocate_sitetype(csite,onsp + n_dist_types)
        
        call copy_sitetype_mask(tsite,csite,disturb_mask(1:onsp), &
             sum(disturb_mask),sum(disturb_mask))
       
        call deallocate_sitetype(tsite)


        ! Initialize all the potential as well as implemented disturbance
        ! patches. Since agb and basal area is calculated based on these
        ! numbers
        do q = onsp+1, onsp+n_dist_types
            call initialize_disturbed_patch_ar(csite,cpoly%met(isi)%atm_tmp,cpoly%met(isi)%atm_shv,q,1)
        enddo

        ! Loop over q, the *destination* landuse type
        do q = 1, n_dist_types
           
           ! First, decide if enough area is disturbed to warrant creating
           ! a new patch.
           area = 0.0

           do ipa=1,onsp

              if( (q == 1 .and. csite%dist_type(ipa) > 1) .or.  & ! conversion to ag
                   (q == 2 .and. csite%dist_type(ipa) == 1) .or.  & ! abandonment
                   (q == 3 .and. csite%dist_type(ipa) > 1 .and.   & ! natural.....
                   (csite%age(ipa) > treefall_age_threshold   &  ! old enough for treefall
                   .or. cpoly%nat_dist_type(isi) == 1)) )then  ! or it's a fire
                 
                 area = area + csite%area(ipa) * (1.0 - exp(   &
                      - (cpoly%disturbance_rates(q,csite%dist_type(ipa),isi)  &
                      + cpoly%disturbance_memory(q,csite%dist_type(ipa),isi))))
                 
              endif
              
           enddo

           
           
           if(area > min_new_patch_area)then
              write(unit=*,fmt='(a,1x,es12.5,1x,a,1x,i5)') &
                  ' ---> Making new patch, with area=',area,' for dist_type=',q

              ! Set the flag that this patch should be kept as a newly created
              ! transition patch.
              disturb_mask(onsp+q) = 1   
 
              csite%dist_type(onsp+q)  = q
              csite%plantation(onsp+q) = 0
              csite%area(onsp+q)       = area
              
              ! Initialize to zero the new trasitioned patches (redundant)
              call initialize_disturbed_patch_ar(csite,cpoly%met(isi)%atm_tmp,cpoly%met(isi)%atm_shv,onsp+q,1)
              
              ! Now go through patches, adding its contribution to the new patch.
              do ipa=1,onsp
                                 
                 if( (q == 1 .and. csite%dist_type(ipa) > 1) .or.  & ! conversion to ag
                      (q == 2 .and. csite%dist_type(ipa) == 1) .or.  & ! abandonment
                      (q == 3 .and. csite%dist_type(ipa) > 1 .and.   & ! natural.....
                      (csite%age(ipa) > treefall_age_threshold   &  ! old enough for treefall
                      .or. cpoly%nat_dist_type(isi) == 1)) )then  ! or it's a fire
                
                    dA = csite%area(ipa) * (1.0 - exp(   &
                         - (cpoly%disturbance_rates(q,csite%dist_type(ipa),isi)  &
                         + cpoly%disturbance_memory(q,csite%dist_type(ipa),isi))))
                    
                    area_fac = dA / csite%area(onsp+q)
                    
                    call increment_patch_vars_ar(csite,q+onsp,ipa,area_fac)
                    call insert_survivors_ar(csite,q+onsp,ipa,q,area_fac,cpoly%nat_dist_type(isi))

                    call accum_dist_litt_ar(csite,q+onsp,ipa,q,area_fac, &
                         cpoly%loss_fraction(q,isi),cpoly%nat_dist_type(isi))
                    
                    ! update patch area
                    csite%area(ipa) = csite%area(ipa) - dA
                    
                 endif
              enddo
              
              ! if the new patch is agriculture, plant it with grasses.
              if(q == 1)call plant_patch_ar(csite,q+onsp,cpoly%agri_stocking_pft(isi),   &
                   cpoly%agri_stocking_density(isi), 1.0, cpoly%lsl(isi))
              
              qpatch => csite%patch(q+onsp)

              !  fuse then terminate cohorts
              if (csite%patch(q+onsp)%ncohorts > 0) then
                 call fuse_cohorts_ar(csite,q+onsp, cpoly%green_leaf_factor(:,isi), cpoly%lsl(isi))
                 call terminate_cohorts_ar(csite,q+onsp)
                 call split_cohorts_ar(qpatch, cpoly%green_leaf_factor(:,isi), cpoly%lsl(isi))
              endif
          
              ! Store AGB, basal area profiles in memory.
              initial_agb(1:n_pft, 1:n_dbh) = cpoly%agb(1:n_pft, 1:n_dbh,isi)
              initial_basal_area(1:n_pft, 1:n_dbh) = cpoly%basal_area(1:n_pft, 1:n_dbh, isi)
              
              ! Update the derived properties including veg_height, lai
              call update_patch_derived_props_ar(csite, cpoly%lsl(isi), cpoly%met(isi)%rhos,q+onsp)

              ! Update soil temp, fracliq, etc.
              call new_patch_sfc_props_ar(csite, q+onsp, cpoly%met(isi)%rhos)
              
              ! Update AGB, basal area.
              ! !!!!!!!!!! SHOULD THIS BE HERE OR OUTSIDE THIS LOOP?? !!!!!!!!
              ! !!!! IS THIS BECAUSE UPDATED SITE PROPS ARE REQUIRED ON EVERY
              ! !!!! NEW DISTURBANCE PATCH CREATION ATTEMPT? I WILL LEAVE IT
              ! !!!! FOR NOW.  IT DOES NOT INCREMEMNT ANY VARIABLES SO IT CANT HURT

              call update_site_derived_props_ar(cpoly,1,isi)
              
              ! Update either cut or mortality
              if(q /= 3)then
                 cpoly%agb_cut(1:n_pft, 1:n_dbh,isi) = cpoly%agb_cut(1:n_pft, 1:n_dbh,isi) +  &
                      initial_agb(1:n_pft, 1:n_dbh) - cpoly%agb(1:n_pft, 1:n_dbh,isi)
                 cpoly%basal_area_cut(1:n_pft, 1:n_dbh,isi) =   &
                      cpoly%basal_area_cut(1:n_pft, 1:n_dbh,isi) +  &
                      initial_basal_area(1:n_pft, 1:n_dbh) -   &
                      cpoly%basal_area(1:n_pft, 1:n_dbh,isi)
              else
                 cpoly%agb_mort(1:n_pft, 1:n_dbh,isi) = cpoly%agb_mort(1:n_pft, 1:n_dbh,isi) +  &
                      initial_agb(1:n_pft, 1:n_dbh) - cpoly%agb(1:n_pft, 1:n_dbh,isi)
                 cpoly%basal_area_mort(1:n_pft, 1:n_dbh,isi) =   &
                      cpoly%basal_area_mort(1:n_pft, 1:n_dbh,isi) +  &
                      initial_basal_area(1:n_pft, 1:n_dbh) -   &
                      cpoly%basal_area(1:n_pft, 1:n_dbh,isi)
              endif
              
              !  clear the disturbance memory for this disturbance type
              cpoly%disturbance_memory(q,1:n_dist_types,isi) = 0.0
              
           else
              if(area > 0.0)then
                 ! the patch creation has been skipped because 
                 !    the area was too small 
                 ! put the current disturbance rates in memory to be 
                 !    added at the next timestep
                 cpoly%disturbance_memory(q,1:n_dist_types,isi) =   &
                      cpoly%disturbance_memory(q,1:n_dist_types,isi) +   &
                      cpoly%disturbance_rates(q,1:n_dist_types,isi)
              endif
           endif
           
        enddo

        ! Reallocate the current site to fit the original patches and whatever
        ! was generated in disturbance (ie, make it non sparse)
        ! Populate the original site with both the modified 
        ! original patches, and the newly created patches.  The index of all of these
        ! are disturb_mask.  This mask should be ones for all original patches, and 
        ! sparse from there after.
        ! --------------------------------------------------------------------------
        call allocate_sitetype(tsite,sum(disturb_mask))
        
        call copy_sitetype_mask(csite,tsite,disturb_mask,size(disturb_mask),sum(disturb_mask))
        
        call deallocate_sitetype(csite)
        
        call allocate_sitetype(csite,sum(disturb_mask))
        
        disturb_mask = 0
        disturb_mask(1:csite%npatches) = 1
        call copy_sitetype_mask(tsite,csite,disturb_mask(1:csite%npatches),sum(disturb_mask),sum(disturb_mask))
        
        call deallocate_sitetype(tsite)
        deallocate(disturb_mask)
        ! --------------------------------------------------------------------------

        
     enddo
  enddo

  deallocate(tsite)
  
  return
end subroutine apply_disturbances_ar

  !====================================================================

  subroutine site_disturbance_rates_ar(month, year, cgrid)

    use ed_state_vars,only:edtype,polygontype,sitetype
    use disturb_coms, only: treefall_disturbance_rate, lutime, include_fire
    use pft_coms, only: agf_bs

    implicit none

    type(edtype),target       :: cgrid
    type(polygontype),pointer :: cpoly
    type(sitetype),pointer    :: csite
    integer :: ipy,isi

    integer, intent(in) :: month
    integer, intent(in) :: year

    integer :: im,iyear,useyear
    type(lutime), pointer :: clutime
    real :: fire_dist_rate

    ! Loop over sites and polygons.
    
    do ipy = 1,cgrid%npolygons
       
       cpoly => cgrid%polygon(ipy)
       
       do isi = 1,cpoly%nsites
          
          csite => cpoly%site(isi)

          !  Calculate fire disturbance rates only if fire is on.
          if(include_fire == 1)then
             fire_dist_rate = sum(cpoly%lambda_fire(1:12,isi)) / 12.0
          else
             fire_dist_rate = 0.0
          endif

          cpoly%fire_disturbance_rate = fire_dist_rate
    
          !  treefall disturbance is currently spatiotemporally constant, 
          !  from ED2IN

          ! For natural disturbance, use largest disturbance mode.
          if(fire_dist_rate > treefall_disturbance_rate)then
             cpoly%nat_disturbance_rate(isi) = fire_dist_rate
             cpoly%nat_dist_type(isi) = 1
          else
             cpoly%nat_disturbance_rate(isi) = treefall_disturbance_rate
             cpoly%nat_dist_type(isi) = 0
          endif
    
          ! Set disturbance rates assuming only natural disturbance
          cpoly%disturbance_rates(1:2,1:3,isi)= 0.0
          cpoly%disturbance_rates(3,1,isi)= 0.0
          cpoly%disturbance_rates(3,2:3,isi)= cpoly%nat_disturbance_rate(isi)
          
          ! Now it is time for anthropogenic disturbance rates.
             
          ! Only consider years after the start of the dataset.
          ! For years after the end of the dataset, repeat the last year.

          if (cpoly%num_landuse_years(isi)>0) then

             if(year >= cpoly%clutimes(1,isi)%landuse_year)then

                useyear = cpoly%num_landuse_years(isi)
                ! Loop over years
                find_lu_year: do iyear = 1,cpoly%num_landuse_years(isi)
                   
                   if (year.eq.cpoly%clutimes(iyear,isi)%landuse_year) then
                      useyear = iyear
                      exit find_lu_year
                   endif
                   
                enddo find_lu_year

                !    update land-use transition matrix
                
                clutime => cpoly%clutimes(useyear,isi)

                ! AA2AA, agriculture to agriculture
                cpoly%disturbance_rates(1,1,isi)= 0.0
                
                ! NA2AA, secondary forest to agriculture
                cpoly%disturbance_rates(1,2,isi)= clutime%landuse(7) +   &
                     clutime%landuse(9)
                
                ! NN2AA, primary forest to agriculture
                cpoly%disturbance_rates(1,3,isi)= clutime%landuse(4) +   &
                     clutime%landuse(5)
                
                ! AA2NA, agriculture to secondary forest
                cpoly%disturbance_rates(2,1,isi)= clutime%landuse(8) +   &
                     clutime%landuse(10)
                
                ! NA2NA, secondary forest to secondary forest (this is taken 
                ! care of in the harvesting.)
                cpoly%disturbance_rates(2,2,isi)= 0.0
                
                ! NN2NA, primary forest to secondary forest (zero here 
                ! because we meet a biomass target instead.)
                cpoly%disturbance_rates(2,3,isi)= 0.0 ! clutime%landuse(11)
                
                ! AA2NN, agriculture to primary forest (should be zero)
                cpoly%disturbance_rates(3,1,isi)= 0.0 
                !cpoly%disturbance_rates(3,1,isi) = clutime%landuse(3) + clutime%landuse(6)
                
                ! NA2NN, secondary forest to primary forest
                cpoly%disturbance_rates(3,2,isi)= cpoly%nat_disturbance_rate(isi)
                
                ! NN2NN, primary forest to primary forest
                cpoly%disturbance_rates(3,3,isi)= cpoly%nat_disturbance_rate(isi)

             endif ! after first year
          endif ! landuse exists

          ! fraction of above ground litter from disturbance that is 
          ! removed from patch
          cpoly%loss_fraction(1,isi) = agf_bs
          cpoly%loss_fraction(2,isi) = agf_bs
          cpoly%loss_fraction(3,isi) = 0.0

       enddo
       
    enddo

    return
  end subroutine site_disturbance_rates_ar

  !=====================================================================

  subroutine initialize_disturbed_patch_ar(csite,atm_tmp,atm_shv,np,dp)
    
    use ed_state_vars, only: sitetype,patchtype
    use consts_coms, only : t00
    use grid_coms, only: nzs, nzg

    implicit none
    type(sitetype),target    :: csite
    type(patchtype),pointer  :: cpatch
    real, intent(in)         :: atm_tmp,atm_shv
    integer, intent(in)      :: np,dp
    integer                  :: ipa,k


    ! SHOULD WE CONSIDER USING A DONOR PATCH
    ! TO INITIALIZE MORE OF THE TRANSIENT STATE VARIABLES
    ! PRIOR TO THIS MESSAGE, NTEXT_SOIL WAS NOT BEING SET PROPERLY


    csite%patch(np)%ncohorts = 0

    ! Initializing some water/energy properties 
    csite%can_temp(np) = atm_tmp
    csite%can_shv(np)  = atm_shv
    ! For now, choose heat/vapor capacities for stability
    csite%can_depth(np) = 30.0
    do k=1,nzs
       csite%sfcwater_tempk(k,np) = t00   ! Set canopy temp to 0 C
       csite%sfcwater_fracliq(k,np) = 1.0 ! Set to 100% liquid
    end do


    ! Included the initialization for all fast variables.
    call init_ed_patch_vars_array(csite,np,np)

    ! dist_type is not set here.

    csite%age(np) = 0.0

    ! area is not set here.

    csite%fast_soil_C(np) = 0.0

    csite%slow_soil_C(np) = 0.0

    csite%structural_soil_C(np) = 0.0

    csite%structural_soil_L(np) = 0.0

    csite%mineralized_soil_N(np) = 0.0

    csite%fast_soil_N(np) = 0.0

    csite%sum_dgd(np) = 0.0

    csite%sum_chd(np) = 0.0

    ! plantation is not set here.

    csite%can_depth(np) = 0.0


    !--------------------------------------------------------------------------------------!
    !     I think the following variables need to receive the donor patch value, otherwise !
    ! it will try to use ed_grndvap with non-sense values.                                 !
    !--------------------------------------------------------------------------------------!
    csite%ntext_soil(1:nzg,np) = csite%ntext_soil(1:nzg,dp)

    csite%can_temp(np) = 0.0

    csite%can_shv(np) = 0.0

    csite%soil_energy(1:nzg,np) = 0.0

    csite%soil_water(1:nzg,np) = 0.0

    csite%sfcwater_mass(1:nzs,np) =  0.0

    csite%sfcwater_energy(1:nzs,np) = 0.0

    csite%sfcwater_depth(1:nzs,np) = 0.0

    !--------------------------------------------------------------------------------------!
    csite%rough(np) = 0.0

    ! this 
    call init_ed_patch_vars_array(csite,np,np)

    csite%fsc_in(np) = 0.0

    csite%ssc_in(np) = 0.0

    csite%ssl_in(np) = 0.0

    csite%fsn_in(np) = 0.0

    csite%total_plant_nitrogen_uptake(np) = 0.0

    return
  end subroutine initialize_disturbed_patch_ar

  !======================================================================

  subroutine increment_patch_vars_ar(csite,np, cp, area_fac)

    use ed_state_vars, only: sitetype,patchtype
    use max_dims, only: n_pft
    use grid_coms, only: nzg

    implicit none
    type(sitetype),target    :: csite
    type(patchtype),pointer  :: cpatch
    integer :: np,cp

    real :: area_fac
    integer :: k

    csite%fast_soil_C(np) = csite%fast_soil_C(np) + csite%fast_soil_C(cp) * area_fac

    csite%slow_soil_C(np) = csite%slow_soil_C(np) + csite%slow_soil_C(cp) * area_fac

    csite%structural_soil_C(np) = csite%structural_soil_C(np) + csite%structural_soil_C(cp) *   &
         area_fac

    csite%structural_soil_L(np) = csite%structural_soil_L(np) + csite%structural_soil_L(cp) *   &
         area_fac

    csite%mineralized_soil_N(np) = csite%mineralized_soil_N(np) + csite%mineralized_soil_N(cp) *   &
         area_fac

    csite%fast_soil_N(np) = csite%fast_soil_N(np) + csite%fast_soil_N(cp) * area_fac

    csite%sum_dgd(np) = csite%sum_dgd(np) + csite%sum_dgd(cp) * area_fac

    csite%sum_chd(np) = csite%sum_chd(np) + csite%sum_chd(cp) * area_fac

    csite%can_temp(np) = csite%can_temp(np) + csite%can_temp(cp) * area_fac

    csite%can_shv(np) = csite%can_shv(np) + csite%can_shv(cp) * area_fac

    csite%can_depth(np) = csite%can_depth(np) + csite%can_depth(cp) * area_fac

    do k = 1, csite%nlev_sfcwater(cp)
       csite%sfcwater_mass(k,np) = csite%sfcwater_mass(k,np) + csite%sfcwater_mass(k,cp) *   &
            area_fac

       !!!! IS THE FOLLOWING LINE CORRECT? IT SEEMS TO BE ADDING AN EXTRA AND UNECESARY
       !!!! MASS TERM FROM THE DONOR PATCH  !!!

       csite%sfcwater_energy(k,np) = csite%sfcwater_energy(k,np) +   &
            csite%sfcwater_energy(k,np) * csite%sfcwater_mass(k,cp) * area_fac
       csite%sfcwater_depth(k,np) = csite%sfcwater_depth(k,np) + csite%sfcwater_depth(k,cp) *   &
            area_fac
    enddo

    do k = 1, nzg
       csite%ntext_soil(k,np) = csite%ntext_soil(k,np)
       csite%soil_energy(k,np) = csite%soil_energy(k,np) +csite%soil_energy(k,cp) * area_fac
       csite%soil_water(k,np) = csite%soil_water(k,np) + csite%soil_water(k,cp) * area_fac
    enddo

    csite%rough(np) = csite%rough(np) + csite%rough(cp) * area_fac

    csite%mean_rh(np) = csite%mean_rh(np) + csite%mean_rh(cp) * area_fac

    csite%dmean_A_decomp(np) = csite%dmean_A_decomp(np) + csite%dmean_A_decomp(cp) * area_fac

    csite%dmean_Af_decomp(np) = csite%dmean_Af_decomp(np) + csite%dmean_Af_decomp(cp) * area_fac

    csite%repro(1:n_pft,np) = csite%repro(1:n_pft,np) + csite%repro(1:n_pft,cp) * area_fac

    csite%fsc_in(np) = csite%fsc_in(np) + csite%fsc_in(cp) * area_fac

    csite%ssc_in(np) = csite%ssc_in(np) + csite%ssc_in(cp) * area_fac

    csite%ssl_in(np) = csite%ssl_in(np) + csite%ssl_in(cp) * area_fac

    csite%fsn_in(np) = csite%fsn_in(np) + csite%fsn_in(cp) * area_fac

    csite%total_plant_nitrogen_uptake(np) = csite%total_plant_nitrogen_uptake(np) +   &
         csite%total_plant_nitrogen_uptake(cp) * area_fac

    return
  end subroutine increment_patch_vars_ar

  !=======================================================================

  subroutine insert_survivors_ar(csite, np, cp, q, area_fac,nat_dist_type)

    use ed_state_vars, only: sitetype,patchtype
    use therm_lib,only: update_veg_energy_ct
    
    implicit none
    type(sitetype),target    :: csite
    type(patchtype),pointer  :: cpatch !The old patches we are looping
    type(patchtype),pointer  :: npatch !New transition patch
    type(patchtype),pointer  :: tpatch !Temp patch
    integer :: np,cp,ico,nco
    
    integer, intent(in) :: q
    integer, intent(in) :: nat_dist_type
    real, intent(in) :: area_fac
    real :: n_survivors

    real :: cohort_area_fac

    integer,allocatable :: mask(:)

    cpatch => csite%patch(cp)
    npatch => csite%patch(np)
    
    allocate(tpatch)

    allocate(mask(cpatch%ncohorts))
    mask = 0
    
    do ico = 1,cpatch%ncohorts
       
       cohort_area_fac = survivorship_ar(q,nat_dist_type, csite, cp, ico) * area_fac
       n_survivors = cpatch%nplant(ico) * cohort_area_fac
       
       ! If something survived, make a new cohort
       if(n_survivors > 0.0) mask(ico) = 1
       
    enddo

    if (npatch%ncohorts > 0) then
       
       ! If the new patch has received survivors from a donor
       ! already, then it should have cohorts. So the temp patch
       ! vector will be the sum of the new cohorts found here, plus
       ! those already applied previously in the loop calling this
       
       nco = npatch%ncohorts
       call allocate_patchtype(tpatch,sum(mask) + npatch%ncohorts)
       call copy_patchtype(npatch,tpatch,1,npatch%ncohorts,1,npatch%ncohorts)
       call deallocate_patchtype(npatch)
       
    else
       
       nco = 0
       call allocate_patchtype(tpatch,sum(mask))

    endif


    do ico = 1,cpatch%ncohorts
       
       cohort_area_fac = survivorship_ar(q,nat_dist_type, csite, cp, ico) * area_fac
       n_survivors = cpatch%nplant(ico) * cohort_area_fac
       
       if(mask(ico) > 0.0) then

          nco = nco + 1
          
          call copy_patchtype(cpatch,tpatch,ico,ico,nco,nco)
          
          ! Adjust area-based variables
          tpatch%nplant(nco)     = tpatch%nplant(nco)    * cohort_area_fac
          tpatch%lai(nco)        = tpatch%lai(nco)       * cohort_area_fac       
          tpatch%veg_water(nco)  = tpatch%veg_water(nco) * cohort_area_fac
          tpatch%mean_gpp(nco)   = tpatch%mean_gpp(nco)  * cohort_area_fac
          tpatch%mean_leaf_resp(nco)      = tpatch%mean_leaf_resp(nco) * cohort_area_fac
          tpatch%mean_root_resp(nco)      = tpatch%mean_root_resp(nco) * cohort_area_fac
          tpatch%growth_respiration(nco)  = tpatch%growth_respiration(nco) * cohort_area_fac
          tpatch%storage_respiration(nco) = tpatch%storage_respiration(nco) * cohort_area_fac
          tpatch%vleaf_respiration(nco)   = tpatch%vleaf_respiration(nco) * cohort_area_fac
          tpatch%Psi_open(nco)            = tpatch%Psi_open(nco) * cohort_area_fac

          ! Adjust the vegetation energy - RGK 11-26-2008

          call update_veg_energy_ct(cpatch,ico)

          
       endif
          
    enddo  ! end loop over cohorts

    ! Copy the temporary patch into the newpatch
    
    call allocate_patchtype(npatch,tpatch%ncohorts)
    call copy_patchtype(tpatch,npatch,1,tpatch%ncohorts,1,tpatch%ncohorts)
    call deallocate_patchtype(tpatch)
    deallocate(tpatch)
    
    deallocate(mask)


    return
  end subroutine insert_survivors_ar

  !======================================================================
  subroutine accum_dist_litt_ar(csite,np,cp, q, area_fac,loss_fraction,nat_dist_type)
    
    use ed_state_vars, only: sitetype,patchtype,polygontype
    use decomp_coms, only: f_labile
    use max_dims, only: n_pft
    use pft_coms, only: c2n_storage, c2n_leaf, c2n_recruit, c2n_stem, l2n_stem
    use grid_coms, only: nzg

    implicit none

    type(sitetype),target    :: csite
    type(patchtype),pointer  :: cpatch
    type(patchtype),pointer  :: npatch
    integer :: np,cp,ico,nco,isi
    real,intent(in) :: loss_fraction

    integer, intent(in) :: q
    real, intent(in) :: area_fac
    integer,intent(in) :: nat_dist_type
    real :: fast_litter
    real :: struct_litter
    real :: fast_litter_n

    fast_litter = 0.0
    struct_litter = 0.0
    fast_litter_n = 0.0

    cpatch => csite%patch(cp)
    npatch => csite%patch(np)

    do ico = 1,cpatch%ncohorts

       fast_litter = fast_litter + (1.0 - &
            survivorship_ar(q,nat_dist_type,csite, cp, ico)) * &
            (f_labile(cpatch%pft(ico)) * cpatch%balive(ico) &
            + cpatch%bstorage(ico)) * cpatch%nplant(ico)

       fast_litter_n = fast_litter_n + (1.0 - &
            survivorship_ar(q,nat_dist_type, csite, cp, ico)) * &
            (f_labile(cpatch%pft(ico)) * cpatch%balive(ico) &
            / c2n_leaf(cpatch%pft(ico)) + cpatch%bstorage(ico) /  &
            c2n_storage ) * cpatch%nplant(ico)

       struct_litter = struct_litter + cpatch%nplant(ico) *   &
            (1.0 - survivorship_ar(q,nat_dist_type, csite, cp, ico)) * ( (1.0 -   &
            loss_fraction ) * cpatch%bdead(ico) +   & ! DOUBLE CHECK THIS, IS THIS RIGHT??
            (1.0 - f_labile(cpatch%pft(ico))) * cpatch%balive(ico))
       
    enddo

    !  Load disturbance litter directly into carbon and N pools
    csite%fast_soil_C(np) = csite%fast_soil_C(np) + fast_litter * area_fac

    csite%structural_soil_C(np) = csite%structural_soil_C(np) + struct_litter * area_fac

    csite%structural_soil_L(np) = csite%structural_soil_L(np) + l2n_stem / c2n_stem *   &
         struct_litter * area_fac

    csite%fast_soil_N(np) = csite%fast_soil_N(np) + fast_litter_n * area_fac

    return
  end subroutine accum_dist_litt_ar

  !============================================================================

  real function survivorship_ar(dest_type,poly_dest_type,csite, ipa, ico)

    use ed_state_vars, only: patchtype,sitetype
    use disturb_coms, only: treefall_hite_threshold
    use pft_coms, only: treefall_s_ltht, treefall_s_gtht
    
    implicit none

    type(patchtype),pointer  :: cpatch
    type(sitetype),target   :: csite
    integer :: ico,ipa
    integer, intent(in) :: dest_type
    integer, intent(in) :: poly_dest_type

    cpatch => csite%patch(ipa)

    if(dest_type == 1)then
       survivorship_ar = 0.0  !  agric
    else
       if(dest_type == 2)then
          survivorship_ar = 0.0   ! secondary land
       else  !   natural land 
          if(dest_type == 3)then
             if(poly_dest_type == 1)then
                survivorship_ar = 0.0 ! fire
             elseif(poly_dest_type == 0)then
                if(cpatch%hite(ico) < treefall_hite_threshold)then
                   survivorship_ar =  treefall_s_ltht(cpatch%pft(ico)) 
                else 
                   survivorship_ar = treefall_s_gtht(cpatch%pft(ico))
                endif
             endif
          endif
       endif
    endif
    
    return
  end function survivorship_ar

  !=============================================================

  subroutine plant_patch_ar(csite, np, pft, density, height_factor, lsl)

    use ed_state_vars,only : sitetype,patchtype

    use pft_coms, only: q, qsw, sla, hgt_min
    use misc_coms, only: dtlsm
    use fuse_fiss_utils_ar, only : sort_cohorts_ar
!    use canopy_air_coms, only: hcapveg_ref,heathite_min
    use therm_lib,only : calc_hcapveg
    use consts_coms, only: t3ple

    implicit none

    type(patchtype),pointer  :: cpatch
    type(patchtype),pointer  :: tpatch
    type(sitetype),target   :: csite

    integer, intent(in) :: lsl
    integer :: np,nc
    integer, intent(in) :: pft
    real, intent(in) :: density
    real, intent(in) :: height_factor
    real :: h2dbh
    real :: dbh2bd
    real :: dbh2bl
    real :: hcapveg

    
    cpatch => csite%patch(np)

    ! Reallocate the current cohort with an extra space for the
    ! planted plantation cohort.  It is unlikely but perhaps possible
    ! For this patch to have no cohorts yet.
    ! ------------------------------------------------------------
    nc = cpatch%ncohorts + 1
    if ( cpatch%ncohorts>0) then
       
       allocate(tpatch)
       call allocate_patchtype(tpatch,cpatch%ncohorts)
       call copy_patchtype(cpatch,tpatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
       call deallocate_patchtype(cpatch)
       call allocate_patchtype(cpatch,tpatch%ncohorts + 1)
       call copy_patchtype(tpatch,cpatch,1,tpatch%ncohorts,1,tpatch%ncohorts)
       call deallocate_patchtype(tpatch)
       deallocate(tpatch)
    else
       call allocate_patchtype(cpatch,1)
    endif
       
    cpatch%ncohorts = nc
    csite%paco_n(np)= nc

    ! Just make one cohort.  It will soon be split at the splitting call of
    ! apply_disturbances().  Place this cohort in the insertion point.


    cpatch%pft(nc) = pft
    cpatch%nplant(nc) = density
    cpatch%hite(nc) = hgt_min(cpatch%pft(nc)) * height_factor
    cpatch%dbh(nc) = h2dbh(cpatch%hite,cpatch%pft)
    cpatch%bdead(nc) = dbh2bd(cpatch%dbh(nc),cpatch%hite(nc),cpatch%pft(nc))
    cpatch%bleaf(nc) = dbh2bl(cpatch%dbh(nc),cpatch%pft(nc))
    print*,"add",cpatch%hite(nc),cpatch%dbh(nc),cpatch%bdead(nc),cpatch%bleaf(nc)
    cpatch%phenology_status = 0
    cpatch%balive(nc) = cpatch%bleaf(nc) * &
         (1.0 + q(cpatch%pft(nc)) + qsw(cpatch%pft(nc)) * cpatch%hite(nc))
    cpatch%lai(nc) = cpatch%bleaf(nc) * cpatch%nplant(nc) * sla(cpatch%pft(nc))
    cpatch%bstorage(nc) = 0.0

    cpatch%veg_temp(nc) = csite%can_temp(np)
    cpatch%veg_water(nc) = 0.0

    !----- Because we assigned no water, the internal energy is simply hcapveg*(T-T3)

    cpatch%hcapveg(nc) = calc_hcapveg(cpatch%bleaf(nc),cpatch%bdead(nc), &
         cpatch%nplant(nc),cpatch%pft(nc))
    
    cpatch%veg_energy(nc) = cpatch%hcapveg(nc) * (cpatch%veg_temp(nc)-t3ple)
    
    call init_ed_cohort_vars_array(cpatch, nc, lsl)
    
    cpatch%new_recruit_flag(nc) = 1 ! should plantations be considered recruits?

    ! Sort the cohorts so that the new cohort is at the correct height bin

    call sort_cohorts_ar(cpatch)
    
    return
  end subroutine plant_patch_ar
  
end module disturbance_utils_ar
