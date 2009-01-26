!========================================================================
! NOTICE:  These subroutines have not been thoroughly tested. Please
!          report problems to David Medvigy, medvigy@post.harvard.edu.
!========================================================================

subroutine apply_forestry_ar(cpoly, isi, year, rhos)

  use ed_state_vars,only : polygontype,sitetype,allocate_sitetype,deallocate_sitetype      &
                          ,copy_sitetype_mask

  use disturb_coms, only: lutime, min_new_patch_area, plantation_year
  use disturbance_utils_ar, only: initialize_disturbed_patch_ar, plant_patch_ar
  use fuse_fiss_utils_ar,only : terminate_patches_ar

  use max_dims, only: n_pft, n_dbh

  implicit none

  type(polygontype),target :: cpoly
  type(sitetype),pointer :: csite
  type(sitetype),pointer :: tempsite
  integer :: isi,ipft,idbh
  integer :: np
  real, intent(in) :: rhos
  integer :: iyear,useyear
  integer, intent(in) :: year
  real :: primary_harvest_target
  real :: secondary_harvest_target
  real :: total_harvest_target
  type(lutime), pointer :: clutime
  real :: total_site_biomass

  real :: area_mature_primary
  real :: agb_mature_primary
  real :: area_mature_secondary
  real :: agb_mature_secondary
  real :: area_mature_plantation
  real :: agb_mature_plantation
  real :: total_harvested_area
  real :: lambda_mature_primary
  real :: lambda_mature_secondary
  real :: lambda_mature_plantation
  real :: harvest_deficit
  integer,allocatable :: mask(:)

  csite => cpoly%site(isi)

  ! The site is harvested up to a target biomass.  Patches are harvested 
  ! first from those above the harvest_age, with equal rates.  If the 
  ! biomass target is not met, the remaining biomass is harvested from the 
  ! patches below the minimum age, starting with the oldest.
  ! Harvest rates are taken from George Hurtt's GLU, Global landuse files. 
  ! Elements 12 and 16 are secondary harvesting and 14 and 18 are primary.

  ! First, find the right year in the clutimes vector

  useyear = cpoly%num_landuse_years(isi)

  if(year >= cpoly%clutimes(1,isi)%landuse_year)then
     ! Loop over years
     find_lu_year: do iyear = 1,cpoly%num_landuse_years(isi)
        
        if (year.eq.cpoly%clutimes(iyear,isi)%landuse_year) then
           useyear = iyear
           exit find_lu_year
        endif
        
     enddo find_lu_year
  else
     print*,"The simulation is in a year occuring prior"
     print*,"to the time period of valid land-use defined years."
     print*,"stopping the simulation in apply_forestry"
     stop
  endif

  clutime => cpoly%clutimes(useyear,isi)
  
  ! Set primary and secondary targets based on current rates and unapplied
  ! harvest from previous years (memory)
  primary_harvest_target = clutime%landuse(14) + clutime%landuse(18) +  &
       cpoly%primary_harvest_memory(isi)
  secondary_harvest_target = clutime%landuse(12) + clutime%landuse(16) +  &
       cpoly%secondary_harvest_memory(isi)
  total_harvest_target = primary_harvest_target + secondary_harvest_target

  ! Decide whether or not to create a harvest patch: 
  ! (a) must have site agb > 0 and (b) harvest must exceed some minimum 
  ! threshold.  Note: this is not necessarily the total harvest area.  
  
  ! WILL THIS SUM OVER ALL DIMENSIONS INTO A SCALAR VALUE?
  total_site_biomass = 0.
  do ipft=1,n_pft
     do idbh=1,n_dbh
        total_site_biomass = total_site_biomass + cpoly%agb(ipft, idbh, isi)
     end do
  end do

  if(total_site_biomass == 0.0 .or.  &
       total_harvest_target <= total_site_biomass * min_new_patch_area)then
     ! Update memory and return
     cpoly%primary_harvest_memory(isi) = primary_harvest_target
     cpoly%secondary_harvest_memory(isi) = secondary_harvest_target
     return
  endif

  ! The following routines extends the allocation of the current site
  ! by one patch, and preserves the original patches
  ! ------------------------------------------------
  allocate(tempsite)
  allocate(mask(csite%npatches))
  call allocate_sitetype(tempsite,csite%npatches)
  mask(:) = 1
  call copy_sitetype_mask(csite,tempsite,mask,csite%npatches,csite%npatches)
  call deallocate_sitetype(csite)
  call allocate_sitetype(csite,tempsite%npatches + 1)
  call copy_sitetype_mask(tempsite,csite,mask,tempsite%npatches,tempsite%npatches)
  call deallocate_sitetype(tempsite)
  deallocate(tempsite)
  deallocate(mask)
  ! -------------------------------------------------


  
  np = csite%npatches

  ! Initialize the new patch (np) in the last position
  csite%dist_type(np) = 2
  call initialize_disturbed_patch_ar(csite,cpoly%met(isi)%atm_tmp,cpoly%met(isi)%atm_shv,np,1)

  ! Compute current stocks of agb in mature forests.
  call inventory_mat_forests_ar(cpoly,isi,   &
       area_mature_primary, agb_mature_primary, area_mature_secondary,  &
       agb_mature_secondary, area_mature_plantation, agb_mature_plantation)

  ! Compute the mature-forest harvest rates
  call mat_forest_harv_rates_ar(agb_mature_primary,  &
     agb_mature_secondary, agb_mature_plantation, primary_harvest_target,   &
     secondary_harvest_target, lambda_mature_primary,   &
     lambda_mature_secondary, lambda_mature_plantation, harvest_deficit)

  ! Apply harvesting to the mature stands
  call harv_mat_patches_ar(cpoly,isi, np, lambda_mature_primary,   &
       lambda_mature_secondary, lambda_mature_plantation)
  
  ! Compute harvested area from mature patches.  This is also updated
  ! in harvest_immature_patches().
  total_harvested_area = lambda_mature_primary * area_mature_primary +  &
       lambda_mature_secondary * area_mature_secondary +  &
       lambda_mature_plantation * area_mature_plantation
  call harv_immat_patches_ar(cpoly, isi, np, harvest_deficit, total_harvested_area)

  ! Now we know the area of the new patch, and can normalize the averaged
  ! patch quantities. But this is done only when the new patch area is significant
  ! otherwise, just terminate it.

  csite%area(np) = total_harvested_area
  if (total_harvested_area > min_new_patch_area) then
     write(unit=*,fmt='(a,1x,es12.5,1x,a,1x,i5)') &
         ' ---> Making new patch (harvesting), with area=',csite%area(np),' for dist_type=',2
     call norm_harv_patch_ar(csite,np)

     ! Plant the patch if it is a plantation.
     if(cpoly%plantation(isi) == 1 .and. year > plantation_year)then
        call plant_patch_ar(csite,np, cpoly%plantation_stocking_pft(isi),   &
             cpoly%plantation_stocking_density(isi), 2.0, cpoly%lsl(isi))
        csite%plantation(np) = 1
     else
     
     ! Including some near bare ground state there, otherwise it will be with no cohorts...
     !   David, do you see any problem by doing that? I included this just because the fast
     !   time dynamics was crashing with the zero-cohort patch. 
        write (unit=*,fmt='(a)') ' ----> Including a nearly-bare ground state for a patch with no cohorts...'
        call init_bare_ground_patchtype(.false.,csite,cpoly%lsl(isi),cpoly%met(isi)%atm_tmp,np,np)
     end if

     call update_patch_derived_props_ar(csite, cpoly%lsl(isi), rhos, np)
  
     call new_patch_sfc_props_ar(csite,np, rhos) 
  end if
  

  ! DAVID COULD YOU COMPARE THIS WITH THE LEGACY CODE?
  ! PREVIOUSLY, TERMINATE_PATCHES WAS CALLED PRIOR TO PLANT PATCH
  ! WHICH COULD HAVE POSSIBLY CHANGED THE NUMBER OF PATCHES IN THE
  ! THE SITE< AND POSSIBLY REMOVING THE NEW PATCH (NP)
  ! SO I MOVED IT TO THE END, SO THAT THE INDEXES WOULD BE PRESERVED
  ! FOR PLANT_PATCH, AND UPDATE_PATCH_DERIVED_PROPS. ALSO THE CODE FOR
  ! THE TERMINATION WAS LOOPING SITES, BUT BEING CALLED AT THE SITE LEVEL
  ! AND IT WAS NOT LOOPING OVER PATCHES, IT ONLY EVALUATED THE OLDEST PATCH
  ! IS THAT RIGHT? MAYBE IT WAS. NOW IT IS LOOPING OVER ALL PATCHES IN THE SITE,
  ! BUT ONLY CONSIDERS THE CURRENT SITE'S PATCHES (CSITE)

  ! If any patches now have zero area, terminate them.
  call terminate_patches_ar(csite)
  

  ! DAVID, SHOULD WE CHANGE THESE VARIABLES EVEN WHEN THE AREA WAS TINY?
  ! FROM DISTURBANCE.F90 IT SEEMS THAT WE SHOULD, BUT NOT SURE...
  
  ! Clear out the primary harvest memory.
  cpoly%primary_harvest_memory(isi) = 0.0
  
  ! There still may be a deficit if we have harvested all of the patch agb.
  cpoly%secondary_harvest_memory(isi) = harvest_deficit
  
  return
end subroutine apply_forestry_ar

!=================================================================================

subroutine inventory_mat_forests_ar(cpoly,isi, area_mature_primary,   &
     agb_mature_primary, area_mature_secondary, agb_mature_secondary,  &
     area_mature_plantation, agb_mature_plantation)
  
  use ed_state_vars,only:polygontype,sitetype,patchtype
  use disturb_coms, only: plantation_rotation, mature_harvest_age
  use allometry, only: ed_biomass
  implicit none


  type(polygontype),target :: cpoly
  type(sitetype),pointer   :: csite
  type(patchtype),pointer  :: cpatch
  integer :: isi,ipa,ico

  real, intent(out) :: area_mature_primary
  real, intent(out) :: agb_mature_primary
  real, intent(out) :: area_mature_secondary
  real, intent(out) :: agb_mature_secondary
  real, intent(out) :: area_mature_plantation
  real, intent(out) :: agb_mature_plantation


  ! Initialize inventory
  area_mature_primary = 0.0
  area_mature_secondary = 0.0
  area_mature_plantation = 0.0
  agb_mature_primary = 0.0
  agb_mature_secondary = 0.0
  agb_mature_plantation = 0.0

  csite => cpoly%site(isi)
  
  do ipa=1,csite%npatches
     
     cpatch => csite%patch(ipa)
     
     ! Compute the patch agb
     csite%plant_ag_biomass(ipa) = 0.0
     
     do ico=1,cpatch%ncohorts
        csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa) + ed_biomass(cpatch%bdead(ico),   &
             cpatch%balive(ico), cpatch%bleaf(ico), cpatch%pft(ico), cpatch%hite(ico), &
             cpatch%bstorage(ico)) * cpatch%nplant(ico)
     enddo

     if(csite%plant_ag_biomass(ipa) < 0.01)cycle

     ! Increment appropriate counter
     if(csite%plantation(ipa) == 1 .and. csite%age(ipa) > plantation_rotation)then
        
        ! Mature plantation
        area_mature_plantation = area_mature_plantation + csite%area(ipa)
        agb_mature_plantation = agb_mature_plantation +   &
             csite%plant_ag_biomass(ipa) * csite%area(ipa)

     elseif(csite%dist_type(ipa) == 2 .and. csite%plantation(ipa) /= 1 .and.   &
          csite%age(ipa) > mature_harvest_age)then

        ! Mature secondary
        area_mature_secondary = area_mature_secondary + csite%area(ipa)
        agb_mature_secondary = agb_mature_secondary +   &
             csite%plant_ag_biomass(ipa) * csite%area(ipa)

     elseif(csite%dist_type(ipa) == 3 .and. csite%age(ipa) > mature_harvest_age)then

        ! Mature primary
        area_mature_primary = area_mature_primary + csite%area(ipa)
        agb_mature_primary = agb_mature_primary +   &
             csite%plant_ag_biomass(ipa) * csite%area(ipa)

     endif

  enddo

  return
end subroutine inventory_mat_forests_ar

!======================================================================

subroutine mat_forest_harv_rates_ar(agb_mature_primary,  &
     agb_mature_secondary, agb_mature_plantation, primary_harvest_target,   &
     secondary_harvest_target, lambda_mature_primary,   &
     lambda_mature_secondary, lambda_mature_plantation, harvest_deficit)

  implicit none

  real, intent(in) :: agb_mature_primary
  real, intent(in) :: agb_mature_secondary
  real, intent(in) :: agb_mature_plantation
  real, intent(in) :: primary_harvest_target
  real, intent(inout) :: secondary_harvest_target
  real, intent(out) :: lambda_mature_primary
  real, intent(out) :: lambda_mature_plantation
  real, intent(out) :: lambda_mature_secondary
  real, intent(out) :: harvest_deficit

  ! Compute harvesting rate in mature primary forest.  If there is 
  ! not enough biomass to harvest, harvest what is possible from primary
  ! and attempt to harvest the remainder of the target from secondary.
  if(agb_mature_primary > primary_harvest_target)then
     lambda_mature_primary = primary_harvest_target / agb_mature_primary
  else
     lambda_mature_primary = 1.0
     harvest_deficit = primary_harvest_target - agb_mature_primary
     secondary_harvest_target = secondary_harvest_target + harvest_deficit
  endif

  ! Compute harvesting rate in mature plantations and mature secondary
  ! forests.  First try to remove all biomass from plantations.  If this
  ! is not possible, remove what you could and try to remove the remainder
  ! from mature secondary forests.  If this is again not possible, store
  ! what is leaf in harvest_deficit.
  if(agb_mature_plantation > secondary_harvest_target)then
     lambda_mature_plantation = secondary_harvest_target /   &
          agb_mature_plantation
     lambda_mature_secondary = 0.0
     harvest_deficit = 0.0
  else
     lambda_mature_plantation = 1.0
     harvest_deficit = secondary_harvest_target - agb_mature_plantation
     if(agb_mature_secondary > harvest_deficit)then
        lambda_mature_secondary = harvest_deficit / agb_mature_secondary
        harvest_deficit = 0.0
     else
        lambda_mature_secondary = 1.0
        harvest_deficit = harvest_deficit - agb_mature_secondary
     endif
  endif

  return
end subroutine mat_forest_harv_rates_ar

!======================================================================

subroutine harv_mat_patches_ar(cpoly,isi,np, lambda_mature_primary,   &
     lambda_mature_secondary, lambda_mature_plantation)
  
  use ed_state_vars, only: polygontype,sitetype,patchtype
  use disturb_coms, only: mature_harvest_age, plantation_rotation
  use disturbance_utils_ar, only: accum_dist_litt_ar, &
       increment_patch_vars_ar

  implicit none

  type(polygontype),target :: cpoly
  type(sitetype),pointer   :: csite
  type(patchtype),pointer  :: cpatch

  integer :: isi,ipa,np
  real :: dA
  real, intent(in) :: lambda_mature_plantation
  real, intent(in) :: lambda_mature_secondary
  real, intent(in) :: lambda_mature_primary

  ! Loop over patches
  csite => cpoly%site(isi)

  do ipa=1,csite%npatches

     cpatch => csite%patch(ipa)
     
     if(csite%plantation(ipa) == 1 .and. csite%age(ipa) > plantation_rotation)then
        ! Harvest mature plantations
        dA = csite%area(ipa) * lambda_mature_plantation
     elseif(csite%dist_type(ipa) == 2 .and. csite%plantation(ipa) /= 1 .and.   &
          csite%age(ipa) > mature_harvest_age)then
        ! Harvest mature secondary
        dA = csite%area(ipa) * lambda_mature_secondary
     elseif(csite%dist_type(ipa) == 3 .and. csite%age(ipa) > mature_harvest_age)then
        ! Harvest mature primary
        dA = csite%area(ipa) * lambda_mature_primary
     else
        dA = 0.0  ! Immature patches not harvested here.
     endif
     
     if(dA > 0.0 .and. csite%plant_ag_biomass(ipa) >= 0.01)then
        csite%area(ipa) = csite%area(ipa) - dA
        call increment_patch_vars_ar(csite,np,ipa, dA)
        call accum_dist_litt_ar(csite,np,ipa, 1, dA, &
             cpoly%loss_fraction(1,isi),cpoly%nat_dist_type(isi))
     endif

  enddo


  return
end subroutine harv_mat_patches_ar

!====================================================================

subroutine harv_immat_patches_ar(cpoly,isi, np, harvest_deficit,   &
     total_harvest_area)
  
  use ed_state_vars, only: polygontype,sitetype,patchtype
  use disturb_coms, only: plantation_rotation, mature_harvest_age
  use disturbance_utils_ar, only: accum_dist_litt_ar, &
       increment_patch_vars_ar

  implicit none

  type(polygontype),target :: cpoly
  type(sitetype),pointer   :: csite
  integer :: isi,ipa,np
  real, intent(inout) :: harvest_deficit

  real :: lambda
  real :: dA
  real, intent(inout) :: total_harvest_area

  ! Loop over patches
 
  csite => cpoly%site(isi)

  !  David, csite%patch(npatches) is the patch that has just been created due to harvesting. 
  !  It doesn't have area assigned yet, which means that it has 1E+34. In any case, I don't 
  !  think it falls into the "immature secondary" or "immature plantation" categories since 
  !  it is already harvested... so I switched the loop to stop at csite%npatches-1 rather 
  !  than csite%npatches. 

  do ipa=1,csite%npatches-1

     if(csite%plant_ag_biomass(ipa) < 0.01)cycle

     ! First harvest the immature secondary
     if(harvest_deficit > 0.0 .and.        &  ! There is still a deficit
          csite%dist_type(ipa) == 2 .and.  &  ! Secondary forest
          ( (csite%plantation(ipa) == 1 .and. csite%age(ipa) < plantation_rotation) .or.  &
          ! either immature plantation or immature secondary
          (csite%plantation(ipa) /= 1 .and. csite%age(ipa) < mature_harvest_age) ) ) then

        if( (csite%area(ipa) * csite%plant_ag_biomass(ipa)) > harvest_deficit)then

           ! Patch is not totally harvested
           lambda = harvest_deficit / (csite%area(ipa) * csite%plant_ag_biomass(ipa))
           dA = csite%area(ipa) * lambda
           harvest_deficit = 0.0

        else

           ! Patch is totally harvested
           dA = csite%area(ipa)
           harvest_deficit = harvest_deficit - csite%area(ipa) * csite%plant_ag_biomass(ipa)

        endif

        total_harvest_area = total_harvest_area + dA
        csite%area(ipa) = csite%area(ipa) - dA

        call increment_patch_vars_ar(csite,np,ipa, dA)
        call accum_dist_litt_ar(csite,np,ipa, 2, dA, &
             cpoly%loss_fraction(2,isi),cpoly%nat_dist_type(isi))
 
     endif


  enddo

  ! Return if we have reached our harvest target.
  if(harvest_deficit <= 0.0)return

  ! If we did not reach our target, loop again through patches, this time 
  ! harvesting from immature primary.
  do ipa=1,csite%npatches
  
     if(csite%plant_ag_biomass(ipa) < 0.01)cycle

     ! If necessary, harvest the immature primary
     if(harvest_deficit > 0.0 .and.  &  ! There is still a deficit
          csite%dist_type(ipa) == 3 .and.    &  ! and this is primary forest
          csite%age(ipa) < mature_harvest_age)then     ! and it is immature

        if( (csite%area(ipa) * csite%plant_ag_biomass(ipa)) > harvest_deficit)then
           
           ! Patch is not totally harvested
           lambda = harvest_deficit / (csite%area(ipa) * csite%plant_ag_biomass(ipa))
           dA = csite%area(ipa) * lambda
           harvest_deficit = 0.0

        else

           ! Patch is totally harvested
           dA = csite%area(ipa)
           harvest_deficit = harvest_deficit - csite%area(ipa) * csite%plant_ag_biomass(ipa)
           
        endif

        total_harvest_area = total_harvest_area + dA
        csite%area(ipa) = csite%area(ipa) - dA

        call increment_patch_vars_ar(csite,np,ipa, dA)
        call accum_dist_litt_ar(csite,np,ipa, 2, dA, &
             cpoly%loss_fraction(2,isi),cpoly%nat_dist_type(isi))

     endif

  enddo

  return
end subroutine harv_immat_patches_ar

!==================================================================

subroutine norm_harv_patch_ar(csite,np)

  use ed_state_vars, only: sitetype,patchtype
  use disturb_coms, only : min_new_patch_area

  use max_dims, only: n_pft
  use grid_coms, only: nzg, nzs

  implicit none
  type(sitetype),target   :: csite
  integer :: np
  real :: area_fac
  integer :: k

  !    This is to prevent FPE errors when the harvested area is insignificant or zero. 
  ! Also, if the area is small, the patch will be terminated soon, so I will skip the 
  ! normalization.
  
  if (csite%area(np) < min_new_patch_area) then
     return
  else
     area_fac = 1.0 / csite%area(np)
  end if

  csite%fast_soil_C(np) = csite%fast_soil_C(np) * area_fac

  csite%slow_soil_C(np) = csite%slow_soil_C(np) * area_fac

  csite%structural_soil_C(np) = csite%structural_soil_C(np) * area_fac

  csite%structural_soil_L(np) = csite%structural_soil_L(np) * area_fac

  csite%mineralized_soil_N(np) = csite%mineralized_soil_N(np) * area_fac

  csite%fast_soil_N(np) = csite%fast_soil_N(np) * area_fac

  csite%sum_dgd(np) = csite%sum_dgd(np) * area_fac

  csite%sum_chd(np) = csite%sum_chd(np) * area_fac

  csite%can_temp(np) = csite%can_temp(np) * area_fac

  csite%can_shv(np) = csite%can_shv(np) * area_fac

  csite%can_depth(np) = csite%can_depth(np) * area_fac

  do k = 1, nzs
     csite%sfcwater_mass(k,np) = csite%sfcwater_mass(k,np) * area_fac
     csite%sfcwater_energy(k,np) = csite%sfcwater_energy(k,np) * area_fac
     csite%sfcwater_depth(k,np) = csite%sfcwater_depth(k,np) * area_fac
  enddo

  do k = 1, nzg
     csite%soil_energy(k,np) = csite%soil_energy(k,np) * area_fac
     csite%soil_water(k,np) = csite%soil_water(k,np) * dble(area_fac)
  enddo

  csite%rough(np) = csite%rough(np) * area_fac

  csite%mean_rh(np) = csite%mean_rh(np) * area_fac

  csite%dmean_A_decomp(np) = csite%dmean_A_decomp(np) * area_fac

  csite%dmean_Af_decomp(np) = csite%dmean_Af_decomp(np) * area_fac

  csite%repro(1:n_pft,np) = csite%repro(1:n_pft,np) * area_fac

  csite%fsc_in(np) = csite%fsc_in(np) * area_fac

  csite%ssc_in(np) = csite%ssc_in(np) * area_fac

  csite%ssl_in(np) = csite%ssl_in(np) * area_fac

  csite%fsn_in(np) = csite%fsn_in(np) * area_fac

  csite%total_plant_nitrogen_uptake(np) = csite%total_plant_nitrogen_uptake(np) * area_fac

  return
end subroutine norm_harv_patch_ar

