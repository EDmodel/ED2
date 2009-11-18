!==========================================================================================!
!==========================================================================================!
!    This module contains subroutines and functions that will apply disturbances to        !
! patches.  This is usually done once a year, and the main disturbance driver will decide  !
! which kind of disturbance should be applied.                                             !
!------------------------------------------------------------------------------------------!
module disturbance_utils

   use ed_state_vars   , only : allocate_patchtype    & ! subroutine
                              , copy_patchtype        & ! subroutine
                              , deallocate_patchtype  & ! subroutine
                              , allocate_sitetype     & ! subroutine
                              , deallocate_sitetype   & ! subroutine
                              , copy_sitetype_mask    ! ! subroutine
   use fuse_fiss_utils , only : fuse_cohorts          & ! subroutine
                              , terminate_cohorts     & ! subroutine
                              , split_cohorts         ! ! subroutine
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !     This is the main disturbance driver.  It will be called every New Year day, and   !
   ! it will decide whether a new patch should be created.  Disturbances can be natural    !
   ! and anthropogenic, and the new patch will be always assigned an identification flag   !
   ! that will tell how that patch was created.  Three categories are currently possible:  !
   ! 1 - agriculture: conversion to agriculture by land clearing;                          !
   ! 2 - secondary forest: logging, land abandonment, and harvest create this patch;       !
   ! 3 - primary forest: natural disturbances (treefall or fire).                          !
   !---------------------------------------------------------------------------------------!
   subroutine apply_disturbances(cgrid)

      use ed_state_vars, only : edtype                  & ! structure
                              , polygontype             & ! structure
                              , sitetype                & ! structure
                              , patchtype               ! ! structure
      use ed_misc_coms , only : current_time            ! ! intent(in)
      use disturb_coms , only : treefall_age_threshold  & ! intent(in)
                              , min_new_patch_area      ! ! intent(in)
      use ed_max_dims  , only : n_dist_types            & ! intent(in)
                              , n_pft                   & ! intent(in)
                              , n_dbh                   ! ! intent(in)
      use mem_sites    , only : maxcohort               ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)                   , target      :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype)              , pointer     :: cpoly
      type(sitetype)                 , pointer     :: csite
      type(sitetype)                 , pointer     :: tsite
      type(patchtype)                , pointer     :: qpatch
      real   , dimension(n_pft,n_dbh)              :: initial_agb
      real   , dimension(n_pft,n_dbh)              :: initial_basal_area
      logical, dimension(:)          , allocatable :: disturb_mask
      integer                                      :: ipy
      integer                                      :: isi
      integer                                      :: ipa
      integer                                      :: onsp
      integer                                      :: q
      logical                                      :: ploughed
      logical                                      :: abandoned
      logical                                      :: natural
      real                                         :: area
      real                                         :: area_fac
      real                                         :: dA
      real                                         :: elim_nplant
      real                                         :: elim_lai
      !------------------------------------------------------------------------------------!

      !----- Allocating the temporary site that will host the original patches. -----------!
      allocate(tsite)

      polyloop: do ipy = 1,cgrid%npolygons
         
         cpoly => cgrid%polygon(ipy)
         siteloop: do isi = 1,cpoly%nsites

            csite => cpoly%site(isi)

            !----- Store AGB, basal area profiles in memory. ------------------------------!
            call update_site_derived_props(cpoly, 1,isi)
            initial_agb(1:n_pft,1:n_dbh) = cpoly%agb(1:n_pft,1:n_dbh,isi)
            initial_basal_area(1:n_pft,1:n_dbh) = cpoly%basal_area(1:n_pft,1:n_dbh,isi)
          
            !------------------------------------------------------------------------------!
            !      First take care of harvesting, i.e., secondary -> secondary and         !
            ! primary -> secondary.                                                        !
            !------------------------------------------------------------------------------!
            call apply_forestry(cpoly,isi, current_time%year)
          
            !----- Update the cut output variables. ---------------------------------------!
            call update_site_derived_props(cpoly, 1,isi)
            cpoly%agb_cut(1:n_pft,1:n_dbh,isi) = cpoly%agb_cut(1:n_pft, 1:n_dbh,isi)       &
                                               + initial_agb(1:n_pft, 1:n_dbh)             &
                                               - cpoly%agb(1:n_pft, 1:n_dbh,isi)

            cpoly%basal_area_cut(1:n_pft,1:n_dbh,isi) =                                    &
                                                 cpoly%basal_area_cut(1:n_pft,1:n_dbh,isi) &
                                               + initial_basal_area(1:n_pft,1:n_dbh)       &
                                               - cpoly%basal_area(1:n_pft,1:n_dbh,isi)

            !----- Save the Original Number (of) Site Patches, onsp... --------------------!
            onsp = csite%npatches


            !------------------------------------------------------------------------------!
            !     Create a temporary site with vectors containing all current patches as   !
            ! well as n_dist_types patches.  Create the newly disturbed patches in here,   !
            ! and depending on how many are created, repopulate the existing site's patch  !
            ! vectors.                                                                     !
            !------------------------------------------------------------------------------!
            call allocate_sitetype(tsite,onsp)

            allocate(disturb_mask(onsp + n_dist_types))
            disturb_mask         = .false.
            disturb_mask(1:onsp) = .true.

            !------------------------------------------------------------------------------!
            !     Transfer the origial patch values into the front end of the temp's       !
            ! space.                                                                       !
            !------------------------------------------------------------------------------!
            call copy_sitetype_mask(csite,tsite,disturb_mask(1:onsp),count(disturb_mask)   &
                                   ,count(disturb_mask))

            !----- Reallocate and transfer them back. -------------------------------------!
            call deallocate_sitetype(csite)
            call allocate_sitetype(csite,onsp + n_dist_types)
            call copy_sitetype_mask(tsite,csite,disturb_mask(1:onsp),count(disturb_mask)   &
                                   ,count(disturb_mask))
            call deallocate_sitetype(tsite)


            !------------------------------------------------------------------------------!
            !      Initialize all the potential as well as implemented disturbance         !
            ! patches.  n_dist_types new patches will be created, each one containing a    !
            ! different patch type.  In case no conversion to that kind of patch has       !
            ! happened, or if the newly created patch is tiny, it will be removed soon.    !                                                                   !
            !------------------------------------------------------------------------------!
            do q = onsp+1, onsp+n_dist_types
                call initialize_disturbed_patch(csite,cpoly%met(isi)%atm_tmp,q,1           &
                                               ,cpoly%lsl(isi))
            end do

            !----- Loop over q, the *destination* landuse type. ---------------------------!
            do q = 1, n_dist_types
               !----- Set up area to zero, in case no conversion happens. -----------------!
               area = 0.0

               do ipa=1,onsp
                  !------------------------------------------------------------------------!
                  !    Now we add the area associated with each kind of possible disturb-  !
                  ! ance that can happen.                                                  !
                  !------------------------------------------------------------------------!
                  ploughed  = q == 1 .and. csite%dist_type(ipa) /= 1  ! Conv. to agriculture
                  abandoned = q == 2 .and. csite%dist_type(ipa) == 1  ! Abandonment
                  !----- Natural disturbance, either trees are old or there is a fire. ----!
                  natural   = q == 3 .and. csite%dist_type(ipa) /= 1 .and.                 &
                              ( csite%age(ipa) > treefall_age_threshold .or.               &
                                cpoly%nat_dist_type(isi) == 1)
                  
                  if (ploughed .or. abandoned .or. natural) then
                     area = area                                                           &
                          + csite%area(ipa) * (1.0 - exp(                                  &
                          - (cpoly%disturbance_rates(q,csite%dist_type(ipa),isi)           &
                          + cpoly%disturbance_memory(q,csite%dist_type(ipa),isi))))
                  end if
               end do

               if (area > min_new_patch_area) then
                  write(unit=*,fmt='(a,1x,es12.5,1x,a,1x,i5)')                             &
                      ' ---> Making new patch, with area=',area,' for dist_type=',q

                  !------------------------------------------------------------------------!
                  !     Set the flag that this patch should be kept as a newly created     !
                  ! transition patch.                                                      !
                  !------------------------------------------------------------------------!
                  disturb_mask(onsp+q)     = .true.
  
                  csite%dist_type(onsp+q)  = q
                  csite%plantation(onsp+q) = 0
                  csite%area(onsp+q)       = area
                  
                  !----- Initialize to zero the new trasitioned patches. ------------------!
                  call initialize_disturbed_patch(csite,cpoly%met(isi)%atm_tmp,onsp+q,1    &
                                                 ,cpoly%lsl(isi))

                  !----- Now go through patches, adding its contribution to the new patch. !
                  do ipa=1,onsp
                     ploughed  = q == 1 .and. csite%dist_type(ipa) /= 1
                     abandoned = q == 2 .and. csite%dist_type(ipa) == 1
                     natural   = q == 3 .and. csite%dist_type(ipa) /= 1 .and.              &
                                 ( csite%age(ipa) > treefall_age_threshold .or.            &
                                   cpoly%nat_dist_type(isi) == 1)
                                     
                     if (ploughed .or. abandoned .or. natural) then
                        dA = csite%area(ipa) * (1.0 - exp(                                 &
                           - (cpoly%disturbance_rates(q,csite%dist_type(ipa),isi)          &
                           + cpoly%disturbance_memory(q,csite%dist_type(ipa),isi))))
                        area_fac = dA / csite%area(onsp+q)

                        call increment_patch_vars(csite,q+onsp,ipa,area_fac)
                        call insert_survivors(csite,q+onsp,ipa,q,area_fac                  &
                                             ,cpoly%nat_dist_type(isi))
                        call accum_dist_litt(csite,q+onsp,ipa,q,area_fac, &
                             cpoly%loss_fraction(q,isi),cpoly%nat_dist_type(isi))

                        !----- Update patch area. -----------------------------------------!
                        csite%area(ipa) = csite%area(ipa) - dA
                     end if
                  end do

                  !------------------------------------------------------------------------!
                  !      Update temperature and density.  This must be done before plant-  !
                  ! ing, since the leaf temperature is initially assigned as the canopy    !
                  ! air temperature.                                                       !
                  !------------------------------------------------------------------------!
                  call update_patch_thermo_props(csite,q+onsp,q+onsp)

                  !----- If the new patch is agriculture, plant it with grasses. ----------!
                  if (q == 1) then 
                     call plant_patch(csite,q+onsp,cpoly%agri_stocking_pft(isi)            &
                                     ,cpoly%agri_stocking_density(isi)                     &
                                     ,cpoly%green_leaf_factor(:,isi), 1.0, cpoly%lsl(isi))
                  end if

                  qpatch => csite%patch(q+onsp)

                  !----- Fuse then terminate cohorts. -------------------------------------!
                  if (csite%patch(q+onsp)%ncohorts > 0 .and. maxcohort >= 0) then
                     call fuse_cohorts(csite,q+onsp,cpoly%green_leaf_factor(:,isi)         &
                                      ,cpoly%lsl(isi))
                     call terminate_cohorts(csite,q+onsp,elim_nplant,elim_lai)
                     call split_cohorts(qpatch,cpoly%green_leaf_factor(:,isi)              &
                                       ,cpoly%lsl(isi))
                  end if
              
                  !----- Store AGB, basal area profiles in memory. ------------------------!
                  initial_agb(1:n_pft,1:n_dbh) = cpoly%agb(1:n_pft, 1:n_dbh,isi)
                  initial_basal_area(1:n_pft,1:n_dbh) =                                    &
                                                      cpoly%basal_area(1:n_pft,1:n_dbh,isi)

                  !------------------------------------------------------------------------!
                  !     Update the derived properties including veg_height, patch hcapveg, !
                  ! patch-level LAI, WAI, WPA.                                             !
                  !------------------------------------------------------------------------!
                  call update_patch_derived_props(csite,cpoly%lsl(isi),cpoly%met(isi)%prss &
                                                 ,q+onsp)
                  !----- Update soil temperature, liquid fraction, etc. -------------------!
                  call new_patch_sfc_props(csite, q+onsp)
                  !----- Update budget properties. ----------------------------------------!
                  call update_budget(csite,cpoly%lsl(isi),q+onsp,q+onsp)

                  !----- Update AGB, basal area. ------------------------------------------!
                  call update_site_derived_props(cpoly,1,isi)

                  !----- Update either cut or mortality. ----------------------------------!
                  if (q /= 3) then
                     cpoly%agb_cut(1:n_pft,1:n_dbh,isi) =                                  &
                            cpoly%agb_cut(1:n_pft, 1:n_dbh,isi)                            &
                          + initial_agb(1:n_pft, 1:n_dbh)                                  &
                          - cpoly%agb(1:n_pft, 1:n_dbh,isi)
                     cpoly%basal_area_cut(1:n_pft, 1:n_dbh,isi) =                          &
                            cpoly%basal_area_cut(1:n_pft, 1:n_dbh,isi)                     &
                          + initial_basal_area(1:n_pft, 1:n_dbh)                           &
                          - cpoly%basal_area(1:n_pft, 1:n_dbh,isi)
                  else
                     cpoly%agb_mort(1:n_pft,1:n_dbh,isi) =                                 &
                            cpoly%agb_mort(1:n_pft,1:n_dbh,isi)                            &
                          + initial_agb(1:n_pft,1:n_dbh)                                   &
                          - cpoly%agb(1:n_pft,1:n_dbh,isi)
                     cpoly%basal_area_mort(1:n_pft, 1:n_dbh,isi) =                         &
                            cpoly%basal_area_mort(1:n_pft, 1:n_dbh,isi)                    &
                          + initial_basal_area(1:n_pft, 1:n_dbh)                           &
                          - cpoly%basal_area(1:n_pft, 1:n_dbh,isi)
                  endif
                  
                  !----- Clear the disturbance memory for this disturbance type. ----------!
                  cpoly%disturbance_memory(q,1:n_dist_types,isi) = 0.0

               elseif(area > 0.0)then
                  !------------------------------------------------------------------------!
                  !     The patch creation has been skipped because the area was too       !
                  ! small.  Put the current disturbance rates in memory to be added at the !
                  ! next timestep.                                                         !
                  !------------------------------------------------------------------------!
                  cpoly%disturbance_memory(q,1:n_dist_types,isi) =                         &
                         cpoly%disturbance_memory(q,1:n_dist_types,isi)                    &
                       + cpoly%disturbance_rates(q,1:n_dist_types,isi)
               end if
            end do

            !------------------------------------------------------------------------------!
            !      Reallocate the current site to fit the original patches and whatever    !
            ! was generated in disturbance (ie, make it non sparse).  Populate the         !
            ! original site with both the modified original patches, and the newly created !
            ! patches.  The index of all of these are disturb_mask.  This mask should be   !
            ! ones for all original patches, and sparse from there after.                  !
            !------------------------------------------------------------------------------!
            call allocate_sitetype(tsite,count(disturb_mask))
            call copy_sitetype_mask(csite,tsite,disturb_mask,size(disturb_mask)            &
                                   ,count(disturb_mask))
            call deallocate_sitetype(csite)
            call allocate_sitetype(csite,count(disturb_mask))
          
            disturb_mask = .false.
            disturb_mask(1:csite%npatches) = .true.
            call copy_sitetype_mask(tsite,csite,disturb_mask(1:csite%npatches)             &
                                   ,count(disturb_mask),count(disturb_mask))
          
            call deallocate_sitetype(tsite)
            deallocate(disturb_mask)
            !------------------------------------------------------------------------------!

          
         end do siteloop
      end do polyloop

      !----- Free memory before leaving... ------------------------------------------------!
      deallocate(tsite)

      return
   end subroutine apply_disturbances
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
  subroutine site_disturbance_rates(month, year, cgrid)

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

    integer :: iyear,useyear
    type(lutime), pointer :: clutime
    real :: fire_dist_rate

    ! Loop over sites and polygons.
    
    do ipy = 1,cgrid%npolygons
       
       cpoly => cgrid%polygon(ipy)
       
       do isi = 1,cpoly%nsites
          
          csite => cpoly%site(isi)

          !  Calculate fire disturbance rates only if fire is on.
          select case (include_fire)
          case (0) 
             fire_dist_rate = 0.0
          case (1,2)
             fire_dist_rate = sum(cpoly%lambda_fire(1:12,isi)) / 12.0
          end select

          cpoly%fire_disturbance_rate(isi) = fire_dist_rate
    
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
  end subroutine site_disturbance_rates

  !=====================================================================

  subroutine initialize_disturbed_patch(csite,atm_tmp,np,dp,lsl)
    
    use ed_state_vars, only: sitetype,patchtype
    use consts_coms, only : t00
    use grid_coms, only: nzs, nzg

    implicit none
    type(sitetype),target    :: csite
    real   , intent(in)      :: atm_tmp
    integer, intent(in)      :: np,dp,lsl
    integer                  :: k


    ! SHOULD WE CONSIDER USING A DONOR PATCH
    ! TO INITIALIZE MORE OF THE TRANSIENT STATE VARIABLES
    ! PRIOR TO THIS MESSAGE, NTEXT_SOIL WAS NOT BEING SET PROPERLY


    csite%patch(np)%ncohorts = 0

    ! For now, choose heat/vapor capacities for stability
    do k=1,nzs
       csite%sfcwater_tempk(k,np) = atm_tmp   ! Set canopy temp to 0 C
       csite%sfcwater_fracliq(k,np) = 1.0     ! Set to 100% liquid
    end do


    ! Included the initialization for all fast variables.
    call init_ed_patch_vars(csite,np,np,lsl)

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

    csite%can_theta(np) = 0.0

    csite%can_prss(np) = 0.0

    csite%can_shv(np) = 0.0

    csite%can_co2(np) = 0.0

    csite%soil_energy(1:nzg,np) = 0.0

    csite%soil_water(1:nzg,np) = 0.0

    csite%sfcwater_mass(1:nzs,np) =  0.0

    csite%sfcwater_energy(1:nzs,np) = 0.0

    csite%sfcwater_depth(1:nzs,np) = 0.0
    
    csite%hcapveg(np) = 0.0

    !--------------------------------------------------------------------------------------!
    csite%rough(np) = 0.0

    ! this 
    call init_ed_patch_vars(csite,np,np,lsl)

    csite%fsc_in(np) = 0.0

    csite%ssc_in(np) = 0.0

    csite%ssl_in(np) = 0.0

    csite%fsn_in(np) = 0.0

    csite%total_plant_nitrogen_uptake(np) = 0.0

    return
  end subroutine initialize_disturbed_patch

  !======================================================================

  subroutine increment_patch_vars(csite,np, cp, area_fac)

    use ed_state_vars, only: sitetype,patchtype
    use ed_max_dims, only: n_pft
    use grid_coms, only: nzg

    implicit none
    type(sitetype),target    :: csite
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

    csite%can_theta(np) = csite%can_theta(np) + csite%can_theta(cp) * area_fac

    csite%can_prss(np) = csite%can_prss(np) + csite%can_prss(cp) * area_fac

    csite%can_shv(np) = csite%can_shv(np) + csite%can_shv(cp) * area_fac

    csite%can_co2(np) = csite%can_co2(np) + csite%can_co2(cp) * area_fac

    csite%can_depth(np) = csite%can_depth(np) + csite%can_depth(cp) * area_fac

    do k = 1, csite%nlev_sfcwater(cp)
       csite%sfcwater_mass(k,np) = csite%sfcwater_mass(k,np) + csite%sfcwater_mass(k,cp) *   &
            area_fac

       csite%sfcwater_energy(k,np) = csite%sfcwater_energy(k,np) +   &
            csite%sfcwater_energy(k,cp) * csite%sfcwater_mass(k,cp) * area_fac
       csite%sfcwater_depth(k,np) = csite%sfcwater_depth(k,np) + csite%sfcwater_depth(k,cp) *   &
            area_fac
    enddo

    do k = 1, nzg
       csite%ntext_soil(k,np)  = csite%ntext_soil(k,np)
       csite%soil_energy(k,np) = csite%soil_energy(k,np) + csite%soil_energy(k,cp) * area_fac
       csite%soil_water(k,np)  = csite%soil_water(k,np)  + csite%soil_water(k,cp)  * area_fac
    enddo

    csite%rough(np) = csite%rough(np) + csite%rough(cp) * area_fac

    csite%mean_rh(np) = csite%mean_rh(np) + csite%mean_rh(cp) * area_fac

    csite%today_A_decomp(np) = csite%today_A_decomp(np) + csite%today_A_decomp(cp) * area_fac

    csite%today_Af_decomp(np) = csite%today_Af_decomp(np) + csite%today_Af_decomp(cp) * area_fac

    csite%repro(1:n_pft,np) = csite%repro(1:n_pft,np) + csite%repro(1:n_pft,cp) * area_fac

    csite%fsc_in(np) = csite%fsc_in(np) + csite%fsc_in(cp) * area_fac

    csite%ssc_in(np) = csite%ssc_in(np) + csite%ssc_in(cp) * area_fac

    csite%ssl_in(np) = csite%ssl_in(np) + csite%ssl_in(cp) * area_fac

    csite%fsn_in(np) = csite%fsn_in(np) + csite%fsn_in(cp) * area_fac

    csite%total_plant_nitrogen_uptake(np) = csite%total_plant_nitrogen_uptake(np) +   &
         csite%total_plant_nitrogen_uptake(cp) * area_fac

    return
  end subroutine increment_patch_vars

  !=======================================================================

  subroutine insert_survivors(csite, np, cp, q, area_fac,nat_dist_type)

    use ed_state_vars, only: sitetype,patchtype
    use ed_misc_coms , only: idoutput,imoutput
    
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

    logical, dimension(:), allocatable :: mask

    cpatch => csite%patch(cp)
    npatch => csite%patch(np)
    
    allocate(tpatch)

    allocate(mask(cpatch%ncohorts))
    mask(:) = .false.
    
    do ico = 1,cpatch%ncohorts
       
       cohort_area_fac = survivorship(q,nat_dist_type, csite, cp, ico) * area_fac
       n_survivors = cpatch%nplant(ico) * cohort_area_fac
       
       ! If something survived, make a new cohort
       mask(ico) = n_survivors > 0.0
       
    enddo

    if (npatch%ncohorts > 0) then
       
       ! If the new patch has received survivors from a donor
       ! already, then it should have cohorts. So the temp patch
       ! vector will be the sum of the new cohorts found here, plus
       ! those already applied previously in the loop calling this
       
       nco = npatch%ncohorts
       call allocate_patchtype(tpatch,count(mask) + npatch%ncohorts)
       call copy_patchtype(npatch,tpatch,1,npatch%ncohorts,1,npatch%ncohorts)
       call deallocate_patchtype(npatch)
       
    else
       
       nco = 0
       call allocate_patchtype(tpatch,count(mask))

    endif


    do ico = 1,cpatch%ncohorts
       
       cohort_area_fac = survivorship(q,nat_dist_type, csite, cp, ico) * area_fac
       n_survivors = cpatch%nplant(ico) * cohort_area_fac
       
       if(mask(ico)) then

          nco = nco + 1
          
          call copy_patchtype(cpatch,tpatch,ico,ico,nco,nco)

          ! Adjust area-based variables
          tpatch%lai(nco)                 = tpatch%lai(nco)                  * cohort_area_fac
          tpatch%wai(nco)                 = tpatch%wai(nco)                  * cohort_area_fac
          tpatch%wpa(nco)                 = tpatch%wpa(nco)                  * cohort_area_fac
          tpatch%nplant(nco)              = tpatch%nplant(nco)               * cohort_area_fac
          tpatch%mean_gpp(nco)            = tpatch%mean_gpp(nco)             * cohort_area_fac
          tpatch%mean_leaf_resp(nco)      = tpatch%mean_leaf_resp(nco)       * cohort_area_fac
          tpatch%mean_root_resp(nco)      = tpatch%mean_root_resp(nco)       * cohort_area_fac
          tpatch%today_gpp(nco)           = tpatch%today_gpp(nco)            * cohort_area_fac
          tpatch%today_gpp_pot(nco)       = tpatch%today_gpp_pot(nco)        * cohort_area_fac
          tpatch%today_gpp_max(nco)       = tpatch%today_gpp_max(nco)        * cohort_area_fac
          tpatch%today_leaf_resp(nco)     = tpatch%today_leaf_resp(nco)      * cohort_area_fac
          tpatch%today_root_resp(nco)     = tpatch%today_root_resp(nco)      * cohort_area_fac
          tpatch%Psi_open(nco)            = tpatch%Psi_open(nco)             * cohort_area_fac
          tpatch%gpp(nco)                 = tpatch%gpp(nco)                  * cohort_area_fac
          tpatch%leaf_respiration(nco)    = tpatch%leaf_respiration(nco)     * cohort_area_fac
          tpatch%root_respiration(nco)    = tpatch%root_respiration(nco)     * cohort_area_fac
          tpatch%monthly_dndt(nco)        = tpatch%monthly_dndt(nco)         * cohort_area_fac
          tpatch%veg_water(nco)           = tpatch%veg_water(nco)            * cohort_area_fac
          tpatch%hcapveg(nco)             = tpatch%hcapveg(nco)              * cohort_area_fac
          tpatch%veg_energy(nco)          = tpatch%veg_energy(nco)           * cohort_area_fac
          !----- Carbon flux monthly means are extensive, we must convert them. ------------!
          if (idoutput > 0 .or. imoutput > 0) then
             tpatch%dmean_par_v     (nco) = tpatch%dmean_par_v     (nco)     * cohort_area_fac
             tpatch%dmean_par_v_beam(nco) = tpatch%dmean_par_v_beam(nco)     * cohort_area_fac
             tpatch%dmean_par_v_diff(nco) = tpatch%dmean_par_v_diff(nco)     * cohort_area_fac
          end if
          if (imoutput > 0) then
             tpatch%mmean_par_v     (nco) = tpatch%mmean_par_v     (nco)     * cohort_area_fac
             tpatch%mmean_par_v_beam(nco) = tpatch%mmean_par_v_beam(nco)     * cohort_area_fac
             tpatch%mmean_par_v_diff(nco) = tpatch%mmean_par_v_diff(nco)     * cohort_area_fac
          end if
       end if
          
    enddo  ! end loop over cohorts

    ! Copy the temporary patch into the newpatch
    
    call allocate_patchtype(npatch,tpatch%ncohorts)
    call copy_patchtype(tpatch,npatch,1,tpatch%ncohorts,1,tpatch%ncohorts)
    call deallocate_patchtype(tpatch)
    deallocate(tpatch)
    
    deallocate(mask)


    return
  end subroutine insert_survivors

  !======================================================================
  subroutine accum_dist_litt(csite,np,cp, q, area_fac,loss_fraction,nat_dist_type)
    
    use ed_state_vars, only: sitetype,patchtype,polygontype
    use decomp_coms, only: f_labile
    use ed_max_dims, only: n_pft
    use pft_coms, only: c2n_storage, c2n_leaf, c2n_recruit, c2n_stem, l2n_stem
    use grid_coms, only: nzg

    implicit none

    type(sitetype),target    :: csite
    type(patchtype),pointer  :: cpatch
    type(patchtype),pointer  :: npatch
    integer :: np,cp,ico
    real,intent(in) :: loss_fraction

    integer, intent(in) :: q
    real, intent(in) :: area_fac
    integer,intent(in) :: nat_dist_type
    real :: fast_litter
    real :: struct_litter, struct_lignin
    real :: fast_litter_n
    real :: struct_cohort

    fast_litter = 0.0
    struct_litter = 0.0
    struct_lignin = 0.0
    fast_litter_n = 0.0

    cpatch => csite%patch(cp)
    npatch => csite%patch(np)

    do ico = 1,cpatch%ncohorts

       fast_litter = fast_litter + (1.0 - &
            survivorship(q,nat_dist_type,csite, cp, ico)) * &
            (f_labile(cpatch%pft(ico)) * cpatch%balive(ico) &
            + cpatch%bstorage(ico)) * cpatch%nplant(ico)

       fast_litter_n = fast_litter_n + (1.0 - &
            survivorship(q,nat_dist_type, csite, cp, ico)) * &
            (f_labile(cpatch%pft(ico)) * cpatch%balive(ico) &
            / c2n_leaf(cpatch%pft(ico)) + cpatch%bstorage(ico) /  &
            c2n_storage ) * cpatch%nplant(ico)

       struct_cohort = cpatch%nplant(ico) *   &
            (1.0 - survivorship(q,nat_dist_type, csite, cp, ico)) * ( (1.0 -   &
            loss_fraction ) * cpatch%bdead(ico) +   & ! DOUBLE CHECK THIS, IS THIS RIGHT??
            (1.0 - f_labile(cpatch%pft(ico))) * cpatch%balive(ico))

       struct_litter = struct_litter + struct_cohort
       struct_lignin = struct_lignin + struct_cohort * l2n_stem / c2n_stem(cpatch%pft(ico))
       
    enddo

    !  Load disturbance litter directly into carbon and N pools
    csite%fast_soil_C(np) = csite%fast_soil_C(np) + fast_litter * area_fac

    csite%structural_soil_C(np) = csite%structural_soil_C(np) + struct_litter * area_fac

    csite%structural_soil_L(np) = csite%structural_soil_L(np) +  &
         struct_lignin * area_fac

    csite%fast_soil_N(np) = csite%fast_soil_N(np) + fast_litter_n * area_fac

    return
  end subroutine accum_dist_litt

  !============================================================================

  real function survivorship(dest_type,poly_dest_type,csite, ipa, ico)

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
       survivorship = 0.0  !  agric
    else
       if(dest_type == 2)then
          survivorship = 0.0   ! secondary land
       else  !   natural land 
          if(dest_type == 3)then
             if(poly_dest_type == 1)then
                survivorship = 0.0 ! fire
             elseif(poly_dest_type == 0)then
                if(cpatch%hite(ico) < treefall_hite_threshold)then
                   survivorship =  treefall_s_ltht(cpatch%pft(ico)) 
                else 
                   survivorship = treefall_s_gtht(cpatch%pft(ico))
                endif
             endif
          endif
       endif
    endif
    
    return
  end function survivorship

  !=============================================================

  subroutine plant_patch(csite, np, pft, density, green_leaf_factor, height_factor, lsl)

    use ed_state_vars,only : sitetype,patchtype

    use pft_coms, only: q, qsw, sla, hgt_min, max_dbh
    use ed_misc_coms, only: dtlsm
    use fuse_fiss_utils, only : sort_cohorts
    use ed_therm_lib,only : calc_hcapveg
    use consts_coms, only: t3ple, pio4
    use allometry, only : h2dbh, dbh2bd, dbh2bl, dbh2h, area_indices, ed_biomass
    use ed_max_dims, only : n_pft

    implicit none

    type(patchtype),pointer  :: cpatch
    type(patchtype),pointer  :: tpatch
    type(sitetype),target   :: csite

    integer, intent(in) :: lsl
    integer :: np,nc
    integer, intent(in) :: pft
    real, intent(in) :: density
    real, intent(in) :: height_factor
    real, dimension(n_pft), intent(in) :: green_leaf_factor
    real :: salloc, salloci
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
    cpatch%hite(nc) = hgt_min(cpatch%pft(nc)) * min(1.0,height_factor)
    if(.false.) then !! it's a tree
       cpatch%dbh(nc) = h2dbh(cpatch%hite(nc),cpatch%pft(nc))
       cpatch%bdead(nc) = dbh2bd(cpatch%dbh(nc),cpatch%hite(nc),cpatch%pft(nc))
       cpatch%bleaf(nc) = dbh2bl(cpatch%dbh(nc),cpatch%pft(nc))
    else
       !! set actual bleaf
       cpatch%dbh(nc)   = h2dbh(cpatch%hite(nc),cpatch%pft(nc))
       cpatch%bleaf(nc) = dbh2bl(cpatch%dbh(nc),cpatch%pft(nc))

       !! reset allometry to make it grow
       cpatch%dbh(nc)   = max_dbh(pft)
       cpatch%bdead(nc) = dbh2bd(cpatch%dbh(nc),cpatch%hite(nc),cpatch%pft(nc))
       cpatch%hite(nc)  = dbh2h(pft,max_dbh(pft))
    end if

    cpatch%phenology_status(nc) = 0
    salloc = 1.0 + q(cpatch%pft(nc)) + qsw(cpatch%pft(nc)) * cpatch%hite(nc)
    salloci = 1.0 / salloc
    cpatch%balive(nc)   = cpatch%bleaf(nc) * salloc
    cpatch%broot(nc)    = cpatch%balive(nc) * q(cpatch%pft(nc)) * salloci
    cpatch%bsapwood(nc) = cpatch%balive(nc) * qsw(cpatch%pft(nc))                &
                        * cpatch%hite(nc) * salloci
    cpatch%sla(nc)=sla(cpatch%pft(nc))
    call area_indices(cpatch%nplant(nc),cpatch%bleaf(nc),cpatch%bdead(nc)        &
                     ,cpatch%balive(nc),cpatch%dbh(nc), cpatch%hite(nc)          &
                     ,cpatch%pft(nc),cpatch%sla(nc), cpatch%lai(nc)              &
                     ,cpatch%wpa(nc),cpatch%wai(nc))

    cpatch%bstorage(nc) = 1.0*(cpatch%balive(nc)) !! changed by MCD, was 0.0


    !----- Finding the new basal area and above-ground biomass. ---------------------------!
    cpatch%basarea(nc) = pio4 * cpatch%dbh(nc) * cpatch%dbh(nc)                
    cpatch%agb(nc)     = ed_biomass(cpatch%bdead(nc),cpatch%balive(nc),cpatch%bleaf(nc)    &
                                   ,cpatch%pft(nc),cpatch%hite(nc) ,cpatch%bstorage(nc))     


    call init_ed_cohort_vars(cpatch, nc, lsl)

    cpatch%veg_temp(nc)  = csite%can_temp(np)
    cpatch%veg_water(nc) = 0.0
    cpatch%veg_fliq(nc)  = 0.0

    !----- Because we assigned no water, the internal energy is simply hcapveg*T

    cpatch%hcapveg(nc) = calc_hcapveg(cpatch%bleaf(nc),cpatch%bdead(nc)   &
                                     ,cpatch%balive(nc),cpatch%nplant(nc) &
                                     ,cpatch%hite(nc),cpatch%pft(nc)      &
                                     ,cpatch%phenology_status(nc))
    
    cpatch%veg_energy(nc) = cpatch%hcapveg(nc) * cpatch%veg_temp(nc)
    
    
    cpatch%new_recruit_flag(nc) = 1 ! should plantations be considered recruits?

    ! Sort the cohorts so that the new cohort is at the correct height bin

    call sort_cohorts(cpatch)
    
    return
  end subroutine plant_patch
  
end module disturbance_utils
