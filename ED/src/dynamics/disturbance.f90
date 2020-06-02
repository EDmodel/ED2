!==========================================================================================!
!==========================================================================================!
!    This module contains subroutines and functions that will apply disturbances to        !
! patches.  This is usually done once a year, and the main disturbance driver will decide  !
! which kind of disturbance should be applied.                                             !
!------------------------------------------------------------------------------------------!
module disturbance

   use ed_state_vars   , only : allocate_patchtype    & ! subroutine
                              , copy_patchtype        & ! subroutine
                              , deallocate_patchtype  & ! subroutine
                              , allocate_sitetype     & ! subroutine
                              , deallocate_sitetype   & ! subroutine
                              , copy_sitetype_mask    & ! subroutine
                              , copy_sitetype         ! ! subroutine
   use fuse_fiss_utils , only : sort_cohorts          & ! subroutine
                              , old_fuse_cohorts      & ! subroutine
                              , new_fuse_cohorts      & ! subroutine
                              , terminate_cohorts     & ! subroutine
                              , split_cohorts         & ! subroutine
                              , fuse_2_patches        & ! subroutine
                              , fuse_2_cohorts        ! ! subroutine
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
      use ed_state_vars       , only : edtype                     & ! structure
                                     , polygontype                & ! structure
                                     , sitetype                   & ! structure
                                     , patchtype                  ! ! structure
      use met_driver_coms     , only : met_driv_state             ! ! structure
      use ed_misc_coms        , only : current_time               & ! intent(in)
                                     , ibigleaf                   & ! intent(in)
                                     , lianas_included            ! ! intent(in)
      use disturb_coms        , only : min_patch_area             & ! intent(in)
                                     , plantation_year            & ! intent(in)
                                     , treefall_hite_threshold    & ! intent(in)
                                     , min_oldgrowth              & ! intent(in)
                                     , sl_skid_rel_area           ! ! intent(in)
      use ed_max_dims         , only : n_dist_types               & ! intent(in)
                                     , n_pft                      & ! intent(in)
                                     , n_dbh                      ! ! intent(in)
      use mem_polygons        , only : maxcohort                  ! ! intent(in)
      use grid_coms           , only : nzg                        & ! intent(in)
                                     , nzs                        ! ! intent(in)
      use pft_coms            , only : include_pft                ! ! intent(in)
      use allometry           , only : area_indices               ! ! function
      use mortality           , only : disturbance_mortality      ! ! subroutine
      use consts_coms         , only : lnexp_max                  & ! intent(in)
                                     , tiny_num                   & ! intent(in)
                                     , huge_num                   ! ! intent(in)
      use forestry            , only : find_lambda_harvest        ! ! sub-routine
      use detailed_coms       , only : idetailed                  ! ! intent(in)
      use update_derived_utils, only : update_patch_thermo_props  & ! subroutine
                                     , update_patch_thermo_fmean  & ! subroutine
                                     , update_patch_derived_props & ! subroutine
                                     , update_site_derived_props  ! ! subroutine
      use fusion_fission_coms , only : ifusion                    ! ! intent(in)
      use ed_type_init        , only : new_patch_sfc_props        ! ! subroutine
      use stable_cohorts      , only : is_resolvable              ! ! subroutine
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)                    , target      :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype)               , pointer     :: cpoly
      type(met_driv_state)            , pointer     :: cmet
      type(sitetype)                  , pointer     :: csite
      type(sitetype)                  , pointer     :: tsite
      type(patchtype)                 , pointer     :: cpatch
      type(patchtype)                 , pointer     :: qpatch
      integer, dimension(n_pft)                     :: allind
      integer                                       :: ipy
      integer                                       :: isi
      integer                                       :: ipa
      integer                                       :: npa
      integer                                       :: ico
      integer                                       :: ipft
      integer                                       :: i
      integer                                       :: apa
      integer                                       :: zpa
      integer                                       :: mypfts
      integer                                       :: onsp
      integer                                       :: nnsp_ble
      integer                                       :: old_lu
      integer                                       :: new_lu
      integer                                       :: bfus_ncoh
      integer, dimension(:)           , allocatable :: pfts
      logical, dimension(:)           , allocatable :: disturb_mask
      real   , dimension(:)           , allocatable :: original_area
      real   , dimension(:)           , allocatable :: original_lu
      real   , dimension(:)           , allocatable :: lambda_harvest
      real   , dimension(:,:)         , allocatable :: lambda_arr
      real   , dimension(:,:)         , allocatable :: pot_area_loss
      real   , dimension(:,:)         , allocatable :: act_area_loss
      real   , dimension(n_dist_types)              :: lambda_now
      real   , dimension(n_dist_types)              :: pat_area_loss
      real   , dimension(n_dist_types)              :: pot_area_gain
      real   , dimension(n_dist_types)              :: act_area_gain
      real   , dimension(n_dist_types)              :: one_area_loss
      logical                                       :: disturbed
      logical                                       :: same_pft
      logical                                       :: is_primary
      real   , dimension(n_pft)                     :: mindbh_harvest
      real                                          :: pot_area_remain
      real                                          :: area_loss_tot
      real                                          :: lambda_sum
      real                                          :: area_fac
      real                                          :: orig_area
      real                                          :: dist_area
      real                                          :: elim_nplant
      real                                          :: elim_lai
      real                                          :: new_nplant
      real                                          :: lfactor
      real                                          :: cb_enthalpy_np
      real                                          :: can_exner_np
      real                                          :: cb_mass_np
      real                                          :: cb_molar_np
      real                                          :: cb_water_np
      real                                          :: cb_co2_np
      logical                                       :: is_managed
      logical                                       :: is_oldgrowth
      logical                                       :: is_plantation
      logical                                       :: print_detailed
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Create a temporary vector with the allowed PFTs.  This will be later used      !
      ! to define the new patches in case this is a big leaf simulation.                   !
      !------------------------------------------------------------------------------------!
      mypfts = count(include_pft)
      allocate(pfts(mypfts))
      allind = (/ (i, i=1, n_pft) /)
      pfts   = pack(allind,include_pft)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !       Find out whether to print detailed information on screen.                    !
      !------------------------------------------------------------------------------------!
      print_detailed = btest(idetailed,6)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      nnsp_ble stands for the new number of site patches for big leaf.  Each site   !
      ! may have additional mypfts, except if land use is 1 (cropland/pasture), 2          !
      ! (forest plantation), or 8 (cropland), in which case only one PFT is allowed.       !
      !------------------------------------------------------------------------------------!
      nnsp_ble = 3 + (n_dist_types - 3) * mypfts
      !------------------------------------------------------------------------------------!



      !----- Allocate the temporary site that will host the original patches. -------------!
      nullify (tsite)
      allocate(tsite)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Loop over polygons and sites.                                                  !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)


         !----- Reset site-level harvest variables. ---------------------------------------!
         cgrid%crop_harvest   (ipy) = 0.0
         cgrid%logging_harvest(ipy) = 0.0
         cgrid%combusted_fuel (ipy) = 0.0
         !---------------------------------------------------------------------------------!


         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)
            cmet  => cpoly%met(isi)

            !----- Save the Original Number (of) Site Patches, onsp... --------------------!
            onsp = csite%npatches
            !------------------------------------------------------------------------------!


            !----- Test whether the new explored patch will be plantation or logged. ------!
            is_plantation = cpoly%plantation(isi) == 1               .and.                 &
                            current_time%year     >  plantation_year
            !------------------------------------------------------------------------------!


            !----- Reset site-level harvest variables. ------------------------------------!
            cpoly%crop_harvest   (isi) = 0.0
            cpoly%logging_harvest(isi) = 0.0
            cpoly%combusted_fuel (isi) = 0.0
            !------------------------------------------------------------------------------!

            !----- Copy mindbh for harvesting to a local variable. ------------------------!
            mindbh_harvest  (:) = cpoly%mindbh_harvest  (:,isi)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Allocate the area losses and gains.  We define potential and actual      !
            ! so we can correct patches that would become too small or new patches that    !
            ! would be too small to be generated.                                          !
            !------------------------------------------------------------------------------!
            allocate (original_area  (onsp)             )
            allocate (original_lu    (onsp)             )
            allocate (lambda_harvest (onsp)             )
            allocate (lambda_arr     (onsp,n_dist_types))
            allocate (pot_area_loss  (onsp,n_dist_types))
            allocate (act_area_loss  (onsp,n_dist_types))
            original_area  (:  ) = 0.0
            original_lu    (:  ) = 0.0
            lambda_harvest (:  ) = 0.0
            pot_area_gain  (  :) = 0.0
            act_area_gain  (  :) = 0.0
            lambda_arr     (:,:) = 0.0
            pot_area_loss  (:,:) = 0.0
            act_area_loss  (:,:) = 0.0
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Save the original area, it will be used to correct the transition matrix !
            ! at the end of the sub-routine.                                               !
            !------------------------------------------------------------------------------!
            save_area: do ipa=1,onsp
               original_area(ipa) = csite%area     (ipa)
               original_lu  (ipa) = csite%dist_type(ipa)
            end do save_area
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Find the area to be harvested when biomass targets have been            !
            ! established.                                                                 !
            !------------------------------------------------------------------------------!
            call find_lambda_harvest(cpoly,isi,onsp,lambda_harvest)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Prepare the current site to receive the newly disturbed patches.  This   !
            ! step is done differently depending on whether this is a big leaf or a SAS    !
            ! simulation.                                                                  !
            !------------------------------------------------------------------------------!
            select case (ibigleaf)
            case (0)
               !---------------------------------------------------------------------------!
               !     Size-and-age structure.                                               !
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !     Transfer the original patch values to a temporary patch.  Then re-    !
               ! -allocate csite with room for new patches, and transfer the original      !
               ! patches back.                                                             !
               !---------------------------------------------------------------------------!
               call allocate_sitetype  (tsite,onsp)
               call copy_sitetype      (csite,tsite,1,onsp,1,onsp)
               call deallocate_sitetype(csite)
               call allocate_sitetype  (csite,onsp + n_dist_types)
               call copy_sitetype      (tsite,csite,1,onsp,1,onsp)
               call deallocate_sitetype(tsite)
               !---------------------------------------------------------------------------!


               !----- Allocate and initialise a disturbance mask vector. ------------------!
               allocate(disturb_mask(onsp + n_dist_types))
               disturb_mask(:)      = .false.
               disturb_mask(1:onsp) = .true.
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Initialize all the potential as well as implemented disturbance      !
               ! patches.  n_dist_types new patches will be created, each one containing a !
               ! different patch type.  In case no conversion to that kind of patch has    !
               ! happened, or if the newly created patch is tiny, it will be removed soon. !
               !---------------------------------------------------------------------------!
               init_dist_sas: do new_lu = onsp+1, onsp+n_dist_types
                  call initialize_disturbed_patch(csite,new_lu,cpoly%lsl(isi))
               end do init_dist_sas
               !---------------------------------------------------------------------------!

            case (1)
               !---------------------------------------------------------------------------!
               !     Big-leaf ED.                                                          !
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Transfer the original patch values to a temporary patch.  Then re-    !
               ! -allocate csite with room for new PFT/patches, and transfer the original  !
               ! patches back.                                                             !
               !---------------------------------------------------------------------------!
               call allocate_sitetype  (tsite,onsp)
               call copy_sitetype      (csite,tsite,1,onsp,1,onsp)
               call deallocate_sitetype(csite)
               call allocate_sitetype  (csite,onsp + nnsp_ble)
               call copy_sitetype      (tsite,csite,1,onsp,1,onsp)
               call deallocate_sitetype(tsite)
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Allocate the mask to decide which patches will remain by the end.     !
               !---------------------------------------------------------------------------!
               allocate(disturb_mask(onsp + nnsp_ble))
               disturb_mask         = .false.
               disturb_mask(1:onsp) = .true.
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Initialize all the potential as well as implemented disturbance      !
               ! patches.  n_dist_types new patches will be created, each one containing a !
               ! different patch type.  In case no conversion to that kind of patch has    !
               ! happened, or if the newly created patch is tiny, it will be removed soon. !
               !---------------------------------------------------------------------------!
               init_distpatch_ble: do ipa = onsp+1, onsp+nnsp_ble
                  call initialize_disturbed_patch(csite,ipa,cpoly%lsl(isi))
               end do init_distpatch_ble
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !       Reset mortality rates due to disturbance.                              !
            !------------------------------------------------------------------------------!
            resetdist: do ipa=1,onsp
               cpatch => csite%patch(ipa)
               do ico=1,cpatch%ncohorts
                  cpatch%mort_rate(6,ico) = 0.0
               end do
            end do resetdist
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      First iteration: we determine the area loss to each disturbance type,   !
            ! and fix the rates based on whether there is any biomass-based harvesting.    !
            ! We also ensure that the lost area is equivalent to the sum of disturbance    !
            ! rates.                                                                       !
            !------------------------------------------------------------------------------!
            old_lu_l1st: do ipa=1,onsp

               !---------------------------------------------------------------------------!
               !      Save the old land use in a shorter variable for convenience.         !
               ! "Managed" means a completely managed patch: pasture, cropland, and        !
               ! forest plantation.  Logged forests are not considered managed in the      !
               ! sense that anything can potentially grow in there.                        !
               !---------------------------------------------------------------------------!
               old_lu       = csite%dist_type(ipa)
               is_managed   = old_lu == 1 .or. old_lu == 2 .or. old_lu == 8
               is_oldgrowth = csite%age(ipa) >= min_oldgrowth(old_lu)
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !     Copy the disturbance rates to a temporary array.  We may need to      !
               ! adjust harvesting in case disturbance is based on demand.                 !
               !---------------------------------------------------------------------------!
               lambda_now   (:) = cpoly%disturbance_rates (:,old_lu,isi)                   &
                                + cpoly%disturbance_memory(:,old_lu,isi)
               one_area_loss(:) = 0.0
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Treefall disturbance cannot happen in case there is no tall tree.    !
               ! This is different from the background mortality that is truly density     !
               ! independent (treefall disturbance is in fact applied to small trees, but  !
               ! it does not generate gaps).                                               !
               !---------------------------------------------------------------------------!
               cpatch => csite%patch(ipa)
               if (cpatch%ncohorts == 0) then
                  !----- Empty patch, and we need trees to have tree fall. ----------------!
                  lambda_now(3) = 0.
                  !------------------------------------------------------------------------!
               else
                  if (cpatch%hite(1) < treefall_hite_threshold) then
                     !----- Tallest cohort is too short to create gaps. -------------------!
                     lambda_now(3) = 0.
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Update lambda in case there are biomass-based harvest demands for    !
               ! this patch.                                                               !
               !---------------------------------------------------------------------------!
               if (lambda_harvest(ipa) > 0) then
                  if (is_plantation) then
                     !----- Forest plantation. --------------------------------------------!
                     lambda_now(2) = lambda_harvest(ipa)
                     !---------------------------------------------------------------------!
                  else
                     !----- Logging.  We must also calculate damage. ----------------------!
                     lambda_now(6) = lambda_harvest(ipa)
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!
               else
                  !----- No harvest should occur. -----------------------------------------!
                  lambda_now(2) = 0.
                  lambda_now(6) = 0.
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Collateral damage due to logging.  This is defined as a function of   !
               ! the tree felling disturbance, and can be thought as a way to quantify     !
               ! the logging impact.                                                       !
               !---------------------------------------------------------------------------!
               if (sl_skid_rel_area == 1.0) then
                  !----- Special case: we simply copy the felling disturbance. ------------!
                  lambda_now(7) = lambda_now(6)
                  !------------------------------------------------------------------------!
               else if (lambda_now(6) > 0.0 .and. sl_skid_rel_area > 0.0) then
                  !------------------------------------------------------------------------!
                  !     The collateral damage is given in terms of area, so we must        !
                  ! convert the ratio to disturbance rate.  We use the definition of       !
                  ! disturbed area.                                                        !
                  !------------------------------------------------------------------------!
                  lfactor = sl_skid_rel_area * (1.0 - exp( - lambda_now(6)))
                  if ( lfactor <= 1.0 * (1.0 - epsilon(1.0)) ) then
                     lambda_now(7) = log( 1.0 / (1.0 - lfactor) )
                  else
                     lambda_now(7) = lnexp_max
                  end if
                  !------------------------------------------------------------------------!
               else
                  !----- No collateral damage, or no logging. -----------------------------!
                  lambda_now(7) = 0.0
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Total area loss depends on the sum of all disturbance rates.  The     !
               ! individual loss is then calculated to be proportional to each disturbance !
               ! rate if the rates were the only rates.  We skip this step in case all     !
               ! disturbance rates are zero.                                               !
               !---------------------------------------------------------------------------!
               if (any(lambda_now > 0.)) then
                  lambda_sum           = min(lnexp_max,sum(lambda_now))
                  area_loss_tot        = csite%area(ipa) * (1.0 - exp(-lambda_sum   ))
                  one_area_loss(:)     = csite%area(ipa) * (1.0 - exp(-lambda_now(:)))
                  pot_area_loss(ipa,:) = one_area_loss(:) * area_loss_tot                  &
                                       / sum(one_area_loss)
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Lambda array is used to show the transition matrix.                  !
               !---------------------------------------------------------------------------!
               lambda_arr(ipa,:) = lambda_now
               !---------------------------------------------------------------------------!
            end do old_lu_l1st
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Loop over the new LU patches, and check whether any of them would        !
            ! have enough area to become a patch.                                          !
            !------------------------------------------------------------------------------!
            new_lu_l2nd: do new_lu=1,n_dist_types
               pot_area_gain(new_lu) = sum(pot_area_loss(:,new_lu))
            end do new_lu_l2nd
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Loop over the old patches, and check whether any of them would become    !
            ! too small.  In case so, then we must adjust the area.                        !
            !------------------------------------------------------------------------------!
            old_lu_l2nd: do ipa=1,onsp
               !----- Only the LU transitions with enough area are actually created. ------!
               where (pot_area_gain(:) >= min_patch_area)
                  act_area_loss(ipa,:) = pot_area_loss(ipa,:)
               end where
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    In case the area would become smaller than the minimum area,  the      !
               ! patch must be eliminated. Set disturb_mask to false, so this patch is     !
               ! purged from csite later in this sub-routine.                              !
               !---------------------------------------------------------------------------!
               pot_area_remain = csite%area(ipa) - sum(act_area_loss(ipa,:))
               if ( pot_area_remain < min_patch_area ) then
                  disturb_mask (ipa)   = .false.
               end if
               !---------------------------------------------------------------------------!

            end do old_lu_l2nd
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Another loop over the new patches.  This time we define the patches      !
            ! that will be generated.                                                      !
            !------------------------------------------------------------------------------!
            new_lu_l3rd: do new_lu=1,n_dist_types
               !----- Retrieve the area gain for this new land use type. ------------------!
               act_area_gain(new_lu) = sum(act_area_loss(:,new_lu))
               cb_enthalpy_np        = 0.0
               can_exner_np          = 0.0
               cb_mass_np            = 0.0
               cb_molar_np           = 0.0
               cb_water_np           = 0.0
               cb_co2_np             = 0.0
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Populate the new patch in case there is sufficient disturbed area.   !
               !---------------------------------------------------------------------------!
               if (act_area_gain(new_lu) >= min_patch_area) then
                  write(unit=*,fmt='(a,1x,es12.5,1x,a,1x,i5)')                             &
                      ' ---> Making new patch, with area=',act_area_gain(new_lu)           &
                                             ,' for dist_type=',new_lu

                  !------------------------------------------------------------------------!
                  !    Fix flags for new patches.  Here it matters whether this simulation !
                  ! is big-leaf or SAS.                                                    !
                  !------------------------------------------------------------------------!
                  select case (ibigleaf)
                  case (0)
                     !---------------------------------------------------------------------!
                     !     Size-and-age structure.                                         !
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     Set the flag that this patch should be kept as a newly created  !
                     ! transition patch.                                                   !
                     !---------------------------------------------------------------------!
                     disturb_mask    (onsp+new_lu) = .true.
                     csite%dist_type (onsp+new_lu) = new_lu
                     csite%area      (onsp+new_lu) = act_area_gain(new_lu)
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     Initialize to zero the new trasitioned patches.                 !
                     !---------------------------------------------------------------------!
                     call initialize_disturbed_patch(csite,onsp+new_lu,cpoly%lsl(isi))
                     !---------------------------------------------------------------------!

                  case (1)
                     !---------------------------------------------------------------------!
                     !     Big leaf.                                                       !
                     !---------------------------------------------------------------------!
                     !---------------------------------------------------------------------!
                     !     Set the flag that this patch should be kept as a newly created  !
                     ! transition patch.                                                   !
                     !---------------------------------------------------------------------!
                     select case (new_lu)
                     case (1,2,8)
                        !------------------------------------------------------------------!
                        !     Cropland, pasture, and forest plantation.  Only one PFT is   !
                        ! allowed in such patches.                                         !
                        !------------------------------------------------------------------!
                        select case (new_lu)
                        case (1,2)
                           npa = onsp + new_lu
                        case (8)
                           npa = onsp + 2 + (new_lu - 3)*mypfts + 1
                        end select
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Set the flag that this patch should be kept as a newly       !
                        ! created transition patch.                                        !
                        !------------------------------------------------------------------!
                        disturb_mask    (npa)  = .true.
                        csite%dist_type (npa)  = new_lu
                        csite%area      (npa)  = act_area_gain(new_lu)
                        !------------------------------------------------------------------!



                        !----- Initialize to zero the new trasitioned patches. ------------!
                        call initialize_disturbed_patch(csite,onsp+new_lu,cpoly%lsl(isi))
                        !------------------------------------------------------------------!

                     case default
                        !------------------------------------------------------------------!
                        !     Non-cultivated lands.  Create one patch for each PFT that    !
                        ! may be present.                                                  !
                        !------------------------------------------------------------------!
                        pft_add_loop: do ipft=1,mypfts
                           !---------------------------------------------------------------!
                           !     Set the flag that this patch should be kept as a newly    !
                           ! created transition patch.                                     !
                           !---------------------------------------------------------------!
                           npa                   = onsp + 2 + (new_lu - 3) * mypfts + ipft
                           disturb_mask    (npa) = .true.
                           csite%dist_type (npa) = new_lu
                           csite%area      (npa) = act_area_gain(new_lu) / real(mypfts)
                           !---------------------------------------------------------------!


                           !----- Initialize to zero the new trasitioned patches. ---------!
                           call initialize_disturbed_patch(csite,npa,cpoly%lsl(isi))
                           !---------------------------------------------------------------!
                        end do pft_add_loop
                        !------------------------------------------------------------------!
                     end select
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Now go through patches, adding its contribution to the new          !
                  ! patch.                                                                 !
                  !------------------------------------------------------------------------!
                  old_lu_l3rd: do ipa=1,onsp
                     !----- Point to the current patch. -----------------------------------!
                     cpatch => csite%patch(ipa)
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !      Save the old land use in a shorter variable for convenience.   !
                     !---------------------------------------------------------------------!
                     old_lu       = csite%dist_type(ipa)
                     is_oldgrowth = csite%age(ipa) >= min_oldgrowth(old_lu)
                     is_primary   = old_lu == 3 .or. is_oldgrowth
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !    Check whether to apply disturbance to this patch, and assign     !
                     ! default values for disturbance pathway and minimum DBH for          !
                     ! harvesting.                                                         !
                     !---------------------------------------------------------------------!
                     disturbed               = act_area_loss(ipa,new_lu) > tiny_num
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !    If the patch is going to be disturbed, add the area and          !
                     ! survivors to the new patch.  The disturbed area has been already    !
                     ! deducted (check loop old_lu_l2nd), no need to change the old        !
                     ! patch area here.                                                    !
                     !---------------------------------------------------------------------!
                     if (disturbed) then

                        !------------------------------------------------------------------!
                        !     Add area and survivors to the new patches.  Here again we    !
                        ! must check whether this is a conventional or a  big leaf         !
                        ! simulation.                                                      !
                        !------------------------------------------------------------------!
                        select case (ibigleaf)
                        case (0)
                           !---------------------------------------------------------------!
                           !     Size-and-age structure.                                   !
                           !---------------------------------------------------------------!
                           area_fac = act_area_loss(ipa,new_lu) / csite%area(onsp+new_lu)
                           call increment_patch_vars(csite,onsp+new_lu,ipa,area_fac        &
                                                    ,cb_enthalpy_np,can_exner_np           &
                                                    ,cb_mass_np,cb_molar_np,cb_water_np    &
                                                    ,cb_co2_np)
                           call insert_survivors(csite,onsp+new_lu,ipa,new_lu,area_fac     &
                                                ,mindbh_harvest)
                           call accum_dist_harv_litt(cpoly,isi,1,onsp+new_lu,ipa,new_lu    &
                                                    ,area_fac,mindbh_harvest)
                           !---------------------------------------------------------------!
                        case (1)
                           !---------------------------------------------------------------!
                           !     Big leaf ED.  Here we must also check the type of patch   !
                           ! we are creating.                                              !
                           !---------------------------------------------------------------!
                           select case (new_lu)
                           case (1,2,8)
                              !------------------------------------------------------------!
                              !    Find the correct index.                                 !
                              !------------------------------------------------------------!
                              select case (new_lu)
                              case (1,2)
                                 npa = onsp + new_lu 
                              case (8)
                                 npa = onsp + 2 + (new_lu - 3) * mypfts + 1
                              end select
                              !------------------------------------------------------------!

                              !------------------------------------------------------------!
                              !     Cropland, pasture, or forest plantation.               !
                              !------------------------------------------------------------!
                              area_fac = act_area_loss(ipa,new_lu)/csite%area(npa)
                              call increment_patch_vars(csite,npa,ipa,area_fac             &
                                                       ,cb_enthalpy_np,can_exner_np        &
                                                       ,cb_mass_np,cb_molar_np,cb_water_np &
                                                       ,cb_co2_np)
                              call insert_survivors(csite,npa,ipa,new_lu,area_fac          &
                                                   ,mindbh_harvest)
                              call accum_dist_harv_litt(cpoly,isi,1,npa,ipa,new_lu         &
                                                       ,area_fac,mindbh_harvest)
                              !------------------------------------------------------------!
                           case default
                              !------------------------------------------------------------!
                              !     Non-cultivated lands.  We must check that the donor    !
                              ! patch and receptor patch have the same PFT, or the donor   !
                              ! is empty.                                                  !
                              !------------------------------------------------------------!
                              do ipft=1,mypfts
                                 npa = onsp + 2 + (new_lu - 3) * mypfts + ipft

                                 if (cpatch%ncohorts == 0) then
                                    !------------------------------------------------------!
                                    !    Donor cohort is empty, split disturbed area       !
                                    ! evenly amongst receptor cohorts.                     !
                                    !------------------------------------------------------!
                                    same_pft = ipft == 1
                                    area_fac = act_area_loss(ipa,new_lu)                   &
                                             / ( real(mypfts) * csite%area(npa) )
                                    !------------------------------------------------------!
                                 else
                                    if (ipft == cpatch%pft(1)) then
                                       !---------------------------------------------------!
                                       !    PFTs match, send all disturbed area to this    !
                                       ! patch.                                            !
                                       !---------------------------------------------------!
                                       same_pft = .true.
                                       area_fac = act_area_loss(ipa,new_lu)/csite%area(npa)
                                       !---------------------------------------------------!
                                    else
                                       !---------------------------------------------------!
                                       !    PFTs do not match, don't send anything to the  !
                                       ! new patch.                                        !
                                       !---------------------------------------------------!
                                       same_pft = .false.
                                       area_fac = 0.
                                       !---------------------------------------------------!
                                    end if
                                    !------------------------------------------------------!
                                 end if
                                 !---------------------------------------------------------!


                                 !---------------------------------------------------------!
                                 !     Accumulate survivors to the new patch.              !
                                 !---------------------------------------------------------!
                                 if (same_pft) then
                                    call increment_patch_vars(csite,npa,ipa,area_fac       &
                                                             ,cb_enthalpy_np,can_exner_np  &
                                                             ,cb_mass_np,cb_molar_np       &
                                                             ,cb_water_np,cb_co2_np)
                                    call insert_survivors(csite,npa,ipa,new_lu,area_fac    &
                                                         ,mindbh_harvest)
                                    call accum_dist_harv_litt(cpoly,isi,1,npa,ipa,new_lu   &
                                                             ,area_fac,mindbh_harvest)
                                 end if
                                 !---------------------------------------------------------!
                              end do
                              !------------------------------------------------------------!
                           end select
                           !---------------------------------------------------------------!
                        end select
                        !------------------------------------------------------------------!
                     end if
                     !---------------------------------------------------------------------!
                  end do old_lu_l3rd
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Normalise canopy air properties.  This must be done after         !
                  ! incrementing canopy air space variables.                               !
                  !------------------------------------------------------------------------!
                  call norm_canopy_air_props(csite,onsp+new_lu,cb_enthalpy_np,can_exner_np &
                                            ,cb_mass_np,cb_molar_np,cb_water_np,cb_co2_np)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Update temperature and density.  This must be done before         !
                  ! planting, because leaf temperature is initially assigned based on the  !
                  ! canopy air temperature.                                                !
                  !------------------------------------------------------------------------!
                  call update_patch_thermo_props(csite,onsp+new_lu,onsp+new_lu,nzg,nzs     &
                                                ,cpoly%ntext_soil(:,isi))
                  call update_patch_thermo_fmean(csite,onsp+new_lu,onsp+new_lu,nzg         &
                                                ,cpoly%ntext_soil(:,isi))
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !     In case the new patch is agriculture or a forest plantation, add   !
                  ! the new grasses/trees.                                                 !
                  !------------------------------------------------------------------------!
                  select case (new_lu)
                  case (1)
                     call plant_patch(csite,onsp+new_lu,nzg                                &
                                     ,cpoly%pasture_stocking_pft(isi)                      &
                                     ,cpoly%pasture_stocking_density(isi)                  &
                                     ,cpoly%ntext_soil(:,isi), 1.0                         &
                                     ,cpoly%lsl(isi))
                  case (2)
                     call plant_patch(csite,onsp+new_lu,nzg                                &
                                     ,cpoly%plantation_stocking_pft(isi)                   &
                                     ,cpoly%plantation_stocking_density(isi)               &
                                     ,cpoly%ntext_soil(:,isi), 2.0                         &
                                     ,cpoly%lsl(isi))
                  case (8)
                     call plant_patch(csite,onsp+new_lu,nzg                                &
                                     ,cpoly%agri_stocking_pft(isi)                         &
                                     ,cpoly%agri_stocking_density(isi)                     &
                                     ,cpoly%ntext_soil(:,isi), 1.0                         &
                                     ,cpoly%lsl(isi))
                  end select
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Update the cohorts.  This step is done only in case this is a size- !
                  ! -and-age structure.                                                    !
                  !------------------------------------------------------------------------!
                  qpatch => csite%patch(onsp+new_lu)
                  bfus_ncoh = qpatch%ncohorts
                  if (ibigleaf == 0 .and. qpatch%ncohorts > 0 .and. maxcohort >= 0) then
                     call sort_cohorts(qpatch)
                     select case (ifusion)
                     case (0)
                        call old_fuse_cohorts(csite,onsp+new_lu,cpoly%lsl(isi),.false.)
                     case (1)
                        call new_fuse_cohorts(csite,onsp+new_lu,cpoly%lsl(isi),.false.)
                     end select
                     call terminate_cohorts(csite,onsp+new_lu,cmet,.false.                 &
                                           ,elim_nplant,elim_lai)
                     call split_cohorts(csite,onsp+new_lu, cpoly%green_leaf_factor(:,isi)  &
                                       ,.false.)
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Prune lianas in case they are included in the simulation.          !
                  ! MLO -> Manfredo, this used to be outside the loop, I moved into here   !
                  !        so we don't update derived properties twice, and to make it     !
                  !        easier to track budget.  However, I wonder if this should be    !
                  !        done before the cohort sorting, fusion, fission and extinction  !
                  !        just above this comment.                                        !
                  !------------------------------------------------------------------------!
                  if (lianas_included) then
                     !----------------------- Prune the lianas ----------------------------!
                     call prune_lianas(csite, onsp + new_lu, cpoly%lsl(isi))
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !     Update the derived properties including veg_height, and patch-     !
                  ! -level LAI, WAI.                                                       !
                  !------------------------------------------------------------------------!
                  call update_patch_derived_props(csite,onsp+new_lu,.false.)
                  !----- Update soil temperature, liquid fraction, etc. -------------------!
                  call new_patch_sfc_props(csite,onsp+new_lu,nzg,nzs                       &
                                          ,cpoly%ntext_soil(:,isi))
                 !----- Update AGB, basal area. ------------------------------------------!
                  call update_site_derived_props(cpoly,1,isi)
                  !------------------------------------------------------------------------!


                  !----- Clear the disturbance memory for this disturbance type. ----------!
                  cpoly%disturbance_memory(new_lu,1:n_dist_types,isi) = 0.0
                  !------------------------------------------------------------------------!

               elseif (pot_area_gain(new_lu) > 0.0)then
                  !------------------------------------------------------------------------!
                  !     The patch creation has been skipped because the area was too       !
                  ! small.  Put the current disturbance rates in memory to be added at     !
                  ! the next timestep.                                                     !
                  !------------------------------------------------------------------------!
                  cpoly%disturbance_memory(new_lu,1:n_dist_types,isi) =                    &
                         cpoly%disturbance_memory(new_lu,1:n_dist_types,isi)               &
                       + cpoly%disturbance_rates (new_lu,1:n_dist_types,isi)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do new_lu_l3rd
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Print the entire disturbance transition matrix if debugging mode is on.  !
            !------------------------------------------------------------------------------!
            if (print_detailed) then
               write(unit=*,fmt='(a)')          ' '
               write(unit=*,fmt='(103a)')       ('-',i=1,103)
               write(unit=*,fmt='(2(a,1x,i5))') ' Summary for IPY =',ipy,'; ISI =',isi
               write(unit=*,fmt='(a)')          ' '
               write(unit=*,fmt='(12(1x,a))')   '  IPA','  ILU'                            &
                                              ,'    AREA','  TO_PST','  TO_FPL','  TO_TFL' &
                                              ,'  TO_BRN','  TO_ABN','  TO_LOG','  TO_SKD' &
                                              ,'  TO_CPL','TOT_LOSS'
               write(unit=*,fmt='(103a)')       ('-',i=1,103)
               do ipa=1,onsp
                  write(unit=*,fmt='(2(1x,i5),10(1x,f8.5))')                               &
                                                  ipa,csite%dist_type(ipa),csite%area(ipa) &
                                                   ,(act_area_loss(ipa,new_lu),new_lu=1,8) &
                                                   ,sum(act_area_loss(ipa,:))
               end do
               write(unit=*,fmt='(103a)')       ('-',i=1,103)
               write(unit=*,fmt='(1x,a,10(1x,f8.5))')               '            TOT_GAIN' &
                                                       ,(act_area_gain(new_lu),new_lu=1,8) &
                                                       ,sum(act_area_gain(:))
               write(unit=*,fmt='(103a)')       ('-',i=1,103)
               write(unit=*,fmt='(a)')          ' Transition matrix (disturbance rates)'
               write(unit=*,fmt='(a)')          ' '
               write(unit=*,fmt='(11(1x,a))')   '  IPA','  ILU'                            &
                                              ,'    AREA','  TO_PST','  TO_FPL','  TO_TFL' &
                                              ,'  TO_BRN','  TO_ABN','  TO_LOG','  TO_SKD' &
                                              ,'  TO_CPL'
               write(unit=*,fmt='(94a)')        ('-',i=1,94)
               do ipa=1,onsp
                  write(unit=*,fmt='(2(1x,i5),9(1x,f8.5))')                                &
                                                  ipa,csite%dist_type(ipa),csite%area(ipa) &
                                                   ,(lambda_arr(ipa,new_lu),new_lu=1,8)
               end do
               write(unit=*,fmt='(94a)')        ('-',i=1,94)
               write(unit=*,fmt='(a)')    ' '
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Update mortality rate due to disturbance, and deduct the distubed area  !
            ! from this patch.                                                             !
            !------------------------------------------------------------------------------!
            old_lu_l4th: do ipa=1,onsp
               pat_area_loss = act_area_loss(ipa,:)
               call disturbance_mortality(csite,ipa,pat_area_loss,mindbh_harvest)
               csite%area(ipa) = csite%area(ipa) - sum(pat_area_loss)
            end do old_lu_l4th
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Big-leaf only.  In BigLeaf-ED, we do not track age since last            !
            ! disturbance; instead, each land use type has one patch for each PFT or one   !
            ! patch and one PFT only in case the patch is cultivated.  Here, we search for !
            ! duplicated patches that have the same PFT and same LU, and fuse them.        !
            !------------------------------------------------------------------------------!
            select case (ibigleaf)
            case (1)
               merge_lu_loop: do new_lu = 1, n_dist_types
                  !----- Define search limits for this LU. --------------------------------!
                  select case (new_lu)
                  case (1)
                     apa = onsp + 1
                     zpa = onsp + 1
                  case (2)
                     apa = onsp + 2
                     zpa = onsp + 2
                  case (8)
                     apa = onsp + 2 + ( new_lu - 3 ) * mypfts + 1
                     zpa = onsp + 2 + ( new_lu - 3 ) * mypfts + 1
                  case default
                     apa = onsp + 2 + ( new_lu - 3 ) * mypfts + 1
                     zpa = onsp + 2 + ( new_lu - 3 ) * mypfts + mypfts
                  end select
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !     Loop over the old patches and seek land use matches.               !
                  !------------------------------------------------------------------------!
                  merge_patch_loop: do npa=apa,zpa

                     !----- Check whether the patch has been created. ---------------------!
                     if (disturb_mask(npa)) then

                        !------------------------------------------------------------------!
                        !     Loop over old patches.                                       !
                        !------------------------------------------------------------------!
                        old_lu_search: do ipa = 1,onsp
                           old_lu = csite%dist_type(ipa)

                           !---------------------------------------------------------------!
                           !    If this is a good receptor, merge the patches.             !
                           !                                                               !
                           !    Requirements for a good receptor:                          !
                           ! 1. Same land use                                              !
                           ! 2. Receptor is not slated to be purged.                       !
                           ! 3. Receptor has the same PFT as the donor (or both patches    !
                           !    are empty).                                                !
                           !---------------------------------------------------------------!
                           if     ( csite%patch(ipa)%ncohorts == 0 .and.                   &
                                    csite%patch(npa)%ncohorts == 0 ) then
                              same_pft = .true.
                           elseif ( csite%patch(ipa)%ncohorts == 0 .or.                    &
                                    csite%patch(npa)%ncohorts == 0 ) then
                              same_pft = .false.
                           else
                              same_pft = csite%patch(ipa)%pft(1) == csite%patch(npa)%pft(1)
                           end if

                           if ( old_lu == new_lu .and. disturb_mask(ipa) .and. same_pft)   &
                           then
                              !----- Fuse both patches. -----------------------------------!
                              call fuse_2_patches(csite,npa,ipa,nzg,nzs                    &
                                                 ,cpoly%met(isi)                           &
                                                 ,cpoly%lsl(isi)                           &
                                                 ,cpoly%ntext_soil(:,isi)                  &
                                                 ,cpoly%green_leaf_factor(:,isi)           &
                                                 ,.false.,elim_nplant,elim_lai)
                              !------------------------------------------------------------!



                              !------------------------------------------------------------!
                              !    Make sure that all cohorts are fused.                   !
                              !------------------------------------------------------------!
                              cpatch => csite%patch(ipa)
                              do ico=2,cpatch%ncohorts
                                 new_nplant = cpatch%nplant(ico) + cpatch%nplant(1)
                                 ipft       = cpatch%pft(1)

                                 !---------------------------------------------------------!
                                 !    Temporarily make both cohorts "resolvable".  This is !
                                 ! to avoid energy/water leaks when a non-resolvable       !
                                 ! cohort is fused with a resolvable cohort.               !
                                 !---------------------------------------------------------!
                                 call is_resolvable(csite,ipa,ico,.false.,.true.           &
                                                   ,'apply_disturbances (ico,BLE,before)')
                                 call is_resolvable(csite,ipa,1,.false.,.true.             &
                                                   ,'apply_disturbances (1,BLE,before)')
                                 !---------------------------------------------------------!


                                 !----- Fuse cohorts. -------------------------------------!
                                 call fuse_2_cohorts(cpatch,ico,1, csite%can_prss(ipa)     &
                                                    ,csite%can_shv (ipa),cpoly%lsl(isi)    &
                                                    ,.false.)
                                 !---------------------------------------------------------!



                                 !---------------------------------------------------------!
                                 !     Check whether the final fused cohort is resolvable. !
                                 ! In case it is not, then we must subtract the phenology  !
                                 ! effect that was temporarily added before fusing the     !
                                 ! cohorts.  Unlike the "before" calls, here we do not     !
                                 ! force cohorts to be resolvable.                         !
                                 !---------------------------------------------------------!
                                 call is_resolvable(csite,ipa,1,.false.,.false.            &
                                                   ,'apply_disturbances (1,BLE,after)')
                                 !---------------------------------------------------------!



                                 !---------------------------------------------------------!
                                 !     Set nplant to a tiny number, we will delete this    !
                                 ! cohort soon.                                            !
                                 !---------------------------------------------------------!
                                 cpatch%nplant(ico) = tiny_num
                                 !---------------------------------------------------------!
                              end do
                              !------------------------------------------------------------!


                              !------ Remove emptied cohorts. -----------------------------!
                              call terminate_cohorts(csite,ipa,cmet,.false.                &
                                                    ,elim_nplant,elim_lai)
                              !------------------------------------------------------------!


                              !------ Remove the patch. -----------------------------------!
                              disturb_mask(npa) = .false.
                              !------------------------------------------------------------!


                              !------ Patch is gone, stop searching old patches. ----------!
                              exit old_lu_search
                              !------------------------------------------------------------!
                           end if
                           !---------------------------------------------------------------!
                        end do old_lu_search
                        !------------------------------------------------------------------!
                     end if
                     !---------------------------------------------------------------------!
                  end do merge_patch_loop
                  !------------------------------------------------------------------------!
               end do merge_lu_loop
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !      Reallocate the current site to fit the original patches, except         !
            ! those that lost all area through disturbance, and the new disturbed          !
            ! patches that were actually generated.                                        !
            !------------------------------------------------------------------------------!
            !----- Temporary site, with size equal to the patches that will remain. -------!
            call allocate_sitetype(tsite,count(disturb_mask))
            call copy_sitetype_mask(csite,tsite,disturb_mask,size(disturb_mask)            &
                                   ,count(disturb_mask))
            call deallocate_sitetype(csite)
            !----- Re-allocate the tracked site and copy data from temporary site. --------!
            call allocate_sitetype(csite,tsite%npatches)
            call copy_sitetype(tsite,csite,1,tsite%npatches,1,tsite%npatches)
            call deallocate_sitetype(tsite)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Update the disturbance rates to be consistent with the applied          !
            ! transitions.  This is necessary because some transitions are not calculated  !
            ! directly from disturbance_rates (for example, biomass-based harvesting), and !
            ! some exceptions may occur (for example, when previous memory has been        !
            ! applied, or when there wasn't enough area to create a new patch or when      !
            ! disturbance was extended to the entire patch because the remaining area      !
            ! would be too small otherwise.                                                !
            !------------------------------------------------------------------------------!
            !----- Print the area loss and disturbance transitions. -----------------------!
            if (print_detailed) then
               write(unit=*,fmt='(a)')          ' '
               write(unit=*,fmt='(45a)')        ('-',i=1,45)
               write(unit=*,fmt='(2(a,1x,i5))') ' Transitions for IPY =',ipy,'; ISI =',isi
               write(unit=*,fmt='(a)')          ' '
               write(unit=*,fmt='(5(1x,a))')    'OLD_LU','NEW_LU'                          &
                                               ,'ORIG_AREA','AREA_LOSS','DIST_RATE'
               write(unit=*,fmt='(45a)')        ('-',i=1,45)
            end if
            !------------------------------------------------------------------------------!

            dist_from_loop: do old_lu=1,n_dist_types
               dist_to_loop: do new_lu=1,n_dist_types

                  !----- Reset areas for this transtion (old_lu -> new_lu). ---------------!
                  orig_area = 0.0
                  dist_area = 0.0
                  dist_orig_loop: do ipa = 1, onsp
                     if (original_lu(ipa) == old_lu .and. act_area_loss(ipa,new_lu) > 0.)  &
                     then
                        orig_area = orig_area + original_area(ipa)
                        dist_area = dist_area + act_area_loss(ipa,new_lu)
                     end if
                  end do dist_orig_loop
                  !------------------------------------------------------------------------!


                  !----- Update the disturbance rate. -------------------------------------!
                  if (orig_area == 0.0) then
                     cpoly%disturbance_rates(new_lu,old_lu,isi) = 0.
                  elseif (dist_area <  orig_area ) then
                     cpoly%disturbance_rates(new_lu,old_lu,isi) =                          &
                                                log( orig_area / ( orig_area - dist_area ))
                  else
                     cpoly%disturbance_rates(new_lu,old_lu,isi) = lnexp_max
                  end if
                  !------------------------------------------------------------------------!



                  !----- Print the area loss and disturbance transitions. -----------------!
                  if (print_detailed) then
                     write(unit=*,fmt='(2(1x,i6),3(1x,f9.6))')                             &
                                                 old_lu,new_lu,orig_area,dist_area         &
                                                ,cpoly%disturbance_rates(new_lu,old_lu,isi)
                  end if
                  !------------------------------------------------------------------------!
               end do dist_to_loop
               !---------------------------------------------------------------------------!
            end do dist_from_loop
            !------------------------------------------------------------------------------!




            !----- Print the area loss and disturbance transitions. -----------------------!
            if (print_detailed) then
               write(unit=*,fmt='(45a)')        ('-',i=1,45)
               write(unit=*,fmt='(a)')          ' '
            end if
            !------------------------------------------------------------------------------!


            !------ Update polygon-level crop and logging harvest. ------------------------!
            cgrid%crop_harvest   (ipy) = cgrid%crop_harvest   (ipy)                        &
                                       + cpoly%crop_harvest   (isi) * cpoly%area(isi)
            cgrid%logging_harvest(ipy) = cgrid%logging_harvest(ipy)                        &
                                       + cpoly%logging_harvest(isi) * cpoly%area(isi)
            cgrid%combusted_fuel (ipy) = cgrid%combusted_fuel (ipy)                        &
                                       + cpoly%combusted_fuel (isi) * cpoly%area(isi)
            !------------------------------------------------------------------------------!



            !----- Free memory before re-allocating for the next site... ------------------!
            deallocate(disturb_mask   )
            deallocate(original_area  )
            deallocate(original_lu    )
            deallocate(lambda_harvest )
            deallocate(lambda_arr     )
            deallocate(pot_area_loss  )
            deallocate(act_area_loss  )
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!



         !----- Save the total carbon loss through removal for the C budget bookkeeping. --!
         cgrid%cbudget_removedstorage(ipy) = cgrid%crop_harvest   (ipy)                    &
                                           + cgrid%logging_harvest(ipy)                    &
                                           + cgrid%combusted_fuel (ipy)
         !---------------------------------------------------------------------------------!


      end do polyloop
      !------------------------------------------------------------------------------------!

      !----- Free memory before leaving... ------------------------------------------------!
      deallocate(tsite)
      deallocate(pfts )
      !------------------------------------------------------------------------------------!

      return
   end subroutine apply_disturbances
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine fills in the disturbance transition matrix.                      !
   !---------------------------------------------------------------------------------------!
   subroutine site_disturbance_rates(year, cgrid)

      use ed_state_vars, only : edtype                    & ! structure
                              , polygontype               & ! structure
                              , sitetype                  & ! structure
                              , patchtype                 ! ! structure
      use disturb_coms , only : treefall_disturbance_rate & ! intent(in)
                              , lutime                    & ! structure
                              , ianth_disturb             & ! intent(in)
                              , include_fire              & ! intent(in)
                              , plantation_year           & ! intent(in)
                              , plantation_rotation       & ! intent(in)
                              , mature_harvest_age        & ! intent(in)
                              , min_oldgrowth             & ! intent(in)
                              , sl_biomass_harvest        ! ! intent(in)
      use ed_max_dims  , only : n_pft                     & ! intent(in)
                              , n_dist_types              ! ! intent(in)
      use ed_misc_coms , only : current_time              ! ! intent(in)
      use consts_coms  , only : lnexp_max                 ! ! intent(in)

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target                  :: cgrid
      integer          , intent(in)              :: year
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer                 :: cpoly
      type(sitetype)   , pointer                 :: csite
      type(patchtype)  , pointer                 :: cpatch
      type(lutime)     , pointer                 :: clutime
      integer                                    :: ipy
      integer                                    :: isi
      integer                                    :: ipa
      integer                                    :: ico
      integer                                    :: ipft
      integer                                    :: ilu
      integer                                    :: iyear
      integer                                    :: useyear
      real                                       :: bharvest
      real                                       :: pot_harvest_target
      real                                       :: fire_disturbance_rate
      logical                                    :: is_plantation
      logical                                    :: is_rotation
      logical                                    :: is_mature
      logical                                    :: is_oldgrowth
      logical                                    :: harv_secondary
      logical                                    :: harv_primary
      logical                                    :: find_target
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Loop over sites and polygons.                                                 !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !----- Test whether the new explored patch will be plantation or logged. ------!
            is_plantation = cpoly%plantation(isi) == 1               .and.                 &
                            current_time%year     >  plantation_year
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !      Reset the disturbance rates.  We only populate the transitions that are !
            ! may be positive.  The first index is the new land use type, and the second   !
            ! index is the old land use type.  Both use the following convention:          !
            !                                                                              !
            !  1.  Pasture                                                                 !
            !  2.  Forest plantation                                                       !
            !  3.  Tree fall                                                               !
            !  4.  Burnt                                                                   !
            !  5.  Abandoned                                                               !
            !  6.  Logged (felling)                                                        !
            !  7.  Logged (skid trail / road)                                              !
            !  8.  Agricultural lands (cropland / pasture)                                 !
            !------------------------------------------------------------------------------!
            cpoly%disturbance_rates(:,:,isi) = 0.0
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !    Calculate fire disturbance rates only when fire is on.                    !
            !------------------------------------------------------------------------------!
            select case (include_fire)
            case (0)
               fire_disturbance_rate = 0.0
            case default
               fire_disturbance_rate = sum(cpoly%lambda_fire(1:12,isi)) / 12.0
            end select
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Tree fall and fires also occur in plantations.  Treefall is             !
            ! incorporated to mortality, whereas fire is assumed to cause abandonment.     !
            ! Not ideal, but this is the only way out of a forest plantation right now.    !
            !------------------------------------------------------------------------------!
            cpoly%disturbance_rates(5,2,isi) = fire_disturbance_rate
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Disturbance that creates new "tree fall patches".  Only non-cultivated  !
            ! lands may suffer this disturbance.                                           !
            !------------------------------------------------------------------------------!
            cpoly%disturbance_rates(3,3,isi) = treefall_disturbance_rate
            cpoly%disturbance_rates(3,4,isi) = treefall_disturbance_rate
            cpoly%disturbance_rates(3,5,isi) = treefall_disturbance_rate
            cpoly%disturbance_rates(3,6,isi) = treefall_disturbance_rate
            cpoly%disturbance_rates(3,7,isi) = treefall_disturbance_rate
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Disturbance that creates new "burnt patches".  Only non-cultivated      !
            ! lands may suffer this disturbance.                                           !
            !------------------------------------------------------------------------------!
            cpoly%disturbance_rates(4,3,isi) = fire_disturbance_rate
            cpoly%disturbance_rates(4,4,isi) = fire_disturbance_rate
            cpoly%disturbance_rates(4,5,isi) = fire_disturbance_rate
            cpoly%disturbance_rates(4,6,isi) = fire_disturbance_rate
            cpoly%disturbance_rates(4,7,isi) = fire_disturbance_rate
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Check for anthropogenic disturbance rates.  We no longer need  to check !
            ! the year because all years in this simulation are assigned prescribed        !
            ! disturbance rates (even if that means zero disturbance).                     !
            !------------------------------------------------------------------------------!
            useyear = cpoly%num_landuse_years(isi)
            !----- Loop over years. -------------------------------------------------------!
            find_lu_year: do iyear = 1,cpoly%num_landuse_years(isi)
               if (year == cpoly%clutimes(iyear,isi)%landuse_year) then
                  useyear = iyear
                  exit find_lu_year
               end if
            end do find_lu_year
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Link to the current land use change year.                                 !
            !------------------------------------------------------------------------------!
            clutime => cpoly%clutimes(useyear,isi)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Transition rates are based on the input file, formatted as the original !
            ! data set from George Hurtt et al. (2006).  This will be eventually replaced  !
            ! by the most current version.                                                 !
            !                                                                              !
            ! Hurtt, G. et al., 2006: The underpinnings of land-use history: three         !
            !    centuries of global gridded land-use transitions, wood-harvest activity,  !
            !    and resulting secondary forests.  Glob. Change Biol., 12(7), 1208-1229,   !
            !    doi:10.1111/j.1365-2486.2006.01150.x                                      !
            !                                                                              !
            ! The elements of vector clutime%landuse are:                                  !
            ! ..... Disturbance rates .................................................... !
            !   1 - Cropland to pasture                                            [1/yr]  !
            !   2 - Pasture to cropland                                            [1/yr]  !
            !   3 - Pasture to primary forest (in ED this goes to abandonment)     [1/yr]  !
            !   4 - Primary forest to pasture                                      [1/yr]  !
            !   5 - Primary forest to cropland                                     [1/yr]  !
            !   6 - Cropland to primary forest (in ED this goes to abandonment)    [1/yr]  !
            !   7 - Secondary forest to cropland                                   [1/yr]  !
            !   8 - Cropland to secondary forest (in ED this goes to abandonment)  [1/yr]  !
            !   9 - Secondary forest to pasture                                    [1/yr]  !
            !  10 - Pasture to secondary forest (in ED this goes to abandonment)   [1/yr]  !
            !  11 - Primary forest to secondary forest (assumed logging in ED)     [1/yr]  !
            ! ..... Harvested area ....................................................... !
            !  12 - Wood harvest on mature secondary forest land          [          kgC]  !
            !  13 - Wood harvest on mature secondary forest land          [grid fraction]  !
            !  14 - Wood harvest on primary forest land                   [          kgC]  !
            !  15 - Wood harvest on primary forest land                   [grid fraction]  !
            !  16 - Wood harvest on young secondary forest land           [          kgC]  !
            !  17 - Wood harvest on young secondary forest land           [grid fraction]  !
            !  18 - Wood harvest on primary non-forest land               [          kgC]  !
            !  19 - Wood harvest on primary non-forest land               [grid fraction]  !
            ! ............................................................................ !
            !                                                                              !
            !   The definition of primary and secondary forests can be quite confusing.    !
            ! For ED purposes, we assume that patches previously disturbed by treefall     !
            ! are primary, whereas patches previously generated by fire, abandonment, and  !
            ! logging are considered secondary.  Fire is the one that may be the hardest   !
            ! to reconcile: in tropical forests fires are almost exclusively               !
            ! anthropogenic, but they can be mostly natural in boreal forests, for         !
            ! example.                                                                     !
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Apply disturbance rates that generate pasture (see table above).        !
            !------------------------------------------------------------------------------!
            cpoly%disturbance_rates(1,2,isi) = clutime%landuse(9)
            cpoly%disturbance_rates(1,3,isi) = clutime%landuse(4)
            cpoly%disturbance_rates(1,4,isi) = clutime%landuse(9)
            cpoly%disturbance_rates(1,5,isi) = clutime%landuse(9)
            cpoly%disturbance_rates(1,6,isi) = clutime%landuse(9)
            cpoly%disturbance_rates(1,7,isi) = clutime%landuse(9)
            cpoly%disturbance_rates(1,8,isi) = clutime%landuse(1)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      Apply disturbance rates that generate cropland (see table above).        !
            !------------------------------------------------------------------------------!
            cpoly%disturbance_rates(8,1,isi) = clutime%landuse(2)
            cpoly%disturbance_rates(8,2,isi) = clutime%landuse(7)
            cpoly%disturbance_rates(8,3,isi) = clutime%landuse(5)
            cpoly%disturbance_rates(8,4,isi) = clutime%landuse(7)
            cpoly%disturbance_rates(8,5,isi) = clutime%landuse(7)
            cpoly%disturbance_rates(8,6,isi) = clutime%landuse(7)
            cpoly%disturbance_rates(8,7,isi) = clutime%landuse(7)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Secondary forests.   Here the transition depends on whether to establish !
            ! plantations or to assume secondary forests = abandoned lands.                !
            !------------------------------------------------------------------------------!
            select case (cpoly%plantation(isi))
            case (0)
               !----- Abandoned lands only. -----------------------------------------------!
               cpoly%disturbance_rates(5,1,isi) = clutime%landuse( 3) + clutime%landuse(10)
               cpoly%disturbance_rates(5,8,isi) = clutime%landuse( 6) + clutime%landuse( 8)
               !---------------------------------------------------------------------------!
            case (1)
               !----- "Secondary" sends area to plantation, "Primary" to abandoned. -------!
               cpoly%disturbance_rates(2,1,isi) = clutime%landuse(10)
               cpoly%disturbance_rates(2,8,isi) = clutime%landuse( 8)
               cpoly%disturbance_rates(5,1,isi) = clutime%landuse( 3)
               cpoly%disturbance_rates(5,8,isi) = clutime%landuse( 6)
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Harvesting (either plantation -> plantation or logging) when a biomass   !
            ! target does not exist (e.g. SimAmazonia)).  Convert the harvest probability  !
            ! of being cut given that the DBH exceeds the minimum DBH.  This is done only  !
            ! when anthropogenic disturbance is on and we are not seeking the biomass      !
            ! target, otherwise we set it to zero.                                         !
            !------------------------------------------------------------------------------!
            select case (ianth_disturb)
            case (0)
               !----- Anthropogenic disturbance is off. No logging should occur. ----------!
               find_target                         = .false.
               cpoly%primary_harvest_target  (isi) = 0.
               cpoly%secondary_harvest_target(isi) = 0.
               !---------------------------------------------------------------------------!
            case default
               !------ Read anthropogenic disturbance from external data set. -------------!
               if (clutime%landuse(12) < 0 .or. clutime%landuse(14) < 0) then
                  find_target                         = .true.
                  cpoly%primary_harvest_target  (isi) = 0.
                  cpoly%secondary_harvest_target(isi) = 0.
               else
                  find_target                         = .false.
                  cpoly%primary_harvest_target  (isi) = clutime%landuse(14)                &
                                                      + clutime%landuse(18)
                  cpoly%secondary_harvest_target(isi) = clutime%landuse(12)                &
                                                      + clutime%landuse(16)
               end if
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find target biomass as all the biomass that can be harvested.            !
            !------------------------------------------------------------------------------!
            if (find_target) then
               pharv_loop: do ipa=1,csite%npatches
                  cpatch => csite%patch(ipa)


                  !----- Check age and type of forest. ------------------------------------!
                  ilu          = csite%dist_type(ipa)
                  is_rotation  = is_plantation .and. csite%age(ipa) >= plantation_rotation
                  is_mature    = csite%age(ipa) >= mature_harvest_age
                  is_oldgrowth = csite%age(ipa) >= min_oldgrowth (ilu)
                  !------------------------------------------------------------------------!


                  !----- Reset harvest contribution from patch. ---------------------------!
                  bharvest     = 0.
                  !------------------------------------------------------------------------!



                  !----- Select the minimum DBH depending on the forest category. ---------!
                  select case(ilu)
                  case (1,8)
                     !----- Pasture and agriculture.  No timber harvesting here. ----------!
                     harv_primary      = .false.
                     harv_secondary    = .false.
                     !---------------------------------------------------------------------!
                  case (2)
                     !----- Forest plantation.  Usually all biomass is cleared. -----------!
                     harv_primary      = .false.
                     harv_secondary    = is_rotation
                     !---------------------------------------------------------------------!
                  case (3)
                     !----- Treefall.  Assume primary vegetation. -------------------------!
                     harv_primary      = is_mature
                     harv_secondary    = .false.
                     !---------------------------------------------------------------------!
                  case (4:7)
                     !---------------------------------------------------------------------!
                     !      Other disturbances.  Assume "primary" in case this is old-     !
                     ! growth, otherwise assume it secondary vegetation.                   !
                     !---------------------------------------------------------------------!
                     harv_primary   = is_mature .and. is_oldgrowth
                     harv_secondary = is_mature .and. (.not. is_oldgrowth)
                     !---------------------------------------------------------------------!
                  end select
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Accumulate harvestable agb.                                        !
                  !------------------------------------------------------------------------!
                  if (harv_primary .or. harv_secondary) then
                     do ico=1,cpatch%ncohorts
                        ipft = cpatch%pft(ico)
                        if (cpatch%dbh(ico) >= cpoly%mindbh_harvest(ipft,isi)) then
                           bharvest  = bharvest + cpoly%prob_harvest(ipft,isi)             &
                                                * cpatch%nplant(ico) * cpatch%btimber(ico)
                        end if  
                     end do
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Accumulate site-level harvest target.                              !
                  !------------------------------------------------------------------------!
                  if (harv_primary) then
                     cpoly%primary_harvest_target(isi) =                                   &
                        cpoly%primary_harvest_target(isi) + bharvest * csite%area(ipa)
                  else if (harv_secondary) then
                     cpoly%secondary_harvest_target(isi) =                                 &
                        cpoly%secondary_harvest_target(isi) + bharvest * csite%area(ipa)
                  end if
                  !------------------------------------------------------------------------!
               end do pharv_loop
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Check whether to scale targets.                                          !
            !------------------------------------------------------------------------------!
            if (ianth_disturb == 2 .and. sl_biomass_harvest > 0.) then
               pot_harvest_target = cpoly%primary_harvest_target  (isi)                    &
                                  + cpoly%secondary_harvest_target(isi)
               if (pot_harvest_target > 0) then
                  cpoly%primary_harvest_target  (isi) =                                    &
                        cpoly%primary_harvest_target  (isi) * sl_biomass_harvest           &
                                                            / pot_harvest_target
                  cpoly%secondary_harvest_target(isi) =                                    &
                        cpoly%secondary_harvest_target(isi) * sl_biomass_harvest           &
                                                            / pot_harvest_target
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     In case anthropogenic disturbance is turned on, make sure that croplands !
            ! that remain croplands are completely ploughed.  This may need to change for  !
            ! a perennial crop.                                                            !
            !------------------------------------------------------------------------------!
            select case (ianth_disturb)
            case (0)
               !----- Anthropogenic disturbance is off. No logging should occur. ----------!
               continue
               !---------------------------------------------------------------------------!
            case default
               !---------------------------------------------------------------------------!
               !   Anthropogenic disturbance is on.  All area that is not abandoned or     !
               ! converted to pasture should be cleared again for cropland.                !
               !---------------------------------------------------------------------------!
               cpoly%disturbance_rates(8,8,isi) = lnexp_max                                &
                                                - sum(cpoly%disturbance_rates(1:7,8,isi))
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!


         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine site_disturbance_rates
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine assigns initial conditions to a patch that has been disturbed     !
   !---------------------------------------------------------------------------------------!
   subroutine initialize_disturbed_patch(csite,np,lsl)

      use ed_state_vars, only : sitetype           & ! structure
                              , patchtype          ! ! structure
      use ed_max_dims  , only : n_pft              ! ! intent(in)
      use ed_type_init , only : init_ed_patch_vars ! ! sub-routine
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype), target      :: csite
      integer       , intent(in)  :: np
      integer       , intent(in)  :: lsl
      !------------------------------------------------------------------------------------!


      !----- Start with an empty patch.  Surviving cohorts will be added later. -----------!
      csite%patch(np)%ncohorts = 0
      !------------------------------------------------------------------------------------!

      !----- Initialise all fast and long-term variables. ---------------------------------!
      call init_ed_patch_vars(csite,np,np,lsl)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Initialise a few other variables, that should be set to zero or to values      !
      ! other than those assigned in init_ed_patch_vars that has been just called.         !
      !------------------------------------------------------------------------------------!
      csite%age                        (np) = 0.0
      csite%fbeam                      (np) = 1.0
      csite%light_type                 (np) = 1
      csite%ustar                      (np) = 0.0
      csite%htry                       (np) = 0.0
      csite%hprev                      (np) = 0.0
      csite%fast_grnd_C                (np) = 0.0
      csite%fast_soil_C                (np) = 0.0
      csite%structural_grnd_C          (np) = 0.0
      csite%structural_soil_C          (np) = 0.0
      csite%structural_grnd_L          (np) = 0.0
      csite%structural_soil_L          (np) = 0.0
      csite%microbial_soil_C           (np) = 0.0
      csite%passive_soil_C             (np) = 0.0
      csite%slow_soil_C                (np) = 0.0
      csite%fast_grnd_N                (np) = 0.0
      csite%fast_soil_N                (np) = 0.0
      csite%structural_grnd_N          (np) = 0.0
      csite%structural_soil_N          (np) = 0.0
      csite%mineralized_soil_N         (np) = 0.0
      csite%sum_dgd                    (np) = 0.0
      csite%sum_chd                    (np) = 0.0
      !------------------------------------------------------------------------------------!

      return
   end subroutine initialize_disturbed_patch
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will include the contribution of each contributing patch (cp) to  !
   ! the new, disturbed patch (np).                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine increment_patch_vars(csite,np,cp,area_fac,cb_enthalpy,can_exner,cb_mass      &
                                  ,cb_molar,cb_water,cb_co2)
      use ed_state_vars, only : sitetype     & ! structure
                              , patchtype    ! ! structure
      use ed_max_dims  , only : n_pft        ! ! intent(in)
      use grid_coms    , only : nzg          ! ! intent(in)
      use ed_misc_coms , only : writing_long & ! intent(in)
                              , writing_eorq & ! intent(in)
                              , writing_dcyc ! ! intent(in)
      use therm_lib    , only : tq2enthalpy  & ! function
                              , press2exner  ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype), target        :: csite
      integer       , intent(in)    :: np
      integer       , intent(in)    :: cp
      real          , intent(in)    :: area_fac
      real          , intent(inout) :: cb_enthalpy
      real          , intent(inout) :: can_exner
      real          , intent(inout) :: cb_mass
      real          , intent(inout) :: cb_molar
      real          , intent(inout) :: cb_water
      real          , intent(inout) :: cb_co2
      !----- Local variables. -------------------------------------------------------------!
      integer                       :: k
      real                          :: can_enthalpy
      real                          :: can_mass
      real                          :: can_molar
      real                          :: can_water
      real                          :: can_co2
      !------------------------------------------------------------------------------------!

      csite%htry                       (np) = csite%htry                       (np)        &
                                            + csite%htry                       (cp)        &
                                            * area_fac
      csite%hprev                      (np) = csite%hprev                      (np)        &
                                            + csite%hprev                      (cp)        &
                                            * area_fac
      csite%fast_grnd_C                (np) = csite%fast_grnd_C                (np)        &
                                            + csite%fast_grnd_C                (cp)        &
                                            * area_fac
      csite%fast_soil_C                (np) = csite%fast_soil_C                (np)        &
                                            + csite%fast_soil_C                (cp)        &
                                            * area_fac
      csite%structural_grnd_C          (np) = csite%structural_grnd_C          (np)        &
                                            + csite%structural_grnd_C          (cp)        &
                                            * area_fac
      csite%structural_soil_C          (np) = csite%structural_soil_C          (np)        &
                                            + csite%structural_soil_C          (cp)        &
                                            * area_fac
      csite%structural_grnd_L          (np) = csite%structural_grnd_L          (np)        &
                                            + csite%structural_grnd_L          (cp)        &
                                            * area_fac
      csite%structural_soil_L          (np) = csite%structural_soil_L          (np)        &
                                            + csite%structural_soil_L          (cp)        &
                                            * area_fac
      csite%microbial_soil_C           (np) = csite%microbial_soil_C           (np)        &
                                            + csite%microbial_soil_C           (cp)        &
                                            * area_fac
      csite%slow_soil_C                (np) = csite%slow_soil_C                (np)        &
                                            + csite%slow_soil_C                (cp)        &
                                            * area_fac
      csite%passive_soil_C             (np) = csite%passive_soil_C             (np)        &
                                            + csite%passive_soil_C             (cp)        &
                                            * area_fac
      csite%fast_grnd_N                (np) = csite%fast_grnd_N                (np)        &
                                            + csite%fast_grnd_N                (cp)        &
                                            * area_fac
      csite%fast_soil_N                (np) = csite%fast_soil_N                (np)        &
                                            + csite%fast_soil_N                (cp)        &
                                            * area_fac
      csite%structural_grnd_N          (np) = csite%structural_grnd_N          (np)        &
                                            + csite%structural_grnd_N          (cp)        &
                                            * area_fac
      csite%structural_soil_N          (np) = csite%structural_soil_N          (np)        &
                                            + csite%structural_soil_N          (cp)        &
                                            * area_fac
      csite%mineralized_soil_N         (np) = csite%mineralized_soil_N         (np)        &
                                            + csite%mineralized_soil_N         (cp)        &
                                            * area_fac
      csite%sum_dgd                    (np) = csite%sum_dgd                    (np)        &
                                            + csite%sum_dgd                    (cp)        &
                                            * area_fac
      csite%sum_chd                    (np) = csite%sum_chd                    (np)        &
                                            + csite%sum_chd                    (cp)        &
                                            * area_fac
      csite%can_temp_pv                (np) = csite%can_temp_pv                (np)        &
                                            + csite%can_temp_pv                (cp)        &
                                            * area_fac
      csite%can_depth                  (np) = csite%can_depth                  (np)        &
                                            + csite%can_depth                  (cp)        &
                                            * area_fac
      csite%ustar                      (np) = csite%ustar                      (np)        &
                                            + csite%ustar                      (cp)        &
                                            * area_fac
      csite%tstar                      (np) = csite%tstar                      (np)        &
                                            + csite%tstar                      (cp)        &
                                            * area_fac
      csite%qstar                      (np) = csite%qstar                      (np)        &
                                            + csite%qstar                      (cp)        &
                                            * area_fac
      csite%cstar                      (np) = csite%cstar                      (np)        &
                                            + csite%cstar                      (cp)        &
                                            * area_fac
      csite%ggbare                     (np) = csite%ggbare                     (np)        &
                                            + csite%ggbare                     (cp)        &
                                            * area_fac
      csite%ggveg                      (np) = csite%ggveg                      (np)        &
                                            + csite%ggveg                      (cp)        &
                                            * area_fac
      csite%rough                      (np) = csite%rough                      (np)        &
                                            + csite%rough                      (cp)        &
                                            * area_fac
      csite%today_fg_C_loss            (np) = csite%today_fg_C_loss            (np)        &
                                            + csite%today_fg_C_loss            (cp)        &
                                            * area_fac
      csite%today_fs_C_loss            (np) = csite%today_fs_C_loss            (np)        &
                                            + csite%today_fs_C_loss            (cp)        &
                                            * area_fac
      csite%today_fg_N_loss            (np) = csite%today_fg_N_loss            (np)        &
                                            + csite%today_fg_N_loss            (cp)        &
                                            * area_fac
      csite%today_fs_N_loss            (np) = csite%today_fs_N_loss            (np)        &
                                            + csite%today_fs_N_loss            (cp)        &
                                            * area_fac
      csite%today_stg_C_loss           (np) = csite%today_stg_C_loss           (np)        &
                                            + csite%today_stg_C_loss           (cp)        &
                                            * area_fac
      csite%today_sts_C_loss           (np) = csite%today_sts_C_loss           (np)        &
                                            + csite%today_sts_C_loss           (cp)        &
                                            * area_fac
      csite%today_stg_L_loss           (np) = csite%today_stg_L_loss           (np)        &
                                            + csite%today_stg_L_loss           (cp)        &
                                            * area_fac
      csite%today_sts_L_loss           (np) = csite%today_sts_L_loss           (np)        &
                                            + csite%today_sts_L_loss           (cp)        &
                                            * area_fac
      csite%today_stg_N_loss           (np) = csite%today_stg_N_loss           (np)        &
                                            + csite%today_stg_N_loss           (cp)        &
                                            * area_fac
      csite%today_sts_N_loss           (np) = csite%today_sts_N_loss           (np)        &
                                            + csite%today_sts_N_loss           (cp)        &
                                            * area_fac
      csite%today_ms_C_loss            (np) = csite%today_ms_C_loss            (np)        &
                                            + csite%today_ms_C_loss            (cp)        &
                                            * area_fac
      csite%today_ss_C_loss            (np) = csite%today_ss_C_loss            (np)        &
                                            + csite%today_ss_C_loss            (cp)        &
                                            * area_fac
      csite%today_ps_C_loss            (np) = csite%today_ps_C_loss            (np)        &
                                            + csite%today_ps_C_loss            (cp)        &
                                            * area_fac
      csite%today_A_decomp             (np) = csite%today_A_decomp             (np)        &
                                            + csite%today_A_decomp             (cp)        &
                                            * area_fac
      csite%today_B_decomp             (np) = csite%today_B_decomp             (np)        &
                                            + csite%today_B_decomp             (cp)        &
                                            * area_fac
      csite%today_Af_decomp            (np) = csite%today_Af_decomp            (np)        &
                                            + csite%today_Af_decomp            (cp)        &
                                            * area_fac
      csite%today_Bf_decomp            (np) = csite%today_Bf_decomp            (np)        &
                                            + csite%today_Bf_decomp            (cp)        &
                                            * area_fac
      csite%today_rh                   (np) = csite%today_rh                   (np)        &
                                            + csite%today_rh                   (cp)        &
                                            * area_fac
      csite%fgc_in                     (np) = csite%fgc_in                     (np)        &
                                            + csite%fgc_in                     (cp)        &
                                            * area_fac
      csite%fsc_in                     (np) = csite%fsc_in                     (np)        &
                                            + csite%fsc_in                     (cp)        &
                                            * area_fac
      csite%stgc_in                    (np) = csite%stgc_in                    (np)        &
                                            + csite%stgc_in                    (cp)        &
                                            * area_fac
      csite%stsc_in                    (np) = csite%stsc_in                    (np)        &
                                            + csite%stsc_in                    (cp)        &
                                            * area_fac
      csite%stgl_in                    (np) = csite%stgl_in                    (np)        &
                                            + csite%stgl_in                    (cp)        &
                                            * area_fac
      csite%stsl_in                    (np) = csite%stsl_in                    (np)        &
                                            + csite%stsl_in                    (cp)        &
                                            * area_fac
      csite%fgn_in                     (np) = csite%fgn_in                     (np)        &
                                            + csite%fgn_in                     (cp)        &
                                            * area_fac
      csite%fsn_in                     (np) = csite%fsn_in                     (np)        &
                                            + csite%fsn_in                     (cp)        &
                                            * area_fac
      csite%stgn_in                    (np) = csite%stgn_in                    (np)        &
                                            + csite%stgn_in                    (cp)        &
                                            * area_fac
      csite%stsn_in                    (np) = csite%stsn_in                    (np)        &
                                            + csite%stsn_in                    (cp)        &
                                            * area_fac
      csite%total_plant_nitrogen_uptake(np) = csite%total_plant_nitrogen_uptake(np)        &
                                            + csite%total_plant_nitrogen_uptake(cp)        &
                                            * area_fac
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Committed carbon emission pools must be included in the disturbed patch, to    !
      ! ensure that carbon is conserved.                                                   !
      !------------------------------------------------------------------------------------!
      csite%commit_storage_resp        (np) = csite%commit_storage_resp        (np)        &
                                            + csite%commit_storage_resp        (cp)        &
                                            * area_fac
      csite%commit_growth_resp         (np) = csite%commit_growth_resp         (np)        &
                                            + csite%commit_growth_resp         (cp)        &
                                            * area_fac
      !------------------------------------------------------------------------------------!


      !----- Update canopy air properties. ------------------------------------------------!
      can_mass     = csite%can_rhos(cp) * csite%can_depth(cp)
      can_molar    = csite%can_dmol(cp) * csite%can_depth(cp)
      can_enthalpy = can_mass    * tq2enthalpy(csite%can_temp(cp),csite%can_shv(cp),.true.)
      can_water    = can_mass    * csite%can_shv(cp)
      can_co2      = can_molar   * csite%can_co2(cp)
      cb_enthalpy  = cb_enthalpy + can_enthalpy                    * area_fac
      can_exner    = can_exner   + press2exner(csite%can_prss(cp)) * area_fac
      cb_mass      = cb_mass     + can_mass                        * area_fac
      cb_molar     = cb_molar    + can_molar                       * area_fac
      cb_water     = cb_water    + can_water                       * area_fac
      cb_co2       = cb_co2      + can_co2                         * area_fac
      !------------------------------------------------------------------------------------!


      !----- Update all budget variables. -------------------------------------------------!
      csite%co2budget_gpp           (np) = csite%co2budget_gpp           (np)              &
                                         + csite%co2budget_gpp           (cp)              &
                                         * area_fac
      csite%co2budget_rh            (np) = csite%co2budget_rh            (np)              &
                                         + csite%co2budget_rh            (cp)              &
                                         * area_fac
      csite%co2budget_plresp        (np) = csite%co2budget_plresp        (np)              &
                                         + csite%co2budget_plresp        (cp)              &
                                         * area_fac
      csite%co2budget_loss2atm      (np) = csite%co2budget_loss2atm      (np)              &
                                         + csite%co2budget_loss2atm      (cp)              &
                                         * area_fac
      csite%co2budget_denseffect    (np) = csite%co2budget_denseffect    (np)              &
                                         + csite%co2budget_denseffect    (cp)              &
                                         * area_fac
      csite%co2budget_zcaneffect    (np) = csite%co2budget_zcaneffect    (np)              &
                                         + csite%co2budget_zcaneffect    (cp)              &
                                         * area_fac
      csite%co2budget_residual      (np) = csite%co2budget_residual      (np)              &
                                         + csite%co2budget_residual      (cp)              &
                                         * area_fac
      csite%cbudget_loss2atm        (np) = csite%cbudget_loss2atm        (np)              &
                                         + csite%cbudget_loss2atm        (cp)              &
                                         * area_fac
      csite%cbudget_denseffect      (np) = csite%cbudget_denseffect      (np)              &
                                         + csite%cbudget_denseffect      (cp)              &
                                         * area_fac
      csite%cbudget_zcaneffect      (np) = csite%cbudget_zcaneffect      (np)              &
                                         + csite%cbudget_zcaneffect      (cp)              &
                                         * area_fac
      csite%cbudget_seedrain        (np) = csite%cbudget_seedrain        (np)              &
                                         + csite%cbudget_seedrain        (cp)              &
                                         * area_fac
      csite%cbudget_loss2yield      (np) = csite%cbudget_loss2yield      (np)              &
                                         + csite%cbudget_loss2yield      (cp)              &
                                         * area_fac
      csite%cbudget_residual        (np) = csite%cbudget_residual        (np)              &
                                         + csite%cbudget_residual        (cp)              &
                                         * area_fac
      csite%wbudget_precipgain      (np) = csite%wbudget_precipgain      (np)              &
                                         + csite%wbudget_precipgain      (cp)              &
                                         * area_fac
      csite%wbudget_loss2atm        (np) = csite%wbudget_loss2atm        (np)              &
                                         + csite%wbudget_loss2atm        (cp)              &
                                         * area_fac
      csite%wbudget_loss2runoff     (np) = csite%wbudget_loss2runoff     (np)              &
                                         + csite%wbudget_loss2runoff     (cp)              &
                                         * area_fac
      csite%wbudget_loss2drainage   (np) = csite%wbudget_loss2drainage   (np)              &
                                         + csite%wbudget_loss2drainage   (cp)              &
                                         * area_fac
      csite%wbudget_denseffect      (np) = csite%wbudget_denseffect      (np)              &
                                         + csite%wbudget_denseffect      (cp)              &
                                         * area_fac
      csite%wbudget_wcapeffect      (np) = csite%wbudget_wcapeffect      (np)              &
                                         + csite%wbudget_wcapeffect      (cp)              &
                                         * area_fac
      csite%wbudget_zcaneffect      (np) = csite%wbudget_zcaneffect      (np)              &
                                         + csite%wbudget_zcaneffect      (cp)              &
                                         * area_fac
      csite%wbudget_pheneffect      (np) = csite%wbudget_pheneffect      (np)              &
                                         + csite%wbudget_pheneffect      (cp)              &
                                         * area_fac
      csite%wbudget_residual        (np) = csite%wbudget_residual        (np)              &
                                         + csite%wbudget_residual        (cp)              &
                                         * area_fac
      csite%ebudget_precipgain      (np) = csite%ebudget_precipgain      (np)              &
                                         + csite%ebudget_precipgain      (cp)              &
                                         * area_fac
      csite%ebudget_netrad          (np) = csite%ebudget_netrad          (np)              &
                                         + csite%ebudget_netrad          (cp)              &
                                         * area_fac
      csite%ebudget_loss2atm        (np) = csite%ebudget_loss2atm        (np)              &
                                         + csite%ebudget_loss2atm        (cp)              &
                                         * area_fac
      csite%ebudget_loss2runoff     (np) = csite%ebudget_loss2runoff     (np)              &
                                         + csite%ebudget_loss2runoff     (cp)              &
                                         * area_fac
      csite%ebudget_loss2drainage   (np) = csite%ebudget_loss2drainage   (np)              &
                                         + csite%ebudget_loss2drainage   (cp)              &
                                         * area_fac
      csite%ebudget_denseffect      (np) = csite%ebudget_denseffect      (np)              &
                                         + csite%ebudget_denseffect      (cp)              &
                                         * area_fac
      csite%ebudget_prsseffect      (np) = csite%ebudget_prsseffect      (np)              &
                                         + csite%ebudget_prsseffect      (cp)              &
                                         * area_fac
      csite%ebudget_hcapeffect      (np) = csite%ebudget_hcapeffect      (np)              &
                                         + csite%ebudget_hcapeffect      (cp)              &
                                         * area_fac
      csite%ebudget_wcapeffect      (np) = csite%ebudget_wcapeffect      (np)              &
                                         + csite%ebudget_wcapeffect      (cp)              &
                                         * area_fac
      csite%ebudget_zcaneffect      (np) = csite%ebudget_zcaneffect      (np)              &
                                         + csite%ebudget_zcaneffect      (cp)              &
                                         * area_fac
      csite%ebudget_pheneffect      (np) = csite%ebudget_pheneffect      (np)              &
                                         + csite%ebudget_pheneffect      (cp)              &
                                         * area_fac
      csite%ebudget_residual        (np) = csite%ebudget_residual        (np)              &
                                         + csite%ebudget_residual        (cp)              &
                                         * area_fac
      !------------------------------------------------------------------------------------!


      !----- Reproduction array. ----------------------------------------------------------!
      do k=1,n_pft
         csite%repro                 (k,np) = csite%repro                    (k,np)        &
                                            + csite%repro                    (k,cp)        &
                                            * area_fac
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Soil variables.  We merge internal energy and water, then find temperature and !
      ! liquid water fraction once all patches have their contribution.                    !
      !------------------------------------------------------------------------------------!
      do k = 1, nzg
         csite%soil_energy           (k,np) = csite%soil_energy              (k,np)        &
                                            + csite%soil_energy              (k,cp)        &
                                            * area_fac
         csite%soil_water            (k,np) = csite%soil_water               (k,np)        &
                                            + csite%soil_water               (k,cp)        &
                                            * area_fac
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Update temporary surface water layer count in case there is any contributor   !
      ! with more layers.                                                                  !
      !------------------------------------------------------------------------------------!
      if (csite%nlev_sfcwater(cp) > csite%nlev_sfcwater(np)) then
         !----- Increase the number of layers. --------------------------------------------!
         csite%nlev_sfcwater(np) = 1
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Integrate layers.                                                          !
         !---------------------------------------------------------------------------------!
         do k = 1, csite%nlev_sfcwater(cp)
            csite%sfcwater_mass         (1,np) = csite%sfcwater_mass            (1,np)     &
                                               + csite%sfcwater_mass            (k,cp)     &
                                               * area_fac
            csite%sfcwater_energy       (1,np) = csite%sfcwater_energy          (1,np)     &
                                               + csite%sfcwater_energy          (k,cp)     &
                                               * csite%sfcwater_mass            (k,cp)     &
                                               * area_fac
            csite%sfcwater_depth        (1,np) = csite%sfcwater_depth           (1,np)     &
                                               + csite%sfcwater_depth           (k,cp)     &
                                               * area_fac
         end do
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Fast means must be aggregated as well.                                         !
      !------------------------------------------------------------------------------------!
      csite%fmean_rh             (np) = csite%fmean_rh             (np)                    &
                                      + csite%fmean_rh             (cp)                    &
                                      * area_fac
      csite%fmean_fgc_rh         (np) = csite%fmean_fgc_rh         (np)                    &
                                      + csite%fmean_fgc_rh         (cp)                    &
                                      * area_fac
      csite%fmean_fsc_rh         (np) = csite%fmean_fsc_rh         (np)                    &
                                      + csite%fmean_fsc_rh         (cp)                    &
                                      * area_fac
      csite%fmean_stgc_rh        (np) = csite%fmean_stgc_rh        (np)                    &
                                      + csite%fmean_stgc_rh        (cp)                    &
                                      * area_fac
      csite%fmean_stsc_rh        (np) = csite%fmean_stsc_rh        (np)                    &
                                      + csite%fmean_stsc_rh        (cp)                    &
                                      * area_fac
      csite%fmean_msc_rh         (np) = csite%fmean_msc_rh         (np)                    &
                                      + csite%fmean_msc_rh         (cp)                    &
                                      * area_fac
      csite%fmean_ssc_rh         (np) = csite%fmean_ssc_rh         (np)                    &
                                      + csite%fmean_ssc_rh         (cp)                    &
                                      * area_fac
      csite%fmean_psc_rh         (np) = csite%fmean_psc_rh         (np)                    &
                                      + csite%fmean_psc_rh         (cp)                    &
                                      * area_fac
      csite%fmean_nep            (np) = csite%fmean_nep            (np)                    &
                                      + csite%fmean_nep            (cp)                    &
                                      * area_fac
      csite%fmean_rk4step        (np) = csite%fmean_rk4step        (np)                    &
                                      + csite%fmean_rk4step        (cp)                    &
                                      * area_fac
      csite%fmean_available_water(np) = csite%fmean_available_water(np)                    &
                                      + csite%fmean_available_water(cp)                    &
                                      * area_fac 
      csite%fmean_veg_displace   (np) = csite%fmean_veg_displace   (np)                    &
                                      + csite%fmean_veg_displace   (cp)                    &
                                      * area_fac
      csite%fmean_rough          (np) = csite%fmean_rough          (np)                    &
                                      + csite%fmean_rough          (cp)                    &
                                      * area_fac
      csite%fmean_can_theiv      (np) = csite%fmean_can_theiv      (np)                    &
                                      + csite%fmean_can_theiv      (cp)                    &
                                      * area_fac
      csite%fmean_can_theta      (np) = csite%fmean_can_theta      (np)                    &
                                      + csite%fmean_can_theta      (cp)                    &
                                      * area_fac
      csite%fmean_can_vpdef      (np) = csite%fmean_can_vpdef      (np)                    &
                                      + csite%fmean_can_vpdef      (cp)                    &
                                      * area_fac
      csite%fmean_can_shv        (np) = csite%fmean_can_shv        (np)                    &
                                      + csite%fmean_can_shv        (cp)                    &
                                      * area_fac
      csite%fmean_can_co2        (np) = csite%fmean_can_co2        (np)                    &
                                      + csite%fmean_can_co2        (cp)                    &
                                      * area_fac
      csite%fmean_can_prss       (np) = csite%fmean_can_prss       (np)                    &
                                      + csite%fmean_can_prss       (cp)                    &
                                      * area_fac
      csite%fmean_gnd_temp       (np) = csite%fmean_gnd_temp       (np)                    &
                                      + csite%fmean_gnd_temp       (cp)                    &
                                      * area_fac
      csite%fmean_gnd_shv        (np) = csite%fmean_gnd_shv        (np)                    &
                                      + csite%fmean_gnd_shv        (cp)                    &
                                      * area_fac
      csite%fmean_can_ggnd       (np) = csite%fmean_can_ggnd       (np)                    &
                                      + csite%fmean_can_ggnd       (cp)                    &
                                      * area_fac
      csite%fmean_sfcw_depth     (np) = csite%fmean_sfcw_depth     (np)                    &
                                      + csite%fmean_sfcw_depth     (cp)                    &
                                      * area_fac
      !----- Integrate pounding energy in J/m2. -------------------------------------------!
      csite%fmean_sfcw_energy    (np) = csite%fmean_sfcw_energy    (np)                    &
                                      + csite%fmean_sfcw_energy    (cp)                    &
                                      * csite%fmean_sfcw_mass      (cp)                    &
                                      * area_fac
      csite%fmean_sfcw_mass      (np) = csite%fmean_sfcw_mass      (np)                    &
                                      + csite%fmean_sfcw_mass      (cp)                    &
                                      * area_fac
      csite%fmean_rshort_gnd     (np) = csite%fmean_rshort_gnd     (np)                    &
                                      + csite%fmean_rshort_gnd     (cp)                    &
                                      * area_fac
      csite%fmean_par_gnd        (np) = csite%fmean_par_gnd        (np)                    &
                                      + csite%fmean_par_gnd        (cp)                    &
                                      * area_fac
      csite%fmean_rlong_gnd      (np) = csite%fmean_rlong_gnd      (np)                    &
                                      + csite%fmean_rlong_gnd      (cp)                    &
                                      * area_fac
      csite%fmean_rlongup        (np) = csite%fmean_rlongup        (np)                    &
                                      + csite%fmean_rlongup        (cp)                    &
                                      * area_fac
      csite%fmean_parup          (np) = csite%fmean_parup          (np)                    &
                                      + csite%fmean_parup          (cp)                    &
                                      * area_fac
      csite%fmean_nirup          (np) = csite%fmean_nirup          (np)                    &
                                      + csite%fmean_nirup          (cp)                    &
                                      * area_fac
      csite%fmean_rshortup       (np) = csite%fmean_rshortup       (np)                    &
                                      + csite%fmean_rshortup       (cp)                    &
                                      * area_fac
      csite%fmean_rnet           (np) = csite%fmean_rnet           (np)                    &
                                      + csite%fmean_rnet           (cp)                    &
                                      * area_fac
      csite%fmean_albedo         (np) = csite%fmean_albedo         (np)                    &
                                      + csite%fmean_albedo         (cp)                    &
                                      * area_fac
      csite%fmean_albedo_par     (np) = csite%fmean_albedo_par     (np)                    &
                                      + csite%fmean_albedo_par     (cp)                    &
                                      * area_fac
      csite%fmean_albedo_nir     (np) = csite%fmean_albedo_nir     (np)                    &
                                      + csite%fmean_albedo_nir     (cp)                    &
                                      * area_fac
      csite%fmean_rlong_albedo   (np) = csite%fmean_rlong_albedo   (np)                    &
                                      + csite%fmean_rlong_albedo   (cp)                    &
                                      * area_fac
      csite%fmean_ustar          (np) = csite%fmean_ustar          (np)                    &
                                      + csite%fmean_ustar          (cp)                    &
                                      * area_fac
      csite%fmean_tstar          (np) = csite%fmean_tstar          (np)                    &
                                      + csite%fmean_tstar          (cp)                    &
                                      * area_fac
      csite%fmean_qstar          (np) = csite%fmean_qstar          (np)                    &
                                      + csite%fmean_qstar          (cp)                    &
                                      * area_fac
      csite%fmean_cstar          (np) = csite%fmean_cstar          (np)                    &
                                      + csite%fmean_cstar          (cp)                    &
                                      * area_fac
      csite%fmean_carbon_ac      (np) = csite%fmean_carbon_ac      (np)                    &
                                      + csite%fmean_carbon_ac      (cp)                    &
                                      * area_fac
      csite%fmean_carbon_st      (np) = csite%fmean_carbon_st      (np)                    &
                                      + csite%fmean_carbon_st      (cp)                    &
                                      * area_fac
      csite%fmean_vapor_gc       (np) = csite%fmean_vapor_gc       (np)                    &
                                      + csite%fmean_vapor_gc       (cp)                    &
                                      * area_fac
      csite%fmean_vapor_ac       (np) = csite%fmean_vapor_ac       (np)                    &
                                      + csite%fmean_vapor_ac       (cp)                    &
                                      * area_fac
      csite%fmean_throughfall    (np) = csite%fmean_throughfall    (np)                    &
                                      + csite%fmean_throughfall    (cp)                    &
                                      * area_fac
      csite%fmean_runoff         (np) = csite%fmean_runoff         (np)                    &
                                      + csite%fmean_runoff         (cp)                    &
                                      * area_fac
      csite%fmean_drainage       (np) = csite%fmean_drainage       (np)                    &
                                      + csite%fmean_drainage       (cp)                    &
                                      * area_fac
      csite%fmean_sensible_gc    (np) = csite%fmean_sensible_gc    (np)                    &
                                      + csite%fmean_sensible_gc    (cp)                    &
                                      * area_fac
      csite%fmean_sensible_ac    (np) = csite%fmean_sensible_ac    (np)                    &
                                      + csite%fmean_sensible_ac    (cp)                    &
                                      * area_fac
      csite%fmean_qthroughfall   (np) = csite%fmean_qthroughfall   (np)                    &
                                      + csite%fmean_qthroughfall   (cp)                    &
                                      * area_fac
      csite%fmean_qrunoff        (np) = csite%fmean_qrunoff        (np)                    &
                                      + csite%fmean_qrunoff        (cp)                    &
                                      * area_fac
      csite%fmean_qdrainage      (np) = csite%fmean_qdrainage      (np)                    &
                                      + csite%fmean_qdrainage      (cp)                    &
                                      * area_fac

      do k=1, nzg
         csite%fmean_soil_energy(k,np) = csite%fmean_soil_energy(k,np)                     &
                                       + csite%fmean_soil_energy(k,cp)                     &
                                       * area_fac
         csite%fmean_soil_water (k,np) = csite%fmean_soil_water (k,np)                     &
                                       + csite%fmean_soil_water (k,cp)                     &
                                       * area_fac
         csite%fmean_smoist_gg  (k,np) = csite%fmean_smoist_gg  (k,np)                     &
                                       + csite%fmean_smoist_gg  (k,cp)                     &
                                       * area_fac
         csite%fmean_transloss  (k,np) = csite%fmean_transloss  (k,np)                     &
                                       + csite%fmean_transloss  (k,cp)                     &
                                       * area_fac
         csite%fmean_sensible_gg(k,np) = csite%fmean_sensible_gg(k,np)                     &
                                       + csite%fmean_sensible_gg(k,cp)                     &
                                       * area_fac
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Daily means...                                                                 !
      !------------------------------------------------------------------------------------!
      if (writing_long) then
         csite%dmean_A_decomp       (    np) = csite%dmean_A_decomp       (    np)         &
                                             + csite%dmean_A_decomp       (    cp)         &
                                             * area_fac
         csite%dmean_B_decomp       (    np) = csite%dmean_B_decomp       (    np)         &
                                             + csite%dmean_B_decomp       (    cp)         &
                                             * area_fac
         csite%dmean_Af_decomp      (    np) = csite%dmean_Af_decomp      (    np)         &
                                             + csite%dmean_Af_decomp      (    cp)         &
                                             * area_fac
         csite%dmean_Bf_decomp      (    np) = csite%dmean_Bf_decomp      (    np)         &
                                             + csite%dmean_Bf_decomp      (    cp)         &
                                             * area_fac
         csite%dmean_co2_residual   (    np) = csite%dmean_co2_residual   (    np)         &
                                             + csite%dmean_co2_residual   (    cp)         &
                                             * area_fac
         csite%dmean_energy_residual(    np) = csite%dmean_energy_residual(    np)         &
                                             + csite%dmean_energy_residual(    cp)         &
                                             * area_fac
         csite%dmean_water_residual (    np) = csite%dmean_water_residual (    np)         &
                                             + csite%dmean_water_residual (    cp)         &
                                             * area_fac
         csite%dmean_rh             (    np) = csite%dmean_rh             (    np)         &
                                             + csite%dmean_rh             (    cp)         &
                                             * area_fac
         csite%dmean_fgc_rh         (    np) = csite%dmean_fgc_rh         (    np)         &
                                             + csite%dmean_fgc_rh         (    cp)         &
                                             * area_fac
         csite%dmean_fsc_rh         (    np) = csite%dmean_fsc_rh         (    np)         &
                                             + csite%dmean_fsc_rh         (    cp)         &
                                             * area_fac
         csite%dmean_stgc_rh        (    np) = csite%dmean_stgc_rh        (    np)         &
                                             + csite%dmean_stgc_rh        (    cp)         &
                                             * area_fac
         csite%dmean_stsc_rh        (    np) = csite%dmean_stsc_rh        (    np)         &
                                             + csite%dmean_stsc_rh        (    cp)         &
                                             * area_fac
         csite%dmean_msc_rh         (    np) = csite%dmean_msc_rh         (    np)         &
                                             + csite%dmean_msc_rh         (    cp)         &
                                             * area_fac
         csite%dmean_ssc_rh         (    np) = csite%dmean_ssc_rh         (    np)         &
                                             + csite%dmean_ssc_rh         (    cp)         &
                                             * area_fac
         csite%dmean_psc_rh         (    np) = csite%dmean_psc_rh         (    np)         &
                                             + csite%dmean_psc_rh         (    cp)         &
                                             * area_fac
         csite%dmean_nep            (    np) = csite%dmean_nep            (    np)         &
                                             + csite%dmean_nep            (    cp)         &
                                             * area_fac
         csite%dmean_rk4step        (    np) = csite%dmean_rk4step        (    np)         &
                                             + csite%dmean_rk4step        (    cp)         &
                                             * area_fac
         csite%dmean_available_water(    np) = csite%dmean_available_water(    np)         &
                                             + csite%dmean_available_water(    cp)         &
                                             * area_fac
         csite%dmean_veg_displace   (    np) = csite%dmean_veg_displace   (    np)         &
                                             + csite%dmean_veg_displace   (    cp)         &
                                             * area_fac
         csite%dmean_rough          (    np) = csite%dmean_rough          (    np)         &
                                             + csite%dmean_rough          (    cp)         &
                                             * area_fac
         csite%dmean_can_theiv      (    np) = csite%dmean_can_theiv      (    np)         &
                                             + csite%dmean_can_theiv      (    cp)         &
                                             * area_fac
         csite%dmean_can_theta      (    np) = csite%dmean_can_theta      (    np)         &
                                             + csite%dmean_can_theta      (    cp)         &
                                             * area_fac
         csite%dmean_can_vpdef      (    np) = csite%dmean_can_vpdef      (    np)         &
                                             + csite%dmean_can_vpdef      (    cp)         &
                                             * area_fac
         csite%dmean_can_temp       (    np) = csite%dmean_can_temp       (    np)         &
                                             + csite%dmean_can_temp       (    cp)         &
                                             * area_fac
         csite%dmean_can_shv        (    np) = csite%dmean_can_shv        (    np)         &
                                             + csite%dmean_can_shv        (    cp)         &
                                             * area_fac
         csite%dmean_can_co2        (    np) = csite%dmean_can_co2        (    np)         &
                                             + csite%dmean_can_co2        (    cp)         &
                                             * area_fac
         csite%dmean_can_rhos       (    np) = csite%dmean_can_rhos       (    np)         &
                                             + csite%dmean_can_rhos       (    cp)         &
                                             * area_fac
         csite%dmean_can_dmol       (    np) = csite%dmean_can_dmol       (    np)         &
                                             + csite%dmean_can_dmol       (    cp)         &
                                             * area_fac
         csite%dmean_can_prss       (    np) = csite%dmean_can_prss       (    np)         &
                                             + csite%dmean_can_prss       (    cp)         &
                                             * area_fac
         csite%dmean_gnd_temp       (    np) = csite%dmean_gnd_temp       (    np)         &
                                             + csite%dmean_gnd_temp       (    cp)         &
                                             * area_fac
         csite%dmean_gnd_shv        (    np) = csite%dmean_gnd_shv        (    np)         &
                                             + csite%dmean_gnd_shv        (    cp)         &
                                             * area_fac
         csite%dmean_can_ggnd       (    np) = csite%dmean_can_ggnd       (    np)         &
                                             + csite%dmean_can_ggnd       (    cp)         &
                                             * area_fac
         csite%dmean_sfcw_depth     (    np) = csite%dmean_sfcw_depth     (    np)         &
                                             + csite%dmean_sfcw_depth     (    cp)         &
                                             * area_fac
         csite%dmean_sfcw_energy    (    np) = csite%dmean_sfcw_energy    (    np)         &
                                             + csite%dmean_sfcw_energy    (    cp)         &
                                             * area_fac
         csite%dmean_sfcw_mass      (    np) = csite%dmean_sfcw_mass      (    np)         &
                                             + csite%dmean_sfcw_mass      (    cp)         &
                                             * area_fac
         csite%dmean_sfcw_temp      (    np) = csite%dmean_sfcw_temp      (    np)         &
                                             + csite%dmean_sfcw_temp      (    cp)         &
                                             * area_fac
         csite%dmean_sfcw_fliq      (    np) = csite%dmean_sfcw_fliq      (    np)         &
                                             + csite%dmean_sfcw_fliq      (    cp)         &
                                             * area_fac
         csite%dmean_rshort_gnd     (    np) = csite%dmean_rshort_gnd     (    np)         &
                                             + csite%dmean_rshort_gnd     (    cp)         &
                                             * area_fac
         csite%dmean_par_gnd        (    np) = csite%dmean_par_gnd        (    np)         &
                                             + csite%dmean_par_gnd        (    cp)         &
                                             * area_fac
         csite%dmean_rlong_gnd      (    np) = csite%dmean_rlong_gnd      (    np)         &
                                             + csite%dmean_rlong_gnd      (    cp)         &
                                             * area_fac
         csite%dmean_rlongup        (    np) = csite%dmean_rlongup        (    np)         &
                                             + csite%dmean_rlongup        (    cp)         &
                                             * area_fac
         csite%dmean_parup          (    np) = csite%dmean_parup          (    np)         &
                                             + csite%dmean_parup          (    cp)         &
                                             * area_fac
         csite%dmean_nirup          (    np) = csite%dmean_nirup          (    np)         &
                                             + csite%dmean_nirup          (    cp)         &
                                             * area_fac
         csite%dmean_rshortup       (    np) = csite%dmean_rshortup       (    np)         &
                                             + csite%dmean_rshortup       (    cp)         &
                                             * area_fac
         csite%dmean_rnet           (    np) = csite%dmean_rnet           (    np)         &
                                             + csite%dmean_rnet           (    cp)         &
                                             * area_fac
         csite%dmean_albedo         (    np) = csite%dmean_albedo         (    np)         &
                                             + csite%dmean_albedo         (    cp)         &
                                             * area_fac
         csite%dmean_albedo_par     (    np) = csite%dmean_albedo_par     (    np)         &
                                             + csite%dmean_albedo_par     (    cp)         &
                                             * area_fac
         csite%dmean_albedo_nir     (    np) = csite%dmean_albedo_nir     (    np)         &
                                             + csite%dmean_albedo_nir     (    cp)         &
                                             * area_fac
         csite%dmean_rlong_albedo   (    np) = csite%dmean_rlong_albedo   (    np)         &
                                             + csite%dmean_rlong_albedo   (    cp)         &
                                             * area_fac
         csite%dmean_ustar          (    np) = csite%dmean_ustar          (    np)         &
                                             + csite%dmean_ustar          (    cp)         &
                                             * area_fac
         csite%dmean_tstar          (    np) = csite%dmean_tstar          (    np)         &
                                             + csite%dmean_tstar          (    cp)         &
                                             * area_fac
         csite%dmean_qstar          (    np) = csite%dmean_qstar          (    np)         &
                                             + csite%dmean_qstar          (    cp)         &
                                             * area_fac
         csite%dmean_cstar          (    np) = csite%dmean_cstar          (    np)         &
                                             + csite%dmean_cstar          (    cp)         &
                                             * area_fac
         csite%dmean_carbon_ac      (    np) = csite%dmean_carbon_ac      (    np)         &
                                             + csite%dmean_carbon_ac      (    cp)         &
                                             * area_fac
         csite%dmean_carbon_st      (    np) = csite%dmean_carbon_st      (    np)         &
                                             + csite%dmean_carbon_st      (    cp)         &
                                             * area_fac
         csite%dmean_vapor_gc       (    np) = csite%dmean_vapor_gc       (    np)         &
                                             + csite%dmean_vapor_gc       (    cp)         &
                                             * area_fac
         csite%dmean_vapor_ac       (    np) = csite%dmean_vapor_ac       (    np)         &
                                             + csite%dmean_vapor_ac       (    cp)         &
                                             * area_fac
         csite%dmean_throughfall    (    np) = csite%dmean_throughfall    (    np)         &
                                             + csite%dmean_throughfall    (    cp)         &
                                             * area_fac
         csite%dmean_runoff         (    np) = csite%dmean_runoff         (    np)         &
                                             + csite%dmean_runoff         (    cp)         &
                                             * area_fac
         csite%dmean_drainage       (    np) = csite%dmean_drainage       (    np)         &
                                             + csite%dmean_drainage       (    cp)         &
                                             * area_fac
         csite%dmean_sensible_gc    (    np) = csite%dmean_sensible_gc    (    np)         &
                                             + csite%dmean_sensible_gc    (    cp)         &
                                             * area_fac
         csite%dmean_sensible_ac    (    np) = csite%dmean_sensible_ac    (    np)         &
                                             + csite%dmean_sensible_ac    (    cp)         &
                                             * area_fac
         csite%dmean_qthroughfall   (    np) = csite%dmean_qthroughfall   (    np)         &
                                             + csite%dmean_qthroughfall   (    cp)         &
                                             * area_fac
         csite%dmean_qrunoff        (    np) = csite%dmean_qrunoff        (    np)         &
                                             + csite%dmean_qrunoff        (    cp)         &
                                             * area_fac
         csite%dmean_qdrainage      (    np) = csite%dmean_qdrainage      (    np)         &
                                             + csite%dmean_qdrainage      (    cp)         &
                                             * area_fac
         csite%dmean_soil_energy    (  :,np) = csite%dmean_soil_energy    (  :,np)         &
                                             + csite%dmean_soil_energy    (  :,cp)         &
                                             * area_fac
         csite%dmean_soil_mstpot    (  :,np) = csite%dmean_soil_mstpot    (  :,np)         &
                                             + csite%dmean_soil_mstpot    (  :,cp)         &
                                             * area_fac
         csite%dmean_soil_water     (  :,np) = csite%dmean_soil_water     (  :,np)         &
                                             + csite%dmean_soil_water     (  :,cp)         &
                                             * area_fac
         csite%dmean_soil_temp      (  :,np) = csite%dmean_soil_temp      (  :,np)         &
                                             + csite%dmean_soil_temp      (  :,cp)         &
                                             * area_fac
         csite%dmean_soil_fliq      (  :,np) = csite%dmean_soil_fliq      (  :,np)         &
                                             + csite%dmean_soil_fliq      (  :,cp)         &
                                             * area_fac
         csite%dmean_smoist_gg      (  :,np) = csite%dmean_smoist_gg      (  :,np)         &
                                             + csite%dmean_smoist_gg      (  :,cp)         &
                                             * area_fac
         csite%dmean_transloss      (  :,np) = csite%dmean_transloss      (  :,np)         &
                                             + csite%dmean_transloss      (  :,cp)         &
                                             * area_fac
         csite%dmean_sensible_gg    (  :,np) = csite%dmean_sensible_gg    (  :,np)         &
                                             + csite%dmean_sensible_gg    (  :,cp)         &
                                             * area_fac
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Monthly means...                                                               !
      !------------------------------------------------------------------------------------!
      if (writing_eorq) then
         csite%mmean_fast_grnd_c    (    np) = csite%mmean_fast_grnd_c    (    np)         &
                                             + csite%mmean_fast_grnd_c    (    cp)         &
                                             * area_fac
         csite%mmean_fast_soil_c    (    np) = csite%mmean_fast_soil_c    (    np)         &
                                             + csite%mmean_fast_soil_c    (    cp)         &
                                             * area_fac
         csite%mmean_struct_grnd_c  (    np) = csite%mmean_struct_grnd_c  (    np)         &
                                             + csite%mmean_struct_grnd_c  (    cp)         &
                                             * area_fac
         csite%mmean_struct_soil_c  (    np) = csite%mmean_struct_soil_c  (    np)         &
                                             + csite%mmean_struct_soil_c  (    cp)         &
                                             * area_fac
         csite%mmean_struct_grnd_l  (    np) = csite%mmean_struct_grnd_l  (    np)         &
                                             + csite%mmean_struct_grnd_l  (    cp)         &
                                             * area_fac
         csite%mmean_struct_soil_l  (    np) = csite%mmean_struct_soil_l  (    np)         &
                                             + csite%mmean_struct_soil_l  (    cp)         &
                                             * area_fac
         csite%mmean_microbe_soil_c (    np) = csite%mmean_microbe_soil_c (    np)         &
                                             + csite%mmean_microbe_soil_c (    cp)         &
                                             * area_fac
         csite%mmean_slow_soil_c    (    np) = csite%mmean_slow_soil_c    (    np)         &
                                             + csite%mmean_slow_soil_c    (    cp)         &
                                             * area_fac
         csite%mmean_passive_soil_c (    np) = csite%mmean_passive_soil_c (    np)         &
                                             + csite%mmean_passive_soil_c (    cp)         &
                                             * area_fac
         csite%mmean_fast_grnd_n    (    np) = csite%mmean_fast_grnd_n    (    np)         &
                                             + csite%mmean_fast_grnd_n    (    cp)         &
                                             * area_fac
         csite%mmean_fast_soil_n    (    np) = csite%mmean_fast_soil_n    (    np)         &
                                             + csite%mmean_fast_soil_n    (    cp)         &
                                             * area_fac
         csite%mmean_struct_grnd_n  (    np) = csite%mmean_struct_grnd_n  (    np)         &
                                             + csite%mmean_struct_grnd_n  (    cp)         &
                                             * area_fac
         csite%mmean_struct_soil_n  (    np) = csite%mmean_struct_soil_n  (    np)         &
                                             + csite%mmean_struct_soil_n  (    cp)         &
                                             * area_fac
         csite%mmean_mineral_soil_n (    np) = csite%mmean_mineral_soil_n (    np)         &
                                             + csite%mmean_mineral_soil_n (    cp)         &
                                             * area_fac
         csite%mmean_fgc_in         (    np) = csite%mmean_fgc_in         (    np)         &
                                             + csite%mmean_fgc_in         (    cp)         &
                                             * area_fac
         csite%mmean_fsc_in         (    np) = csite%mmean_fsc_in         (    np)         &
                                             + csite%mmean_fsc_in         (    cp)         &
                                             * area_fac
         csite%mmean_stgc_in        (    np) = csite%mmean_stgc_in         (    np)         &
                                             + csite%mmean_stgc_in         (    cp)         &
                                             * area_fac
         csite%mmean_stsc_in        (    np) = csite%mmean_stsc_in         (    np)         &
                                             + csite%mmean_stsc_in         (    cp)         &
                                             * area_fac
         csite%mmean_co2_residual   (    np) = csite%mmean_co2_residual   (    np)         &
                                             + csite%mmean_co2_residual   (    cp)         &
                                             * area_fac
         csite%mmean_energy_residual(    np) = csite%mmean_energy_residual(    np)         &
                                             + csite%mmean_energy_residual(    cp)         &
                                             * area_fac
         csite%mmean_water_residual (    np) = csite%mmean_water_residual (    np)         &
                                             + csite%mmean_water_residual (    cp)         &
                                             * area_fac
         csite%mmean_rh             (    np) = csite%mmean_rh             (    np)         &
                                             + csite%mmean_rh             (    cp)         &
                                             * area_fac
         csite%mmean_fgc_rh         (    np) = csite%mmean_fgc_rh         (    np)         &
                                             + csite%mmean_fgc_rh         (    cp)         &
                                             * area_fac
         csite%mmean_fsc_rh         (    np) = csite%mmean_fsc_rh         (    np)         &
                                             + csite%mmean_fsc_rh         (    cp)         &
                                             * area_fac
         csite%mmean_stgc_rh        (    np) = csite%mmean_stgc_rh        (    np)         &
                                             + csite%mmean_stgc_rh        (    cp)         &
                                             * area_fac
         csite%mmean_stsc_rh        (    np) = csite%mmean_stsc_rh        (    np)         &
                                             + csite%mmean_stsc_rh        (    cp)         &
                                             * area_fac
         csite%mmean_msc_rh         (    np) = csite%mmean_msc_rh         (    np)         &
                                             + csite%mmean_msc_rh         (    cp)         &
                                             * area_fac
         csite%mmean_ssc_rh         (    np) = csite%mmean_ssc_rh         (    np)         &
                                             + csite%mmean_ssc_rh         (    cp)         &
                                             * area_fac
         csite%mmean_psc_rh         (    np) = csite%mmean_psc_rh         (    np)         &
                                             + csite%mmean_psc_rh         (    cp)         &
                                             * area_fac
         csite%mmean_nep            (    np) = csite%mmean_nep            (    np)         &
                                             + csite%mmean_nep            (    cp)         &
                                             * area_fac
         csite%mmean_A_decomp       (    np) = csite%mmean_A_decomp       (    np)         &
                                             + csite%mmean_A_decomp       (    cp)         &
                                             * area_fac
         csite%mmean_B_decomp       (    np) = csite%mmean_B_decomp       (    np)         &
                                             + csite%mmean_B_decomp       (    cp)         &
                                             * area_fac
         csite%mmean_Af_decomp      (    np) = csite%mmean_Af_decomp      (    np)         &
                                             + csite%mmean_Af_decomp      (    cp)         &
                                             * area_fac
         csite%mmean_Bf_decomp      (    np) = csite%mmean_Bf_decomp      (    np)         &
                                             + csite%mmean_Bf_decomp      (    cp)         &
                                             * area_fac
         csite%mmean_rk4step        (    np) = csite%mmean_rk4step        (    np)         &
                                             + csite%mmean_rk4step        (    cp)         &
                                             * area_fac
         csite%mmean_available_water(    np) = csite%mmean_available_water(    np)         &
                                             + csite%mmean_available_water(    cp)         &
                                             * area_fac
         csite%mmean_veg_displace   (    np) = csite%mmean_veg_displace   (    np)         &
                                             + csite%mmean_veg_displace   (    cp)         &
                                             * area_fac
         csite%mmean_rough          (    np) = csite%mmean_rough          (    np)         &
                                             + csite%mmean_rough          (    cp)         &
                                             * area_fac
         csite%mmean_can_theiv      (    np) = csite%mmean_can_theiv      (    np)         &
                                             + csite%mmean_can_theiv      (    cp)         &
                                             * area_fac
         csite%mmean_can_theta      (    np) = csite%mmean_can_theta      (    np)         &
                                             + csite%mmean_can_theta      (    cp)         &
                                             * area_fac
         csite%mmean_can_vpdef      (    np) = csite%mmean_can_vpdef      (    np)         &
                                             + csite%mmean_can_vpdef      (    cp)         &
                                             * area_fac
         csite%mmean_can_temp       (    np) = csite%mmean_can_temp       (    np)         &
                                             + csite%mmean_can_temp       (    cp)         &
                                             * area_fac
         csite%mmean_can_shv        (    np) = csite%mmean_can_shv        (    np)         &
                                             + csite%mmean_can_shv        (    cp)         &
                                             * area_fac
         csite%mmean_can_co2        (    np) = csite%mmean_can_co2        (    np)         &
                                             + csite%mmean_can_co2        (    cp)         &
                                             * area_fac
         csite%mmean_can_rhos       (    np) = csite%mmean_can_rhos       (    np)         &
                                             + csite%mmean_can_rhos       (    cp)         &
                                             * area_fac
         csite%mmean_can_dmol       (    np) = csite%mmean_can_dmol       (    np)         &
                                             + csite%mmean_can_dmol       (    cp)         &
                                             * area_fac
         csite%mmean_can_prss       (    np) = csite%mmean_can_prss       (    np)         &
                                             + csite%mmean_can_prss       (    cp)         &
                                             * area_fac
         csite%mmean_gnd_temp       (    np) = csite%mmean_gnd_temp       (    np)         &
                                             + csite%mmean_gnd_temp       (    cp)         &
                                             * area_fac
         csite%mmean_gnd_shv        (    np) = csite%mmean_gnd_shv        (    np)         &
                                             + csite%mmean_gnd_shv        (    cp)         &
                                             * area_fac
         csite%mmean_can_ggnd       (    np) = csite%mmean_can_ggnd       (    np)         &
                                             + csite%mmean_can_ggnd       (    cp)         &
                                             * area_fac
         csite%mmean_sfcw_depth     (    np) = csite%mmean_sfcw_depth     (    np)         &
                                             + csite%mmean_sfcw_depth     (    cp)         &
                                             * area_fac
         csite%mmean_sfcw_energy    (    np) = csite%mmean_sfcw_energy    (    np)         &
                                             + csite%mmean_sfcw_energy    (    cp)         &
                                             * area_fac
         csite%mmean_sfcw_mass      (    np) = csite%mmean_sfcw_mass      (    np)         &
                                             + csite%mmean_sfcw_mass      (    cp)         &
                                             * area_fac
         csite%mmean_sfcw_temp      (    np) = csite%mmean_sfcw_temp      (    np)         &
                                             + csite%mmean_sfcw_temp      (    cp)         &
                                             * area_fac
         csite%mmean_sfcw_fliq      (    np) = csite%mmean_sfcw_fliq      (    np)         &
                                             + csite%mmean_sfcw_fliq      (    cp)         &
                                             * area_fac
         csite%mmean_rshort_gnd     (    np) = csite%mmean_rshort_gnd     (    np)         &
                                             + csite%mmean_rshort_gnd     (    cp)         &
                                             * area_fac
         csite%mmean_par_gnd        (    np) = csite%mmean_par_gnd        (    np)         &
                                             + csite%mmean_par_gnd        (    cp)         &
                                             * area_fac
         csite%mmean_rlong_gnd      (    np) = csite%mmean_rlong_gnd      (    np)         &
                                             + csite%mmean_rlong_gnd      (    cp)         &
                                             * area_fac
         csite%mmean_rlongup        (    np) = csite%mmean_rlongup        (    np)         &
                                             + csite%mmean_rlongup        (    cp)         &
                                             * area_fac
         csite%mmean_parup          (    np) = csite%mmean_parup          (    np)         &
                                             + csite%mmean_parup          (    cp)         &
                                             * area_fac
         csite%mmean_nirup          (    np) = csite%mmean_nirup          (    np)         &
                                             + csite%mmean_nirup          (    cp)         &
                                             * area_fac
         csite%mmean_rshortup       (    np) = csite%mmean_rshortup       (    np)         &
                                             + csite%mmean_rshortup       (    cp)         &
                                             * area_fac
         csite%mmean_rnet           (    np) = csite%mmean_rnet           (    np)         &
                                             + csite%mmean_rnet           (    cp)         &
                                             * area_fac
         csite%mmean_albedo         (    np) = csite%mmean_albedo         (    np)         &
                                             + csite%mmean_albedo         (    cp)         &
                                             * area_fac
         csite%mmean_albedo_par     (    np) = csite%mmean_albedo_par     (    np)         &
                                             + csite%mmean_albedo_par     (    cp)         &
                                             * area_fac
         csite%mmean_albedo_nir     (    np) = csite%mmean_albedo_nir     (    np)         &
                                             + csite%mmean_albedo_nir     (    cp)         &
                                             * area_fac
         csite%mmean_rlong_albedo   (    np) = csite%mmean_rlong_albedo   (    np)         &
                                             + csite%mmean_rlong_albedo   (    cp)         &
                                             * area_fac
         csite%mmean_ustar          (    np) = csite%mmean_ustar          (    np)         &
                                             + csite%mmean_ustar          (    cp)         &
                                             * area_fac
         csite%mmean_tstar          (    np) = csite%mmean_tstar          (    np)         &
                                             + csite%mmean_tstar          (    cp)         &
                                             * area_fac
         csite%mmean_qstar          (    np) = csite%mmean_qstar          (    np)         &
                                             + csite%mmean_qstar          (    cp)         &
                                             * area_fac
         csite%mmean_cstar          (    np) = csite%mmean_cstar          (    np)         &
                                             + csite%mmean_cstar          (    cp)         &
                                             * area_fac
         csite%mmean_carbon_ac      (    np) = csite%mmean_carbon_ac      (    np)         &
                                             + csite%mmean_carbon_ac      (    cp)         &
                                             * area_fac
         csite%mmean_carbon_st      (    np) = csite%mmean_carbon_st      (    np)         &
                                             + csite%mmean_carbon_st      (    cp)         &
                                             * area_fac
         csite%mmean_vapor_gc       (    np) = csite%mmean_vapor_gc       (    np)         &
                                             + csite%mmean_vapor_gc       (    cp)         &
                                             * area_fac
         csite%mmean_vapor_ac       (    np) = csite%mmean_vapor_ac       (    np)         &
                                             + csite%mmean_vapor_ac       (    cp)         &
                                             * area_fac
         csite%mmean_throughfall    (    np) = csite%mmean_throughfall    (    np)         &
                                             + csite%mmean_throughfall    (    cp)         &
                                             * area_fac
         csite%mmean_runoff         (    np) = csite%mmean_runoff         (    np)         &
                                             + csite%mmean_runoff         (    cp)         &
                                             * area_fac
         csite%mmean_drainage       (    np) = csite%mmean_drainage       (    np)         &
                                             + csite%mmean_drainage       (    cp)         &
                                             * area_fac
         csite%mmean_sensible_gc    (    np) = csite%mmean_sensible_gc    (    np)         &
                                             + csite%mmean_sensible_gc    (    cp)         &
                                             * area_fac
         csite%mmean_sensible_ac    (    np) = csite%mmean_sensible_ac    (    np)         &
                                             + csite%mmean_sensible_ac    (    cp)         &
                                             * area_fac
         csite%mmean_qthroughfall   (    np) = csite%mmean_qthroughfall   (    np)         &
                                             + csite%mmean_qthroughfall   (    cp)         &
                                             * area_fac
         csite%mmean_qrunoff        (    np) = csite%mmean_qrunoff        (    np)         &
                                             + csite%mmean_qrunoff        (    cp)         &
                                             * area_fac
         csite%mmean_qdrainage      (    np) = csite%mmean_qdrainage      (    np)         &
                                             + csite%mmean_qdrainage      (    cp)         &
                                             * area_fac
         csite%mmean_A_decomp       (    np) = csite%mmean_A_decomp       (    np)         &
                                             + csite%mmean_A_decomp       (    cp)         &
                                             * area_fac
         csite%mmean_B_decomp       (    np) = csite%mmean_B_decomp       (    np)         &
                                             + csite%mmean_B_decomp       (    cp)         &
                                             * area_fac
         csite%mmean_Af_decomp      (    np) = csite%mmean_Af_decomp      (    np)         &
                                             + csite%mmean_Af_decomp      (    cp)         &
                                             * area_fac
         csite%mmean_Bf_decomp      (    np) = csite%mmean_Bf_decomp      (    np)         &
                                             + csite%mmean_Bf_decomp      (    cp)         &
                                             * area_fac
         csite%mmean_co2_residual   (    np) = csite%mmean_co2_residual   (    np)         &
                                             + csite%mmean_co2_residual   (    cp)         &
                                             * area_fac
         csite%mmean_energy_residual(    np) = csite%mmean_energy_residual(    np)         &
                                             + csite%mmean_energy_residual(    cp)         &
                                             * area_fac
         csite%mmean_water_residual (    np) = csite%mmean_water_residual (    np)         &
                                             + csite%mmean_water_residual (    cp)         &
                                             * area_fac
         csite%mmsqu_rh             (    np) = csite%mmsqu_rh             (    np)         &
                                             + csite%mmsqu_rh             (    cp)         &
                                             * area_fac
         csite%mmsqu_fgc_rh         (    np) = csite%mmsqu_fgc_rh         (    np)         &
                                             + csite%mmsqu_fgc_rh         (    cp)         &
                                             * area_fac
         csite%mmsqu_fsc_rh         (    np) = csite%mmsqu_fsc_rh         (    np)         &
                                             + csite%mmsqu_fsc_rh         (    cp)         &
                                             * area_fac
         csite%mmsqu_stgc_rh        (    np) = csite%mmsqu_stgc_rh        (    np)         &
                                             + csite%mmsqu_stgc_rh        (    cp)         &
                                             * area_fac
         csite%mmsqu_stsc_rh        (    np) = csite%mmsqu_stsc_rh        (    np)         &
                                             + csite%mmsqu_stsc_rh        (    cp)         &
                                             * area_fac
         csite%mmsqu_msc_rh         (    np) = csite%mmsqu_msc_rh         (    np)         &
                                             + csite%mmsqu_msc_rh         (    cp)         &
                                             * area_fac
         csite%mmsqu_ssc_rh         (    np) = csite%mmsqu_ssc_rh         (    np)         &
                                             + csite%mmsqu_ssc_rh         (    cp)         &
                                             * area_fac
         csite%mmsqu_psc_rh         (    np) = csite%mmsqu_psc_rh         (    np)         &
                                             + csite%mmsqu_psc_rh         (    cp)         &
                                             * area_fac
         csite%mmsqu_nep            (    np) = csite%mmsqu_nep            (    np)         &
                                             + csite%mmsqu_nep            (    cp)         &
                                             * area_fac
         csite%mmsqu_rlongup        (    np) = csite%mmsqu_rlongup        (    np)         &
                                             + csite%mmsqu_rlongup        (    cp)         &
                                             * area_fac
         csite%mmsqu_parup          (    np) = csite%mmsqu_parup          (    np)         &
                                             + csite%mmsqu_parup          (    cp)         &
                                             * area_fac
         csite%mmsqu_nirup          (    np) = csite%mmsqu_nirup          (    np)         &
                                             + csite%mmsqu_nirup          (    cp)         &
                                             * area_fac
         csite%mmsqu_rshortup       (    np) = csite%mmsqu_rshortup       (    np)         &
                                             + csite%mmsqu_rshortup       (    cp)         &
                                             * area_fac
         csite%mmsqu_rnet           (    np) = csite%mmsqu_rnet           (    np)         &
                                             + csite%mmsqu_rnet           (    cp)         &
                                             * area_fac
         csite%mmsqu_albedo         (    np) = csite%mmsqu_albedo         (    np)         &
                                             + csite%mmsqu_albedo         (    cp)         &
                                             * area_fac
         csite%mmsqu_ustar          (    np) = csite%mmsqu_ustar          (    np)         &
                                             + csite%mmsqu_ustar          (    cp)         &
                                             * area_fac
         csite%mmsqu_carbon_ac      (    np) = csite%mmsqu_carbon_ac      (    np)         &
                                             + csite%mmsqu_carbon_ac      (    cp)         &
                                             * area_fac
         csite%mmsqu_carbon_st      (    np) = csite%mmsqu_carbon_st      (    np)         &
                                             + csite%mmsqu_carbon_st      (    cp)         &
                                             * area_fac
         csite%mmsqu_vapor_gc       (    np) = csite%mmsqu_vapor_gc       (    np)         &
                                             + csite%mmsqu_vapor_gc       (    cp)         &
                                             * area_fac
         csite%mmsqu_vapor_ac       (    np) = csite%mmsqu_vapor_ac       (    np)         &
                                             + csite%mmsqu_vapor_ac       (    cp)         &
                                             * area_fac
         csite%mmsqu_sensible_gc    (    np) = csite%mmsqu_sensible_gc    (    np)         &
                                             + csite%mmsqu_sensible_gc    (    cp)         &
                                             * area_fac
         csite%mmsqu_sensible_ac    (    np) = csite%mmsqu_sensible_ac    (    np)         &
                                             + csite%mmsqu_sensible_ac    (    cp)         &
                                             * area_fac
         csite%mmean_soil_energy    (  :,np) = csite%mmean_soil_energy    (  :,np)         &
                                             + csite%mmean_soil_energy    (  :,cp)         &
                                             * area_fac
         csite%mmean_soil_mstpot    (  :,np) = csite%mmean_soil_mstpot    (  :,np)         &
                                             + csite%mmean_soil_mstpot    (  :,cp)         &
                                             * area_fac
         csite%mmean_soil_water     (  :,np) = csite%mmean_soil_water     (  :,np)         &
                                             + csite%mmean_soil_water     (  :,cp)         &
                                             * area_fac
         csite%mmean_soil_temp      (  :,np) = csite%mmean_soil_temp      (  :,np)         &
                                             + csite%mmean_soil_temp      (  :,cp)         &
                                             * area_fac
         csite%mmean_soil_fliq      (  :,np) = csite%mmean_soil_fliq      (  :,np)         &
                                             + csite%mmean_soil_fliq      (  :,cp)         &
                                             * area_fac
         csite%mmean_smoist_gg      (  :,np) = csite%mmean_smoist_gg      (  :,np)         &
                                             + csite%mmean_smoist_gg      (  :,cp)         &
                                             * area_fac
         csite%mmean_transloss      (  :,np) = csite%mmean_transloss      (  :,np)         &
                                             + csite%mmean_transloss      (  :,cp)         &
                                             * area_fac
         csite%mmean_sensible_gg    (  :,np) = csite%mmean_sensible_gg    (  :,np)         &
                                             + csite%mmean_sensible_gg    (  :,cp)         &
                                             * area_fac
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Mean diel...                                                                   !
      !------------------------------------------------------------------------------------!
      if (writing_dcyc) then
         csite%qmean_rh             (  :,np) = csite%qmean_rh             (  :,np)         &
                                             + csite%qmean_rh             (  :,cp)         &
                                             * area_fac
         csite%qmean_fgc_rh         (  :,np) = csite%qmean_fgc_rh         (  :,np)         &
                                             + csite%qmean_fgc_rh         (  :,cp)         &
                                             * area_fac
         csite%qmean_fsc_rh         (  :,np) = csite%qmean_fsc_rh         (  :,np)         &
                                             + csite%qmean_fsc_rh         (  :,cp)         &
                                             * area_fac
         csite%qmean_stgc_rh        (  :,np) = csite%qmean_stgc_rh        (  :,np)         &
                                             + csite%qmean_stgc_rh        (  :,cp)         &
                                             * area_fac
         csite%qmean_stsc_rh        (  :,np) = csite%qmean_stsc_rh        (  :,np)         &
                                             + csite%qmean_stsc_rh        (  :,cp)         &
                                             * area_fac
         csite%qmean_msc_rh         (  :,np) = csite%qmean_msc_rh         (  :,np)         &
                                             + csite%qmean_msc_rh         (  :,cp)         &
                                             * area_fac
         csite%qmean_ssc_rh         (  :,np) = csite%qmean_ssc_rh         (  :,np)         &
                                             + csite%qmean_ssc_rh         (  :,cp)         &
                                             * area_fac
         csite%qmean_psc_rh         (  :,np) = csite%qmean_psc_rh         (  :,np)         &
                                             + csite%qmean_psc_rh         (  :,cp)         &
                                             * area_fac
         csite%qmean_nep            (  :,np) = csite%qmean_nep            (  :,np)         &
                                             + csite%qmean_nep            (  :,cp)         &
                                             * area_fac
         csite%qmean_rk4step        (  :,np) = csite%qmean_rk4step        (  :,np)         &
                                             + csite%qmean_rk4step        (  :,cp)         &
                                             * area_fac
         csite%qmean_available_water(  :,np) = csite%qmean_available_water(  :,np)         &
                                             + csite%qmean_available_water(  :,cp)         &
                                             * area_fac
         csite%qmean_veg_displace   (  :,np) = csite%qmean_veg_displace   (  :,np)         &
                                             + csite%qmean_veg_displace   (  :,cp)         &
                                             * area_fac
         csite%qmean_rough          (  :,np) = csite%qmean_rough          (  :,np)         &
                                             + csite%qmean_rough          (  :,cp)         &
                                             * area_fac
         csite%qmean_can_theiv      (  :,np) = csite%qmean_can_theiv      (  :,np)         &
                                             + csite%qmean_can_theiv      (  :,cp)         &
                                             * area_fac
         csite%qmean_can_theta      (  :,np) = csite%qmean_can_theta      (  :,np)         &
                                             + csite%qmean_can_theta      (  :,cp)         &
                                             * area_fac
         csite%qmean_can_vpdef      (  :,np) = csite%qmean_can_vpdef      (  :,np)         &
                                             + csite%qmean_can_vpdef      (  :,cp)         &
                                             * area_fac
         csite%qmean_can_temp       (  :,np) = csite%qmean_can_temp       (  :,np)         &
                                             + csite%qmean_can_temp       (  :,cp)         &
                                             * area_fac
         csite%qmean_can_shv        (  :,np) = csite%qmean_can_shv        (  :,np)         &
                                             + csite%qmean_can_shv        (  :,cp)         &
                                             * area_fac
         csite%qmean_can_co2        (  :,np) = csite%qmean_can_co2        (  :,np)         &
                                             + csite%qmean_can_co2        (  :,cp)         &
                                             * area_fac
         csite%qmean_can_rhos       (  :,np) = csite%qmean_can_rhos       (  :,np)         &
                                             + csite%qmean_can_rhos       (  :,cp)         &
                                             * area_fac
         csite%qmean_can_dmol       (  :,np) = csite%qmean_can_dmol       (  :,np)         &
                                             + csite%qmean_can_dmol       (  :,cp)         &
                                             * area_fac
         csite%qmean_can_prss       (  :,np) = csite%qmean_can_prss       (  :,np)         &
                                             + csite%qmean_can_prss       (  :,cp)         &
                                             * area_fac
         csite%qmean_gnd_temp       (  :,np) = csite%qmean_gnd_temp       (  :,np)         &
                                             + csite%qmean_gnd_temp       (  :,cp)         &
                                             * area_fac
         csite%qmean_gnd_shv        (  :,np) = csite%qmean_gnd_shv        (  :,np)         &
                                             + csite%qmean_gnd_shv        (  :,cp)         &
                                             * area_fac
         csite%qmean_can_ggnd       (  :,np) = csite%qmean_can_ggnd       (  :,np)         &
                                             + csite%qmean_can_ggnd       (  :,cp)         &
                                             * area_fac
         csite%qmean_sfcw_depth     (  :,np) = csite%qmean_sfcw_depth     (  :,np)         &
                                             + csite%qmean_sfcw_depth     (  :,cp)         &
                                             * area_fac
         csite%qmean_sfcw_energy    (  :,np) = csite%qmean_sfcw_energy    (  :,np)         &
                                             + csite%qmean_sfcw_energy    (  :,cp)         &
                                             * area_fac
         csite%qmean_sfcw_mass      (  :,np) = csite%qmean_sfcw_mass      (  :,np)         &
                                             + csite%qmean_sfcw_mass      (  :,cp)         &
                                             * area_fac
         csite%qmean_sfcw_temp      (  :,np) = csite%qmean_sfcw_temp      (  :,np)         &
                                             + csite%qmean_sfcw_temp      (  :,cp)         &
                                             * area_fac
         csite%qmean_sfcw_fliq      (  :,np) = csite%qmean_sfcw_fliq      (  :,np)         &
                                             + csite%qmean_sfcw_fliq      (  :,cp)         &
                                             * area_fac
         csite%qmean_soil_energy    (:,:,np) = csite%qmean_soil_energy    (:,:,np)         &
                                             + csite%qmean_soil_energy    (:,:,cp)         &
                                             * area_fac
         csite%qmean_soil_mstpot    (:,:,np) = csite%qmean_soil_mstpot    (:,:,np)         &
                                             + csite%qmean_soil_mstpot    (:,:,cp)         &
                                             * area_fac
         csite%qmean_soil_water     (:,:,np) = csite%qmean_soil_water     (:,:,np)         &
                                             + csite%qmean_soil_water     (:,:,cp)         &
                                             * area_fac
         csite%qmean_soil_temp      (:,:,np) = csite%qmean_soil_temp      (:,:,np)         &
                                             + csite%qmean_soil_temp      (:,:,cp)         &
                                             * area_fac
         csite%qmean_soil_fliq      (:,:,np) = csite%qmean_soil_fliq      (:,:,np)         &
                                             + csite%qmean_soil_fliq      (:,:,cp)         &
                                             * area_fac
         csite%qmean_rshort_gnd     (  :,np) = csite%qmean_rshort_gnd     (  :,np)         &
                                             + csite%qmean_rshort_gnd     (  :,cp)         &
                                             * area_fac
         csite%qmean_par_gnd        (  :,np) = csite%qmean_par_gnd        (  :,np)         &
                                             + csite%qmean_par_gnd        (  :,cp)         &
                                             * area_fac
         csite%qmean_rlong_gnd      (  :,np) = csite%qmean_rlong_gnd      (  :,np)         &
                                             + csite%qmean_rlong_gnd      (  :,cp)         &
                                             * area_fac
         csite%qmean_rlongup        (  :,np) = csite%qmean_rlongup        (  :,np)         &
                                             + csite%qmean_rlongup        (  :,cp)         &
                                             * area_fac
         csite%qmean_parup          (  :,np) = csite%qmean_parup          (  :,np)         &
                                             + csite%qmean_parup          (  :,cp)         &
                                             * area_fac
         csite%qmean_nirup          (  :,np) = csite%qmean_nirup          (  :,np)         &
                                             + csite%qmean_nirup          (  :,cp)         &
                                             * area_fac
         csite%qmean_rshortup       (  :,np) = csite%qmean_rshortup       (  :,np)         &
                                             + csite%qmean_rshortup       (  :,cp)         &
                                             * area_fac
         csite%qmean_rnet           (  :,np) = csite%qmean_rnet           (  :,np)         &
                                             + csite%qmean_rnet           (  :,cp)         &
                                             * area_fac
         csite%qmean_albedo         (  :,np) = csite%qmean_albedo         (  :,np)         &
                                             + csite%qmean_albedo         (  :,cp)         &
                                             * area_fac
         csite%qmean_albedo_par     (  :,np) = csite%qmean_albedo_par     (  :,np)         &
                                             + csite%qmean_albedo_par     (  :,cp)         &
                                             * area_fac
         csite%qmean_albedo_nir     (  :,np) = csite%qmean_albedo_nir     (  :,np)         &
                                             + csite%qmean_albedo_nir     (  :,cp)         &
                                             * area_fac
         csite%qmean_rlong_albedo   (  :,np) = csite%qmean_rlong_albedo   (  :,np)         &
                                             + csite%qmean_rlong_albedo   (  :,cp)         &
                                             * area_fac
         csite%qmean_ustar          (  :,np) = csite%qmean_ustar          (  :,np)         &
                                             + csite%qmean_ustar          (  :,cp)         &
                                             * area_fac
         csite%qmean_tstar          (  :,np) = csite%qmean_tstar          (  :,np)         &
                                             + csite%qmean_tstar          (  :,cp)         &
                                             * area_fac
         csite%qmean_qstar          (  :,np) = csite%qmean_qstar          (  :,np)         &
                                             + csite%qmean_qstar          (  :,cp)         &
                                             * area_fac
         csite%qmean_cstar          (  :,np) = csite%qmean_cstar          (  :,np)         &
                                             + csite%qmean_cstar          (  :,cp)         &
                                             * area_fac
         csite%qmean_carbon_ac      (  :,np) = csite%qmean_carbon_ac      (  :,np)         &
                                             + csite%qmean_carbon_ac      (  :,cp)         &
                                             * area_fac
         csite%qmean_carbon_st      (  :,np) = csite%qmean_carbon_st      (  :,np)         &
                                             + csite%qmean_carbon_st      (  :,cp)         &
                                             * area_fac
         csite%qmean_vapor_gc       (  :,np) = csite%qmean_vapor_gc       (  :,np)         &
                                             + csite%qmean_vapor_gc       (  :,cp)         &
                                             * area_fac
         csite%qmean_vapor_ac       (  :,np) = csite%qmean_vapor_ac       (  :,np)         &
                                             + csite%qmean_vapor_ac       (  :,cp)         &
                                             * area_fac
         csite%qmean_smoist_gg      (:,:,np) = csite%qmean_smoist_gg      (:,:,np)         &
                                             + csite%qmean_smoist_gg      (:,:,cp)         &
                                             * area_fac
         csite%qmean_throughfall    (  :,np) = csite%qmean_throughfall    (  :,np)         &
                                             + csite%qmean_throughfall    (  :,cp)         &
                                             * area_fac
         csite%qmean_transloss      (:,:,np) = csite%qmean_transloss      (:,:,np)         &
                                             + csite%qmean_transloss      (:,:,cp)         &
                                             * area_fac
         csite%qmean_runoff         (  :,np) = csite%qmean_runoff         (  :,np)         &
                                             + csite%qmean_runoff         (  :,cp)         &
                                             * area_fac
         csite%qmean_drainage       (  :,np) = csite%qmean_drainage       (  :,np)         &
                                             + csite%qmean_drainage       (  :,cp)         &
                                             * area_fac
         csite%qmean_sensible_gc    (  :,np) = csite%qmean_sensible_gc    (  :,np)         &
                                             + csite%qmean_sensible_gc    (  :,cp)         &
                                             * area_fac
         csite%qmean_sensible_ac    (  :,np) = csite%qmean_sensible_ac    (  :,np)         &
                                             + csite%qmean_sensible_ac    (  :,cp)         &
                                             * area_fac
         csite%qmean_sensible_gg    (:,:,np) = csite%qmean_sensible_gg    (:,:,np)         &
                                             + csite%qmean_sensible_gg    (:,:,cp)         &
                                             * area_fac
         csite%qmean_qthroughfall   (  :,np) = csite%qmean_qthroughfall   (  :,np)         &
                                             + csite%qmean_qthroughfall   (  :,cp)         &
                                             * area_fac
         csite%qmean_qrunoff        (  :,np) = csite%qmean_qrunoff        (  :,np)         &
                                             + csite%qmean_qrunoff        (  :,cp)         &
                                             * area_fac
         csite%qmean_qdrainage      (  :,np) = csite%qmean_qdrainage      (  :,np)         &
                                             + csite%qmean_qdrainage      (  :,cp)         &
                                             * area_fac
         csite%qmsqu_rh             (  :,np) = csite%qmsqu_rh             (  :,np)         &
                                             + csite%qmsqu_rh             (  :,cp)         &
                                             * area_fac
         csite%qmsqu_fgc_rh         (  :,np) = csite%qmsqu_fgc_rh         (  :,np)         &
                                             + csite%qmsqu_fgc_rh         (  :,cp)         &
                                             * area_fac
         csite%qmsqu_fsc_rh         (  :,np) = csite%qmsqu_fsc_rh         (  :,np)         &
                                             + csite%qmsqu_fsc_rh         (  :,cp)         &
                                             * area_fac
         csite%qmsqu_stgc_rh        (  :,np) = csite%qmsqu_stgc_rh        (  :,np)         &
                                             + csite%qmsqu_stgc_rh        (  :,cp)         &
                                             * area_fac
         csite%qmsqu_stsc_rh        (  :,np) = csite%qmsqu_stsc_rh        (  :,np)         &
                                             + csite%qmsqu_stsc_rh        (  :,cp)         &
                                             * area_fac
         csite%qmsqu_msc_rh         (  :,np) = csite%qmsqu_msc_rh         (  :,np)         &
                                             + csite%qmsqu_msc_rh         (  :,cp)         &
                                             * area_fac
         csite%qmsqu_ssc_rh         (  :,np) = csite%qmsqu_ssc_rh         (  :,np)         &
                                             + csite%qmsqu_ssc_rh         (  :,cp)         &
                                             * area_fac
         csite%qmsqu_psc_rh         (  :,np) = csite%qmsqu_psc_rh         (  :,np)         &
                                             + csite%qmsqu_psc_rh         (  :,cp)         &
                                             * area_fac
         csite%qmsqu_nep            (  :,np) = csite%qmsqu_nep            (  :,np)         &
                                             + csite%qmsqu_nep            (  :,cp)         &
                                             * area_fac
         csite%qmsqu_rlongup        (  :,np) = csite%qmsqu_rlongup        (  :,np)         &
                                             + csite%qmsqu_rlongup        (  :,cp)         &
                                             * area_fac
         csite%qmsqu_parup          (  :,np) = csite%qmsqu_parup          (  :,np)         &
                                             + csite%qmsqu_parup          (  :,cp)         &
                                             * area_fac
         csite%qmsqu_nirup          (  :,np) = csite%qmsqu_nirup          (  :,np)         &
                                             + csite%qmsqu_nirup          (  :,cp)         &
                                             * area_fac
         csite%qmsqu_rshortup       (  :,np) = csite%qmsqu_rshortup       (  :,np)         &
                                             + csite%qmsqu_rshortup       (  :,cp)         &
                                             * area_fac
         csite%qmsqu_rnet           (  :,np) = csite%qmsqu_rnet           (  :,np)         &
                                             + csite%qmsqu_rnet           (  :,cp)         &
                                             * area_fac
         csite%qmsqu_albedo         (  :,np) = csite%qmsqu_albedo         (  :,np)         &
                                             + csite%qmsqu_albedo         (  :,cp)         &
                                             * area_fac
         csite%qmsqu_ustar          (  :,np) = csite%qmsqu_ustar          (  :,np)         &
                                             + csite%qmsqu_ustar          (  :,cp)         &
                                             * area_fac
         csite%qmsqu_carbon_ac      (  :,np) = csite%qmsqu_carbon_ac      (  :,np)         &
                                             + csite%qmsqu_carbon_ac      (  :,cp)         &
                                             * area_fac
         csite%qmsqu_carbon_st      (  :,np) = csite%qmsqu_carbon_st      (  :,np)         &
                                             + csite%qmsqu_carbon_st      (  :,cp)         &
                                             * area_fac
         csite%qmsqu_vapor_gc       (  :,np) = csite%qmsqu_vapor_gc       (  :,np)         &
                                             + csite%qmsqu_vapor_gc       (  :,cp)         &
                                             * area_fac
         csite%qmsqu_vapor_ac       (  :,np) = csite%qmsqu_vapor_ac       (  :,np)         &
                                             + csite%qmsqu_vapor_ac       (  :,cp)         &
                                             * area_fac
         csite%qmsqu_sensible_gc    (  :,np) = csite%qmsqu_sensible_gc    (  :,np)         &
                                             + csite%qmsqu_sensible_gc    (  :,cp)         &
                                             * area_fac
         csite%qmsqu_sensible_ac    (  :,np) = csite%qmsqu_sensible_ac    (  :,np)         &
                                             + csite%qmsqu_sensible_ac    (  :,cp)         &
                                             * area_fac
      end if
      !------------------------------------------------------------------------------------!


      return
   end subroutine increment_patch_vars
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine mormalises the canopy air properties integrated across patches.   !
   !---------------------------------------------------------------------------------------!
   subroutine norm_canopy_air_props(csite,ipa,cb_enthalpy,can_exner,cb_mass,cb_molar       &
                                   ,cb_water,cb_co2)
      use ed_state_vars, only : sitetype                 & ! structure
                              , patchtype                ! ! structure
      use therm_lib    , only : exner2press              & ! function
                              , hq2temp                  & ! function
                              , extemp2theta             & ! function
                              , thetaeiv                 & ! function
                              , vpdefil                  & ! function
                              , idealdenssh              & ! function
                              , idealdmolsh              ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype) , target     :: csite
      integer        , intent(in) :: ipa
      real           , intent(in) :: cb_enthalpy
      real           , intent(in) :: can_exner
      real           , intent(in) :: cb_mass
      real           , intent(in) :: cb_molar
      real           , intent(in) :: cb_water
      real           , intent(in) :: cb_co2
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype), pointer    :: cpatch
      real                        :: can_enthalpy
      real                        :: can_rvap
      !------------------------------------------------------------------------------------!


      !----- Alias. -----------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!


      !------ Find air density and the specific properties. -------------------------------!
      can_enthalpy         = cb_enthalpy / cb_mass
      csite%can_shv  (ipa) = cb_water    / cb_mass
      csite%can_co2  (ipa) = cb_co2      / cb_molar
      !------ Find temperature first, as most variables depend upon it. -------------------!
      csite%can_temp (ipa) = hq2temp(can_enthalpy,csite%can_shv(ipa),.true.)
      !------ Find pressure, and ideal-gas densities (mass and molar). --------------------!
      csite%can_prss (ipa) = exner2press(can_exner)
      csite%can_rhos (ipa) = idealdenssh(csite%can_prss(ipa),csite%can_temp(ipa)           &
                                        ,csite%can_shv (ipa),csite%can_shv (ipa) )
      csite%can_dmol (ipa) = idealdmolsh(csite%can_prss(ipa),csite%can_temp(ipa)           &
                                        ,csite%can_shv (ipa))
      !------ Find potential temperatures. ------------------------------------------------!
      can_rvap             = csite%can_shv(ipa) / (1.0 - csite%can_shv(ipa))
      csite%can_theta(ipa) = extemp2theta(can_exner,csite%can_temp(ipa))
      csite%can_theiv(ipa) = thetaeiv( csite%can_theta(ipa), csite%can_prss (ipa)          &
                                     , csite%can_temp (ipa), can_rvap                      &
                                     , can_rvap                                   )
      !------ Update vapour pressure deficit. ---------------------------------------------!
      csite%can_vpdef(ipa) = vpdefil( csite%can_prss(ipa), csite%can_temp(ipa)             &
                                    , csite%can_shv (ipa), .true.              )
      !------------------------------------------------------------------------------------!



      return
   end subroutine norm_canopy_air_props
   !=======================================================================================!
   !=======================================================================================!








   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will populate the disturbed patch with the cohorts that were      !
   ! disturbed but did not go extinct.                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine insert_survivors(csite,np,cp,new_lu,area_fac,mindbh_harvest)
      use ed_state_vars       , only : sitetype                      & ! structure
                                     , patchtype                     ! ! structure
      use ed_max_dims         , only : n_pft                         ! ! intent(in)
      use mortality           , only : survivorship                  ! ! function
      use update_derived_utils, only : update_cohort_extensive_props ! ! subroutine
      use pft_coms            , only : negligible_nplant             ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                  , target        :: csite
      integer                         , intent(in)    :: new_lu
      integer                         , intent(in)    :: np
      integer                         , intent(in)    :: cp
      real          , dimension(n_pft), intent(in)    :: mindbh_harvest
      real                            , intent(in)    :: area_fac
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)                 , pointer       :: cpatch
      type(patchtype)                 , pointer       :: npatch
      type(patchtype)                 , pointer       :: tpatch
      logical        , dimension(:)   , allocatable   :: mask
      real           , dimension(:)   , allocatable   :: survival_fac
      integer                                         :: ipft
      integer                                         :: ico
      integer                                         :: nco
      integer                                         :: addco
      real                                            :: n_survivors
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! cpatch => Old patches we are looping                                               !
      ! npatch => New transition patch                                                     !
      ! tpatch => temporary patch                                                          !
      !------------------------------------------------------------------------------------!
      cpatch => csite%patch(cp)
      npatch => csite%patch(np)
      nullify(tpatch)
      allocate(tpatch)

      !----- Mask: flag to decide whether the cohort survived or not. ---------------------!
      if (cpatch%ncohorts > 0) then
         allocate(mask        (cpatch%ncohorts))
         allocate(survival_fac(cpatch%ncohorts))
         mask(:)         = .false.
         survival_fac(:) = 0.

         survivalloop: do ico = 1,cpatch%ncohorts
            ipft              = cpatch%pft(ico)
            survival_fac(ico) = survivorship(new_lu,csite%dist_type(cp),mindbh_harvest     &
                                            ,cpatch,ico) * area_fac
            n_survivors       = cpatch%nplant(ico) * survival_fac(ico)

            !----- If something survived, make a new cohort. ------------------------------!
            if ( n_survivors >= negligible_nplant(ipft) ) then
               mask(ico)         = .true.
            else 
               survival_fac(ico) = 0.0
            end if
            !------------------------------------------------------------------------------!
         end do survivalloop
         addco = count(mask)
      else
         addco = 0
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     If the new patch has received survivors from a donor already, then it should   !
      ! have cohorts.  So the temporary patch vector will be the sum of the new cohorts    !
      ! found here, plus those already applied previously in the loop calling this sub-    !
      ! routine.                                                                           !
      !------------------------------------------------------------------------------------!
      if (npatch%ncohorts > 0) then
         nco = npatch%ncohorts
         call allocate_patchtype(tpatch,addco + npatch%ncohorts)
         call copy_patchtype(npatch,tpatch,1,npatch%ncohorts,1,npatch%ncohorts)
         call deallocate_patchtype(npatch)
      else
         nco = 0
         call allocate_patchtype(tpatch,addco)
      end if
      !------------------------------------------------------------------------------------!


      cohortloop: do ico = 1,cpatch%ncohorts
         !----- If mask is true, at least some of this cohort survived. -------------------!
         if (mask(ico)) then
            nco = nco + 1
            call copy_patchtype(cpatch,tpatch,ico,ico,nco,nco)

            !------------------------------------------------------------------------------!
            !    Scale the total area based on the new population density and new area.    !
            !------------------------------------------------------------------------------!
            call update_cohort_extensive_props(tpatch,nco,nco,survival_fac(ico))
            !------------------------------------------------------------------------------!

            !----- Make mortality rate due to disturbance zero to avoid double counting. --!
            tpatch%mort_rate(6,nco) = 0.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do cohortloop
      !------------------------------------------------------------------------------------!

      !----- Copy the temporary patch into the newpatch. ----------------------------------!
      call allocate_patchtype(npatch,tpatch%ncohorts)
      call copy_patchtype(tpatch,npatch,1,tpatch%ncohorts,1,tpatch%ncohorts)
      call deallocate_patchtype(tpatch)

      deallocate(tpatch)
      if (allocated(mask        )) deallocate(mask        )
      if (allocated(survival_fac)) deallocate(survival_fac)


      return
   end subroutine insert_survivors
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine updates the harvest and litter pools after a disturbance takes    !
   ! place.                                                                                !
   !---------------------------------------------------------------------------------------!
   subroutine accum_dist_harv_litt(cpoly,isi,census_flag,np,cp,new_lu,area_fac             &
                                  ,mindbh_harvest)
      use ed_state_vars, only : sitetype              & ! structure
                              , patchtype             & ! structure
                              , polygontype           ! ! structure
      use disturb_coms , only : include_fire          & ! intent(in)
                              , cl_fleaf_harvest      & ! intent(in)
                              , cl_fstorage_harvest   & ! intent(in)
                              , f_combusted_fast_c    & ! intent(in)
                              , f_combusted_struct_c  & ! intent(in)
                              , f_combusted_fast_n    & ! intent(in)
                              , f_combusted_struct_n  ! ! intent(in)
      use ed_max_dims  , only : n_pft                 ! ! intent(in)
      use pft_coms     , only : c2n_storage           & ! intent(in)
                              , c2n_leaf              & ! intent(in)
                              , c2n_stem              & ! intent(in)
                              , l2n_stem              & ! intent(in)
                              , agf_bs                & ! intent(in)
                              , f_labile_leaf         & ! intent(in)
                              , f_labile_stem         ! ! intent(in)
      use mortality    , only : survivorship          ! ! function

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(polygontype)                  , target     :: cpoly
      integer                            , intent(in) :: isi
      integer                            , intent(in) :: census_flag
      integer                            , intent(in) :: np
      integer                            , intent(in) :: cp
      real             , dimension(n_pft), intent(in) :: mindbh_harvest
      integer                            , intent(in) :: new_lu
      real                               , intent(in) :: area_fac
      !----- Local variables. -------------------------------------------------------------!
      type(sitetype)                     , pointer    :: csite
      type(patchtype)                    , pointer    :: cpatch
      type(patchtype)                    , pointer    :: npatch
      integer                                         :: ico
      integer                                         :: ipft
      integer                                         :: bdbh
      real                                            :: a_bfast_before
      real                                            :: a_bstruct_before
      real                                            :: a_bstorage_before
      real                                            :: b_bfast_before
      real                                            :: b_bstruct_before
      real                                            :: b_bstorage_before
      real                                            :: a_bfast_remove
      real                                            :: a_bstruct_remove
      real                                            :: a_bstorage_remove
      real                                            :: a_bfast_remain
      real                                            :: a_bstruct_remain
      real                                            :: a_bstorage_remain
      real                                            :: a_bcrop_harvest
      real                                            :: a_blogging_harvest
      real                                            :: a_bcombusted_fuel
      real                                            :: agb_cut
      real                                            :: agb_mort
      real                                            :: ba_cut
      real                                            :: ba_mort
      real                                            :: a_fast_litter
      real                                            :: a_struct_litter
      real                                            :: a_struct_lignin
      real                                            :: a_fast_litter_n
      real                                            :: a_struct_litter_n
      real                                            :: b_fast_litter
      real                                            :: b_struct_litter
      real                                            :: b_struct_lignin
      real                                            :: b_fast_litter_n
      real                                            :: b_struct_litter_n
      real                                            :: a_struct_cohort
      real                                            :: b_struct_cohort
      real                                            :: survival_fac
      real                                            :: a_fast_combusted
      real                                            :: a_struct_combusted
      real                                            :: a_lignin_combusted
      real                                            :: a_fast_combusted_n
      real                                            :: a_struct_combusted_n
      !------------------------------------------------------------------------------------!


      !---- Initialise the non-scaled litter pools. ---------------------------------------!
      a_fast_litter        = 0.0
      b_fast_litter        = 0.0
      a_struct_litter      = 0.0
      b_struct_litter      = 0.0
      a_struct_lignin      = 0.0
      b_struct_lignin      = 0.0
      a_fast_litter_n      = 0.0
      b_fast_litter_n      = 0.0
      a_struct_litter_n    = 0.0
      b_struct_litter_n    = 0.0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! csite  => current site.                                                            !
      ! cpatch => contributing patch.                                                      !
      ! npatch => new patch.                                                               !
      !------------------------------------------------------------------------------------!
      csite  => cpoly%site(isi)
      cpatch => csite%patch(cp)
      npatch => csite%patch(np)
      !------------------------------------------------------------------------------------!


      do ico = 1,cpatch%ncohorts
         ipft = cpatch%pft(ico)
         bdbh = max(0,min( int(cpatch%dbh(ico) * 0.1), 10)) + 1

         !----- Split biomass components that are labile or structural. -------------------!
         a_bfast_before    = f_labile_leaf(ipft) * cpatch%bleaf(ico)                       &
                           + f_labile_stem(ipft)                                           &
                           * (cpatch%bsapwooda(ico)+cpatch%bbarka(ico)+cpatch%bdeada(ico))
         b_bfast_before    = f_labile_leaf(ipft) * cpatch%broot(ico)                       &
                           + f_labile_stem(ipft)                                           &
                           * (cpatch%bsapwoodb(ico)+cpatch%bbarkb(ico)+cpatch%bdeadb(ico))
         a_bstruct_before  = (1.0 - f_labile_leaf(ipft)) * cpatch%bleaf(ico)               &
                           + (1.0 - f_labile_stem(ipft))                                   &
                           * (cpatch%bsapwooda(ico)+cpatch%bbarka(ico)+cpatch%bdeada(ico))
         b_bstruct_before  = (1.0 - f_labile_leaf(ipft)) * cpatch%broot(ico)               &
                           + (1.0 - f_labile_stem(ipft))                                   &
                           * (cpatch%bsapwoodb(ico)+cpatch%bbarkb(ico)+cpatch%bdeadb(ico))
         a_bstorage_before =        agf_bs(ipft)  * cpatch%bstorage(ico)
         b_bstorage_before = (1.0 - agf_bs(ipft)) * cpatch%bstorage(ico)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the biomass that remains in patch upon disturbance (as opposed to      !
         ! being removed from the patch).  In the case of logging or forest plantation, we !
         ! assume that the timber biomass is harvested and removed, and the remaining      !
         ! biomass becomes necromass.  In the case of croplands, we may harvest leaves and !
         ! non-structural carbon.  For the other disturbance types, we assume that every-  !
         ! thing remains in the patch.                                                     !
         !---------------------------------------------------------------------------------!
         select case(new_lu)
         case (8)
            !------ Agriculture. Harvest living tissues and NSC. --------------------------!
            a_bfast_remove     = f_labile_leaf(ipft) * cl_fleaf_harvest * cpatch%bleaf(ico)
            a_bstruct_remove   = ( 1.0 - f_labile_leaf(ipft) )                             &
                               * cl_fleaf_harvest * cpatch%bleaf(ico)
            a_bstorage_remove  = cl_fstorage_harvest * agf_bs(ipft) * cpatch%bstorage(ico)
            a_bcrop_harvest    = cl_fleaf_harvest    * cpatch%bleaf   (ico)                & 
                               + cl_fstorage_harvest * agf_bs(ipft) * cpatch%bstorage(ico)
            a_blogging_harvest = 0.0
            a_bcombusted_fuel  = 0.0
            !------------------------------------------------------------------------------!
         case (2,6)
            !------ Felling / forest plantation.  Harvest commercial timber. --------------!
            a_bfast_remove     =         f_labile_stem(ipft)   * cpatch%btimber(ico)
            a_bstruct_remove   = ( 1.0 - f_labile_stem(ipft) ) * cpatch%btimber(ico)
            a_bstorage_remove  = 0.0
            a_bcrop_harvest    = 0.0
            a_blogging_harvest = cpatch%btimber(ico)
            a_bcombusted_fuel  = 0.0
            !------------------------------------------------------------------------------!
         case (4)
            !------------------------------------------------------------------------------!
            !    Fire.  For the time being, assume all "fast" AG biomass is lost, and      !
            ! that a fraction of the structural biomass is lost through combustion.        !
            !------------------------------------------------------------------------------!
            a_bfast_remove   = ( f_labile_leaf(ipft) * cpatch%bleaf(ico)                   &
                               + f_labile_stem(ipft)                                       &
                               * ( cpatch%bsapwooda(ico)                                   &
                                 + cpatch%bbarka   (ico) + cpatch%bdeada (ico)) )          &
                               * f_combusted_fast_c
            a_bstruct_remove = ( (1.0-f_labile_leaf(ipft)) * cpatch%bleaf(ico)             &
                               + (1.0-f_labile_stem(ipft))                                 &
                               * ( cpatch%bsapwooda(ico)                                   &
                                 + cpatch%bbarka   (ico) + cpatch%bdeada (ico)) )          &
                               * f_combusted_struct_c
            a_bstorage_remove  = agf_bs(ipft) * cpatch%bstorage(ico)
            a_bcrop_harvest    = 0.0
            a_blogging_harvest = 0.0
            a_bcombusted_fuel  = a_bfast_remove + a_bstruct_remove + a_bstorage_remove
         case default
            !------ Other types.  Everything remains in. ----------------------------------!
            a_bfast_remove     = 0.0
            a_bstruct_remove   = 0.0
            a_bstorage_remove  = 0.0
            a_bcrop_harvest    = 0.0
            a_blogging_harvest = 0.0
            a_bcombusted_fuel  = 0.0
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


         !----- Find the remaining biomass that will go to different soil carbon pools. ---!
         a_bfast_remain    = a_bfast_before    - a_bfast_remove
         a_bstruct_remain  = a_bstruct_before  - a_bstruct_remove
         a_bstorage_remain = a_bstorage_before - a_bstorage_remove
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find biomass loss due to cutting (logging) or mortality.  This is slightly  !
         ! different than crop_harvest/logging_harvest because it corresponds to the loss  !
         ! of biomass due to removing or mortality.  In the case of logging, for example,  !
         ! not all biomass is removed from the patch, but both the removed and the remain- !
         ! ing contribute to biomass loss.                                                 !
         !---------------------------------------------------------------------------------!
         select case(new_lu)
         case (2)
            !------ Forest plantation, everything is cut. ---------------------------------!
            agb_cut  = cpatch%agb(ico)
            agb_mort = 0.
            ba_cut   = cpatch%basarea(ico)
            ba_mort  = 0.
            !------------------------------------------------------------------------------!
         case (6)
            !------------------------------------------------------------------------------!
            !     Logging, we must decide whether the losses are due to cutting or due to  !
            ! mortality.                                                                   !
            !------------------------------------------------------------------------------!
            if (cpatch%dbh(ico) >= mindbh_harvest(ipft)) then
               agb_cut          = cpatch%agb(ico)
               agb_mort         = 0.
               ba_cut           = cpatch%basarea(ico)
               ba_mort          = 0.
            else
               agb_cut          = 0.
               agb_mort         = cpatch%agb(ico)
               ba_cut           = 0.
               ba_mort          = cpatch%basarea(ico)
            end if
            !------------------------------------------------------------------------------!
         case default
            !------ Other types.  Everything is assumed mortality. ------------------------! 
            agb_cut          = 0.
            agb_mort         = cpatch%agb(ico)
            ba_cut           = 0.
            ba_mort          = cpatch%basarea(ico)
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


         !----- Find survivorship. --------------------------------------------------------!
         survival_fac  = survivorship(new_lu,csite%dist_type(cp),mindbh_harvest,cpatch,ico)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Update harvest pool.                                                        !
         !---------------------------------------------------------------------------------!
         cpoly%crop_harvest(isi)    = cpoly%crop_harvest(isi)                              &
                                    + (1.-survival_fac) * cpatch%nplant(ico)               &
                                    * a_bcrop_harvest * area_fac * csite%area(np)
         cpoly%logging_harvest(isi) = cpoly%logging_harvest(isi)                           &
                                    + (1.-survival_fac) * cpatch%nplant(ico)               &
                                    * a_blogging_harvest  * area_fac * csite%area(np)
         cpoly%combusted_fuel (isi) = cpoly%combusted_fuel (isi)                           &
                                    + (1.-survival_fac) * cpatch%nplant(ico)               &
                                    * a_bcombusted_fuel   * area_fac * csite%area(np)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Update basal area and above-ground biomass loss due to logging or other    !
         ! types of mortality.                                                             !
         !---------------------------------------------------------------------------------!
         if (census_flag == 0 .or. cpatch%first_census(ico) == 1) then
            cpoly%agb_cut        (ipft,bdbh,isi) = cpoly%agb_cut(ipft,bdbh,isi)            &
                                                 + (1.-survival_fac) * cpatch%nplant(ico)  &
                                                 * agb_cut * area_fac * csite%area(np)
            cpoly%agb_mort       (ipft,bdbh,isi) = cpoly%agb_mort(ipft,bdbh,isi)           &
                                                 + (1.-survival_fac) * cpatch%nplant(ico)  &
                                                 * agb_mort * area_fac * csite%area(np)
            cpoly%basal_area_cut (ipft,bdbh,isi) = cpoly%basal_area_cut(ipft,bdbh,isi)     &
                                                 + (1.-survival_fac) * cpatch%nplant(ico)  &
                                                 * ba_cut * area_fac * csite%area(np)
            cpoly%basal_area_mort(ipft,bdbh,isi) = cpoly%basal_area_mort(ipft,bdbh,isi)    &
                                                 + (1.-survival_fac) * cpatch%nplant(ico)  &
                                                 * ba_mort * area_fac * csite%area(np)
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Update necromass pools (both carbon and nitrogen).  For nitrogen, we assume !
         ! that all inputs go to the fast pool, so nitrogen is conserved.  In the case     !
         ! the new patch is a burnt patch, we also must account for differences in C and N !
         ! volatilisation fractions.                                                       !
         !---------------------------------------------------------------------------------!
         a_fast_litter   = a_fast_litter                                                   &
                         + (1. - survival_fac) * cpatch%nplant(ico)                        &
                         * ( a_bfast_remain + a_bstorage_remain )
         b_fast_litter   = b_fast_litter                                                   &
                         + (1. - survival_fac) * cpatch%nplant(ico)                        &
                         * ( b_bfast_before + b_bstorage_before )
         select case (new_lu)
         case (4)
            !----- Add labile component. --------------------------------------------------!
            if (f_combusted_fast_c > 0.0) then
               a_fast_litter_n = a_fast_litter_n                                           &
                               + (1. - survival_fac) * cpatch%nplant(ico)                  &
                               * ( a_bfast_remain    / c2n_leaf(ipft)                      &
                                 + a_bstorage_remain / c2n_storage        )                &
                               * f_combusted_fast_n / f_combusted_fast_c
            end if
            !------------------------------------------------------------------------------!


            !----- Add lignified component. -----------------------------------------------!
            if (f_combusted_struct_c > 0.0) then
               a_struct_litter_n = a_struct_litter_n                                       &
                                 + (1. - survival_fac) * cpatch%nplant(ico)                &
                                 * a_bstruct_remain     / c2n_stem(ipft)                   &
                                 * f_combusted_struct_n / f_combusted_struct_c
            end if
            !------------------------------------------------------------------------------!
         case default
            !------ Other cases, everything goes to litter. -------------------------------!
            a_fast_litter_n = a_fast_litter_n                                              &
                            + (1. - survival_fac) * cpatch%nplant(ico)                     &
                            * ( a_bfast_remain    / c2n_leaf(ipft)                         &
                              + a_bstorage_remain / c2n_storage       )
            a_struct_litter_n = a_struct_litter_n                                          &
                              + (1. - survival_fac) * cpatch%nplant(ico)                   &
                              *  a_bstruct_remain  / c2n_stem(ipft)
            !------------------------------------------------------------------------------!
         end select
         b_fast_litter_n   = b_fast_litter_n                                               &
                           + (1. - survival_fac) * cpatch%nplant(ico)                      &
                           * ( b_bfast_before    / c2n_leaf(ipft)                          &
                             + b_bstorage_before / c2n_storage        )
         b_struct_litter_n = b_struct_litter_n                                             &
                           + (1. - survival_fac) * cpatch%nplant(ico)                      &
                           * b_bstruct_before   / c2n_stem(ipft) 
         a_struct_cohort   = (1. - survival_fac) * cpatch%nplant(ico) * a_bstruct_remain
         b_struct_cohort   = (1. - survival_fac) * cpatch%nplant(ico) * b_bstruct_before
         a_struct_litter   = a_struct_litter + a_struct_cohort
         b_struct_litter   = b_struct_litter + b_struct_cohort
         a_struct_lignin   = a_struct_lignin + a_struct_cohort * l2n_stem / c2n_stem(ipft)
         b_struct_lignin   = b_struct_lignin + b_struct_cohort * l2n_stem / c2n_stem(ipft)
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Check whether to remove carbon from the pools as combusted fuels.  This is only !
      ! done if this is a burnt patch and if we are using the new fire scheme.  For the    !
      ! time being, we fix a fraction of fast and structural soil carbon that is above     !
      ! ground.  In the future we may split the soil pools into above- and below-ground    !
      ! to make estimates consistent with the contribution of individuals that do not      !
      ! follow the standard tropical allometric parameters.                                !
      !------------------------------------------------------------------------------------!
      if (new_lu == 4 .and. include_fire == 3) then
         a_fast_combusted     = f_combusted_fast_c   * csite%fast_grnd_C      (np)
         a_fast_combusted_n   = f_combusted_fast_n   * csite%fast_grnd_N      (np)
         a_struct_combusted   = f_combusted_struct_c * csite%structural_grnd_C(np)
         a_lignin_combusted   = f_combusted_struct_c * csite%structural_grnd_L(np)
         a_struct_combusted_n = f_combusted_struct_n * csite%structural_grnd_N(np)
      else
         a_fast_combusted     = 0.0
         a_struct_combusted   = 0.0
         a_lignin_combusted   = 0.0
         a_fast_combusted_n   = 0.0
         a_struct_combusted_n = 0.0
      end if
      !------------------------------------------------------------------------------------!



      !----- Load disturbance litter directly into carbon and N pools. --------------------!
      csite%fast_grnd_C      (np) = csite%fast_grnd_C      (np)                            &
                                  + ( a_fast_litter     - a_fast_combusted     ) * area_fac
      csite%fast_soil_C      (np) = csite%fast_soil_C      (np)                            &
                                  + b_fast_litter * area_fac
      csite%structural_grnd_C(np) = csite%structural_grnd_C(np)                            &
                                  + ( a_struct_litter   - a_struct_combusted   ) * area_fac
      csite%structural_soil_C(np) = csite%structural_soil_C(np)                            &
                                  + b_struct_litter * area_fac
      csite%structural_grnd_L(np) = csite%structural_grnd_L(np)                            &
                                  + ( a_struct_lignin   - a_lignin_combusted   ) * area_fac
      csite%structural_soil_L(np) = csite%structural_soil_L(np)                            &
                                  + b_struct_lignin * area_fac
      csite%fast_grnd_N      (np) = csite%fast_grnd_N      (np)                            &
                                  + ( a_fast_litter_n   - a_fast_combusted_n   ) * area_fac
      csite%structural_grnd_N(np) = csite%structural_grnd_N(np)                            &
                                  + ( a_struct_litter_n - a_struct_combusted_n ) * area_fac
      csite%structural_soil_N(np) = csite%structural_soil_N(np)                            &
                                  + b_struct_litter_n * area_fac
      !------------------------------------------------------------------------------------!

      return
   end subroutine accum_dist_harv_litt
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    Add a cohort of the appropriate PFT type to populate a plantation/cropland/pasture !
   ! patch.                                                                                !
   !---------------------------------------------------------------------------------------!
   subroutine plant_patch(csite,np,mzg,pft,density,ntext_soil,height_factor,lsl)
      use ed_state_vars  , only : sitetype                 & ! structure
                                , patchtype                ! ! structure
      use pft_coms       , only : hgt_min                  & ! intent(in)
                                , hgt_max                  & ! intent(in)
                                , agf_bs                   & ! intent(in)
                                , dbh_bigleaf              ! ! intent(in)
      use ed_misc_coms   , only : ibigleaf                 ! ! intent(in)
      use fuse_fiss_utils, only : sort_cohorts             ! ! sub-routine
      use ed_therm_lib   , only : calc_veg_hcap            ! ! function
      use consts_coms    , only : t3ple                    & ! intent(in)
                                , pio4                     ! ! intent(in)
      use allometry      , only : h2dbh                    & ! function
                                , size2bd                   & ! function
                                , size2bt                  & ! function
                                , size2xb                  & ! function
                                , area_indices             & ! function
                                , ed_balive                & ! function
                                , ed_biomass               ! ! function
      use ed_max_dims    , only : n_pft                    ! ! intent(in)
      use phenology_aux  , only : pheninit_balive_bstorage ! ! intent(in)
      use therm_lib      , only : cmtl2uext                ! ! intent(in)
      use stable_cohorts , only : is_resolvable            ! ! sub-routine
      use ed_type_init   , only : init_ed_cohort_vars      ! ! sub-routine
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                  , target     :: csite
      integer                         , intent(in) :: mzg
      integer                         , intent(in) :: np
      integer                         , intent(in) :: pft
      integer                         , intent(in) :: lsl
      integer       , dimension(mzg)  , intent(in) :: ntext_soil
      real                            , intent(in) :: density
      real                            , intent(in) :: height_factor
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)                 , pointer    :: cpatch
      type(patchtype)                 , pointer    :: tpatch
      real                                         :: bdead
      integer                                      :: nc
      !------------------------------------------------------------------------------------!



      cpatch => csite%patch(np)


      !------------------------------------------------------------------------------------!
      !      Reallocate the current cohort with an extra space for the planted plantation  !
      ! cohort.  If the patch was previously empty, we simply create the first cohort.     !
      !------------------------------------------------------------------------------------!
      nc = cpatch%ncohorts + 1
      if (cpatch%ncohorts > 0) then
         nullify(tpatch)
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
      end if

      cpatch%ncohorts = nc
      csite%paco_n(np)= nc

      !------------------------------------------------------------------------------------!
      !      We can add all the biomass needed in a single cohort, even if it exceeds the  !
      ! safe range (i.e., with maximum LAI exceeding the unity).  In case it is a large    !
      ! cohort, it will soon be split at the splitting call of apply_disturbances().  This !
      ! new cohort will always be the last here, they will be sorted afterwards too.       !
      !------------------------------------------------------------------------------------!
      cpatch%pft(nc)    = pft
      cpatch%nplant(nc) = density
      select case (ibigleaf)
      case (0)
         !---------------------------------------------------------------------------------!
         !    SAS approximation, assign height and use it to find DBH and the structural   !
         ! (dead) biomass.                                                                 !
         !---------------------------------------------------------------------------------!
         cpatch%hite  (nc) = hgt_min(cpatch%pft(nc)) * min(1.0,height_factor)
         cpatch%dbh   (nc) = h2dbh(cpatch%hite(nc),cpatch%pft(nc))
         bdead             = size2bd(cpatch%dbh(nc),cpatch%hite(nc),cpatch%pft(nc))
         cpatch%bdeada(nc) =        agf_bs(cpatch%pft(nc))   * bdead
         cpatch%bdeadb(nc) = (1.0 - agf_bs(cpatch%pft(nc)) ) * bdead
         !---------------------------------------------------------------------------------!

      case (1)
         !---------------------------------------------------------------------------------!
         !    Big leaf approximation, assign the typical DBH and height and use them to    !
         ! find height and the structural (dead) biomass.                                  !
         !---------------------------------------------------------------------------------!
         cpatch%hite  (nc) = hgt_max(cpatch%pft(nc))
         cpatch%dbh   (nc) = dbh_bigleaf(cpatch%pft(nc))
         bdead             = size2bd(cpatch%dbh(nc),cpatch%hite(nc),cpatch%pft(nc))
         cpatch%bdeada(nc) =        agf_bs(cpatch%pft(nc))   * bdead
         cpatch%bdeadb(nc) = (1.0 - agf_bs(cpatch%pft(nc)) ) * bdead
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!




      !----- Initialise other cohort-level variables. -------------------------------------!
      call init_ed_cohort_vars(cpatch, nc, lsl,mzg,ntext_soil)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Initialise the active and storage biomass scaled by the leaf drought          !
      ! phenology (or start with 1.0 if the plant doesn't shed their leaves due to water   !
      ! stress.                                                                            !
      !------------------------------------------------------------------------------------!
      call pheninit_balive_bstorage(mzg,cpatch%pft(nc),cpatch%krdepth(nc),cpatch%hite(nc)  &
                                   ,cpatch%dbh(nc),cpatch%sla(nc),csite%soil_water(:,np)   &
                                   ,ntext_soil,cpatch%paw_avg(nc),cpatch%elongf(nc)        &
                                   ,cpatch%phenology_status(nc),cpatch%bleaf(nc)           &
                                   ,cpatch%broot(nc),cpatch%bsapwooda(nc)                  &
                                   ,cpatch%bsapwoodb(nc),cpatch%bbarka(nc)                 &
                                   ,cpatch%bbarkb(nc),cpatch%bstorage(nc),cpatch%cb(:,nc)  &
                                   ,cpatch%cb_lightmax(:,nc),cpatch%cb_moistmax(:,nc)      &
                                   ,cpatch%cb_mlmax(:,nc),cpatch%cbr_bar(nc))
      !------------------------------------------------------------------------------------!



      !----- Update balive. ---------------------------------------------------------------!
      cpatch%balive(nc) = ed_balive(cpatch, nc)
      !------------------------------------------------------------------------------------!



      !----- Compute all area indices needed. ---------------------------------------------!
      call area_indices(cpatch, nc)
      !------------------------------------------------------------------------------------!


      !----- Find the new basal area and above-ground biomass. ----------------------------!
      cpatch%basarea      (nc) = pio4 * cpatch%dbh(nc) * cpatch%dbh(nc)
      cpatch%agb          (nc) = ed_biomass(cpatch, nc)
      cpatch%btimber      (nc) = size2bt(cpatch%dbh(nc),cpatch%hite(nc),cpatch%bdeada(nc)  &
                                        ,cpatch%bsapwooda(nc),cpatch%bbarka(nc)            &
                                        ,cpatch%pft(nc))
      cpatch%thbark       (nc) = size2xb(cpatch%dbh(nc),cpatch%hite(nc),cpatch%bbarka(nc)  &
                                        ,cpatch%bbarkb(nc),cpatch%sla(nc),cpatch%pft(nc))
      cpatch%leaf_temp    (nc) = csite%can_temp  (np)
      cpatch%leaf_temp_pv (nc) = csite%can_temp  (np)
      cpatch%leaf_water   (nc) = 0.0
      cpatch%leaf_vpdef   (nc) = csite%can_vpdef (np)
      cpatch%wood_temp    (nc) = csite%can_temp  (np)
      cpatch%wood_temp_pv (nc) = csite%can_temp  (np)
      cpatch%wood_water   (nc) = 0.0
      !------------------------------------------------------------------------------------!



      !----- Make leaf and wood leaf fractions consistent. --------------------------------!
      if (cpatch%leaf_temp(nc) == t3ple) then
         cpatch%leaf_fliq(nc) = 0.5
      else if (cpatch%leaf_temp(nc) > t3ple) then
         cpatch%leaf_fliq(nc) = 1.0
      else
         cpatch%leaf_fliq(nc) = 0.0
      end if
      if (cpatch%wood_temp(nc) == t3ple) then
         cpatch%wood_fliq(nc) = 0.5
      else if (cpatch%wood_temp(nc) > t3ple) then
         cpatch%wood_fliq(nc) = 1.0
      else
         cpatch%wood_fliq(nc) = 0.0
      end if
      !------------------------------------------------------------------------------------!



      !----- Because we assigned no water, the internal energy is simply hcap*T. ----------!
      call calc_veg_hcap(cpatch%bleaf (nc),cpatch%bdeada(nc),cpatch%bsapwooda(nc)          &
                        ,cpatch%bbarka(nc),cpatch%nplant(nc),cpatch%pft(nc)                &
                        ,cpatch%leaf_hcap(nc),cpatch%wood_hcap(nc))

      cpatch%leaf_energy(nc) = cmtl2uext( cpatch%leaf_hcap     (nc)                        &
                                        , cpatch%leaf_water    (nc)                        &
                                        + cpatch%leaf_water_im2(nc)                        &
                                        , cpatch%leaf_temp     (nc)                        &
                                        , cpatch%leaf_fliq     (nc) )
      cpatch%wood_energy(nc) = cmtl2uext( cpatch%wood_hcap     (nc)                        &
                                        , cpatch%wood_water    (nc)                        &
                                        + cpatch%wood_water_im2(nc)                        &
                                        , cpatch%wood_temp     (nc)                        &
                                        , cpatch%wood_fliq     (nc) )
      call is_resolvable(csite,np,nc,.true.,.false.,'plant_patch')
      !------------------------------------------------------------------------------------!

      !----- Should plantations be considered recruits? -----------------------------------!
      cpatch%new_recruit_flag(nc) = 1
      !------------------------------------------------------------------------------------!

      !----- Sort the cohorts so that the new cohort is at the correct height bin. --------!
      call sort_cohorts(cpatch)
      !------------------------------------------------------------------------------------!

      return
   end subroutine plant_patch
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !      This routine is used to reduce the lianas height (and derived props)             !
   ! in accordance with the current tree population. Each cohort height is rescaled by     !
   ! a factor h_pruning_factor. This way lianas that where hgt_max meters tall will now    !
   ! be as tall as the tallest tree cohort. The other liana cohorts will be scaled         !
   ! proportionally.                                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine prune_lianas(csite, np, lsl)
      use stable_cohorts,  only : is_resolvable            ! ! structure
      use ed_state_vars,   only : patchtype                & ! structure
                                , sitetype                 ! ! structure
      use ed_therm_lib,    only : calc_veg_hcap            & ! function
                                , update_veg_energy_cweh   ! ! function
      use allometry,       only : area_indices             & ! subroutine
                                , ed_biomass               & ! function
                                , h2dbh                    & ! function
                                , size2bl                  & ! function
                                , size2bd                  & ! function
                                , size2krdepth             ! ! function
      use pft_coms,        only : qsw                      & ! intent(in)
                                , qbark                    & ! intent(in)
                                , agf_bs                   & ! intent(in)
                                , hgt_max                  & ! intent(in)
                                , l2n_stem                 & ! intent(in)
                                , c2n_stem                 & ! intent(in)
                                , c2n_leaf                 & ! intent(in)
                                , is_liana                 & ! intent(in)
                                , f_labile_leaf            & ! intent(in)
                                , f_labile_stem            ! ! intent(in)
      use ed_max_dims,     only : n_pft                    ! ! intent(in)
      use consts_coms,     only : pio4                     ! ! intent(in)
      use fuse_fiss_utils, only : sort_cohorts             ! ! sub-routine
      use plant_hydro,     only : rwc2tw                   & ! sub-routine
                                , twi2twe                  ! ! sub-routine

      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                  , target     :: csite
      integer                         , intent(in) :: np
      integer                         , intent(in) :: lsl
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)                 , pointer    :: cpatch
      integer                                      :: ico
      integer                                      :: ipft
      real                                         :: maxh
      real                                         :: h_pruning_factor !< height rescale factor
      real                                         :: a_struct_cohort
      real                                         :: a_fast_litter
      real                                         :: a_struct_litter
      real                                         :: a_struct_lignin
      real                                         :: a_fast_litter_n
      real                                         :: bleaf_in
      real                                         :: bsapa_in
      real                                         :: bdeada_in
      real                                         :: bbarka_in
      real                                         :: old_leaf_hcap
      real                                         :: old_wood_hcap
      real                                         :: old_leaf_water
      real                                         :: old_wood_water
      real                                         :: old_leaf_water_im2
      real                                         :: old_wood_water_im2
      real                                         :: bleaf_max
      real                                         :: delta_blfrt
      real                                         :: delta_bwood
      !------------------------------------------------------------------------------------!


      maxh = 0.0

      cpatch => csite%patch(np)

      cohortloop: do ico=1,cpatch%ncohorts

         !----- Alias for current PFT. ----------------------------------------------------!
         ipft = cpatch%pft(ico)
         !---------------------------------------------------------------------------------!

         !---------- Loop over cohorts to find the maximum height for trees ---------------!
         if (cpatch%hite(ico) > maxh .and. .not. is_liana(ipft)) then
            maxh = cpatch%hite(ico)
         end if

      end do cohortloop

      !------------- pruning_factor, how much should I reduce the height ------------------!
      h_pruning_factor    = maxh / hgt_max(17)

      !-------------------- Initialise the non-scaled litter pools. -----------------------!
      a_fast_litter   = 0.0
      a_struct_litter = 0.0
      a_struct_lignin = 0.0
      a_fast_litter_n = 0.0

      cohortloop2: do ico=1,cpatch%ncohorts

         !-------------------------- Alias for current PFT. -------------------------------!
         ipft = cpatch%pft(ico)
         !---------------------------------------------------------------------------------!

         ! Attention: if maxh turns out to be less than 1 m there's gonna be a problem
         ! because cpatch%hite will be increased instead of reduced
         if (is_liana(ipft) .and. cpatch%hite(ico) > maxh .and. maxh >= 1.0) then

            bleaf_in           = cpatch%bleaf         (ico)
            bsapa_in           = cpatch%bsapwooda     (ico)
            bdeada_in          = cpatch%bdeada        (ico)
            bbarka_in          = cpatch%bbarka        (ico)
            old_leaf_hcap      = cpatch%leaf_hcap     (ico)
            old_wood_hcap      = cpatch%wood_hcap     (ico)
            old_leaf_water     = cpatch%leaf_water    (ico)
            old_wood_water     = cpatch%wood_water    (ico)
            old_leaf_water_im2 = cpatch%leaf_water_im2(ico)
            old_wood_water_im2 = cpatch%wood_water_im2(ico)
            !add the agb_f to bdead
            !if new root depth is smaller keep the old one keep track of the value

            ! Lianas of 35m will be reduced to maxh, all
            cpatch%hite(ico)      = max(cpatch%hite(ico) * h_pruning_factor, 1.0)
            cpatch%dbh(ico)       = h2dbh (cpatch%hite(ico), ipft)
            bleaf_max             = size2bl(cpatch%dbh(ico), cpatch%hite(ico)              &
                                           ,cpatch%sla(ico), ipft)
            cpatch%bleaf(ico)     = bleaf_max * cpatch%elongf(ico)
            cpatch%bdeada(ico)    = agf_bs(ipft)                                           &
                                  * size2bd(cpatch%dbh(ico), cpatch%hite(ico), ipft)
            cpatch%bsapwooda(ico) = agf_bs(ipft) * bleaf_max * qsw(ipft) * cpatch%hite(ico)
            !  (MLO) Manfredo: although qbark is set to zero, check whether these changes
            !  are consistent with your rationale.
            cpatch%bbarka   (ico) = agf_bs(ipft) * bleaf_max * qbark(ipft)                 &
                                  * cpatch%hite(ico)


            !----- Updating LAI, WAI, and CAI. --------------------------------------------!
            call area_indices(cpatch, ico)
            !------------------------------------------------------------------------------!

            !----- Finding the new basal area and above-ground biomass. -------------------!
            cpatch%basarea(ico) = pio4 * cpatch%dbh(ico) * cpatch%dbh(ico)
            cpatch%agb(ico)     = ed_biomass(cpatch, ico)

            !----- Update rooting depth ---------------------------------------------------!
            cpatch%krdepth(ico) = size2krdepth(cpatch%hite(ico),cpatch%dbh(ico),ipft,lsl)
            !if new root depth is smaller keep the old one

            !------------------------------------------------------------------------------!
            !     It is likely that biomass has changed, therefore, update                 !
            ! vegetation energy and heat capacity.                                         !
            !------------------------------------------------------------------------------!
            call calc_veg_hcap(cpatch%bleaf(ico),cpatch%bdeada(ico),cpatch%bsapwooda(ico)  &
                              ,cpatch%bbarka(ico),cpatch%nplant(ico), cpatch%pft(ico)      &
                              ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))
            call rwc2tw(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico)                          &
                       ,cpatch%bleaf(ico),cpatch%bsapwooda(ico),cpatch%bsapwoodb(ico)      &
                       ,cpatch%bdeada(ico),cpatch%bdeadb(ico),cpatch%broot(ico)            &
                       ,cpatch%dbh(ico),cpatch%pft(ico),cpatch%leaf_water_int(ico)         &
                       ,cpatch%wood_water_int(ico))
            call twi2twe(cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico)             &
                        ,cpatch%nplant(ico),cpatch%leaf_water_im2(ico)                     &
                        ,cpatch%wood_water_im2(ico))
            call update_veg_energy_cweh(csite,np,ico,old_leaf_hcap,old_wood_hcap           &
                                       ,old_leaf_water,old_wood_water,old_leaf_water_im2   &
                                       ,old_wood_water_im2,.true.,.false.)
            !----- Update the stability status. -------------------------------------------!
            call is_resolvable(csite,np,ico,.false.,.false.,'prune_lianas')
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Compute the amount of carbon and nitrogen lost due to pruning and send   !
            ! to the necromass pools.                                                      !
            !------------------------------------------------------------------------------!
            delta_blfrt = f_labile_leaf(ipft)                                              &
                        * ( bleaf_in - cpatch%bleaf(ico) )                                 &
                        + f_labile_stem(ipft)                                              &
                        * ( bsapa_in  - cpatch%bsapwooda(ico)                              &
                          + bbarka_in - cpatch%bbarka   (ico)                              &
                          + bdeada_in - cpatch%bdeada   (ico) )
            delta_bwood = ( 1.0 - f_labile_leaf(ipft) )                                    &
                        * ( bleaf_in - cpatch%bleaf(ico) )                                 &
                        + ( 1.0 - f_labile_stem(ipft) )                                    &
                        * ( bsapa_in  - cpatch%bsapwooda(ico)                              &
                          + bbarka_in - cpatch%bbarka   (ico)                              &
                          + bdeada_in - cpatch%bdeada   (ico) )

            a_fast_litter   = a_fast_litter + delta_blfrt * cpatch%nplant(ico)
            a_fast_litter_n = a_fast_litter_n                                              &
                            + delta_blfrt / c2n_leaf(ipft) * cpatch%nplant(ico)            &
                            + delta_bwood / c2n_stem(ipft) * cpatch%nplant(ico)

            a_struct_cohort = delta_bwood * cpatch%nplant(ico)
            a_struct_litter = a_struct_litter + a_struct_cohort
            a_struct_lignin = a_struct_lignin + a_struct_cohort * l2n_stem / c2n_stem(ipft)
            !------------------------------------------------------------------------------!

         end if

      end do cohortloop2
      !------------------------------------------------------------------------------------!


      !--- Sort the cohorts so that the new cohort is at the correct height bin. ----------!
      call sort_cohorts(cpatch)
      !------------------------------------------------------------------------------------!

      !----- Load disturbance litter directly into carbon and N pools. --------------------!
      csite%fast_grnd_C(np)       = csite%fast_grnd_C(np)       + a_fast_litter
      csite%structural_grnd_C(np) = csite%structural_grnd_C(np) + a_struct_litter
      csite%structural_grnd_L(np) = csite%structural_grnd_L(np) + a_struct_lignin
      csite%fast_grnd_N(np)       = csite%fast_grnd_N(np)       + a_fast_litter_n
      !------------------------------------------------------------------------------------!

      return

   end subroutine prune_lianas

end module disturbance
!==========================================================================================!
!==========================================================================================!
