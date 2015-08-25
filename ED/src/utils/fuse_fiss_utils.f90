!==========================================================================================!
!==========================================================================================!
!  Module fuse_fiss_utils.                                                                 !
!  Routines to terminate and fuse patches, and to terminate, fuse and split cohorts.       !
!------------------------------------------------------------------------------------------!
module fuse_fiss_utils

   use ed_state_vars,only :   copy_patchtype        & ! subroutine
                            , deallocate_patchtype  & ! subroutine
                            , allocate_patchtype    & ! subroutine
                            , allocate_sitetype     & ! subroutine
                            , deallocate_sitetype   & ! subroutine
                            , copy_sitetype_mask    & ! subroutine
                            , copy_sitetype         & ! subroutine
                            , copy_patchtype_mask   ! ! subroutine

   contains
   !=======================================================================================!
   !=======================================================================================!
   !      This subroutine will sort the cohorts by size (1st = tallest, last = shortest.)  !
   ! In case there is a tie (for example, when 2 cohorts have reached the maximum possible !
   ! height, then we use DBH for tie breaking, and if they have the exact same DBH, then   !
   ! we simply pick the lowest index (as they are exactly the same).  This could cause     !
   ! some problems when the new grass allometry is implemented, though.                    !
   !---------------------------------------------------------------------------------------!
   subroutine sort_cohorts(cpatch)

      use ed_state_vars,only :  patchtype   ! ! Structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype)              , target      :: cpatch    ! Current patch
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)              , pointer     :: temppatch ! Scratch patch structure
      integer                                    :: ico       ! Counters
      integer                                    :: tallco    ! Index of tallest cohort
      real                                       :: tophgt    ! Maximum height considered
      logical                                    :: sorted    ! Patch is already sorted
      logical        , dimension(:), allocatable :: attop     ! Top cohorts, for tie-break
      !------------------------------------------------------------------------------------!
      
      !----- No need to sort an empty patch or a patch with a single cohort. --------------!
      if (cpatch%ncohorts < 2) return


      !------------------------------------------------------------------------------------!
      !     Check whether this patch is already sorted.   We don't want to do the entire   !
      ! deallocating/copying/allocating thing if it's not needed as this takes up too much !
      ! time.                                                                              !
      !------------------------------------------------------------------------------------!
      sorted = .true.
      sortcheck: do ico=1,cpatch%ncohorts-1
         sorted = cpatch%hite(ico) >= cpatch%dbh(ico+1) .and.                              &
                  cpatch%dbh(ico)  >= cpatch%dbh(ico+1)
         if (.not. sorted) exit sortcheck
      end do sortcheck
      if (sorted) return
      !------------------------------------------------------------------------------------!



      !----- Assign a scratch patch. ------------------------------------------------------!
      nullify(temppatch)
      allocate(temppatch)
      call allocate_patchtype(temppatch,cpatch%ncohorts)

      !----- Allocate the logical flag for tie-breaking. ----------------------------------!
      allocate(attop(cpatch%ncohorts))

      ico = 0
      !---- Loop until all cohorts were sorted. -------------------------------------------!
      do while(ico < cpatch%ncohorts)
         ico = ico + 1

         !----- Find the maximum height. --------------------------------------------------!
         tophgt = maxval(cpatch%hite)

         !----- Find all cohorts that are at this height. ---------------------------------!
         attop  = cpatch%hite == tophgt

         !----- Find the fattest cohort at a given height. --------------------------------!
         tallco = maxloc(cpatch%dbh,dim=1,mask=attop)

         !----- Copy to the scratch structure. --------------------------------------------!
         call copy_patchtype(cpatch,temppatch,tallco,tallco,ico,ico)

         !----- Put a non-sense DBH so this will never "win" again. -----------------------!
         cpatch%hite(tallco) = -huge(1.)
         cpatch%dbh (tallco) = -huge(1.)

      end do

      !------ Copy the scratch patch to the regular one and deallocate it. ----------------!
      call copy_patchtype(temppatch,cpatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
      call deallocate_patchtype(temppatch)
      deallocate(temppatch)

      !----- De-allocate the logical flag. ------------------------------------------------!
      deallocate(attop)

      return

   end subroutine sort_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will eliminate cohorts based on their sizes. This is intended to   !
   ! eliminate cohorts that have little contribution and thus we can speed up the run.     !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_cohorts(csite,ipa,elim_nplant,elim_lai)
      use pft_coms           , only : min_cohort_size  & ! intent(in)
                                    , l2n_stem         & ! intent(in)
                                    , c2n_stem         & ! intent(in)
                                    , c2n_storage      & ! intent(in), lookup table
                                    , c2n_leaf         ! ! intent(in), lookup table

      use decomp_coms        , only : f_labile         ! ! intent(in), lookup table

      use ed_state_vars      , only : patchtype        & ! structure
                                    , sitetype         ! ! structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)       , target      :: csite        ! Current site
      integer              , intent(in)  :: ipa          ! Current patch ID
      real                 , intent(out) :: elim_nplant  ! Nplants eliminated here
      real                 , intent(out) :: elim_lai     ! LAI eliminated here
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)      , pointer     :: cpatch       ! Current patch
      type(patchtype)      , pointer     :: temppatch    ! Scratch patch structure
      logical, dimension(:), allocatable :: remain_table ! Flag: this cohort will remain.
      integer                            :: ico          ! Counter
      integer                            :: ipft         ! PFT size
      real                               :: csize        ! Size of current cohort
      !------------------------------------------------------------------------------------!
      
      cpatch => csite%patch(ipa)
      elim_nplant = 0.
      elim_lai    = 0.

      !----- Initialize the temporary patch structures and the remain/terminate table -----!
      nullify(temppatch)
      allocate(temppatch)
      allocate(remain_table(cpatch%ncohorts))
      remain_table(:) = .true.
     
      !----- Main loop --------------------------------------------------------------------!
      do ico = 1,cpatch%ncohorts

         !----- Save the PFT type in a convenient alias. ----------------------------------!
         ipft = cpatch%pft(ico)

         !----- Checking whether the cohort size is too small -----------------------------!
         csize = cpatch%nplant(ico)                                                        &
               * (cpatch%balive(ico) + cpatch%bdead(ico) + cpatch%bstorage(ico))

         if ( csize < min_cohort_size(ipft) ) then
            !----- Cohort is indeed too small, it won't remain ----------------------------!
            remain_table(ico) = .false.
            elim_nplant = elim_nplant + cpatch%nplant(ico) * csite%area(ipa)
            elim_lai    = elim_lai    + cpatch%lai(ico)    * csite%area(ipa)

            !----- Update litter pools ----------------------------------------------------!
            csite%fsc_in(ipa) = csite%fsc_in(ipa) + cpatch%nplant(ico)                     &
                              * (f_labile(ipft)*cpatch%balive(ico) + cpatch%bstorage(ico))

            csite%fsn_in(ipa) = csite%fsn_in(ipa) + cpatch%nplant(ico)                     &
                              * ( f_labile(ipft) * cpatch%balive(ico) / c2n_leaf(ipft)     &
                                + cpatch%bstorage(ico) / c2n_storage)
            
            csite%ssc_in(ipa) = csite%ssc_in(ipa) + cpatch%nplant(ico)                     &
                              * ( (1.0 - f_labile(ipft)) * cpatch%balive(ico)              &
                                + cpatch%bdead(ico))
            
            csite%ssl_in(ipa) = csite%ssl_in(ipa) + cpatch%nplant(ico)                     &
                              * ( (1.0 - f_labile(ipft)) * cpatch%balive(ico)              &
                                + cpatch%bdead(ico) ) * l2n_stem/c2n_stem(ipft)

         end if
      end do

      !----- Copy the remaining cohorts to a temporary patch ------------------------------!
      call allocate_patchtype(temppatch,count(remain_table))
      call copy_patchtype_mask(cpatch,temppatch,remain_table,size(remain_table)            &
                              ,count(remain_table))

      !----- Reallocate the new patch and populate with the saved cohorts -----------------!
      call deallocate_patchtype(cpatch)
      call allocate_patchtype(cpatch,count(remain_table))
      call copy_patchtype(temppatch,cpatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
      call sort_cohorts(cpatch)
     
      !----- Deallocate the temporary patch -----------------------------------------------!     
      call deallocate_patchtype(temppatch)
      deallocate(temppatch)
      deallocate(remain_table)

      !----- Update the cohort census at the site level -----------------------------------!
      csite%cohort_count(ipa) = cpatch%ncohorts
      if (elim_lai < 0. .or. elim_nplant < 0.) then
         write (unit=*,fmt='(a,1x,es12.5)') 'TERMINATE: ELIM_LAI=',elim_lai
         write (unit=*,fmt='(a,1x,es12.5)') 'TERMINATE: ELIM_NPLANT=',elim_nplant
      end if
      return
   end subroutine terminate_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will eliminate tiny or empty patches. This is intended to          !
   ! eliminate patches that have little contribution and thus we can speed up the run.     !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_patches(csite)

      use ed_state_vars, only : polygontype        & ! Structure
                              , sitetype           & ! Structure
                              , patchtype          ! ! Structure
      use disturb_coms , only : min_patch_area     ! ! intent(in)
      use ed_misc_coms , only : iqoutput           & ! intent(in)
                              , imoutput           & ! intent(in)
                              , idoutput           ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)       , target      :: csite        ! Current site
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)       , pointer     :: tempsite     ! Scratch site
      integer                            :: ipa          ! Counters
      logical, dimension(:), allocatable :: remain_table ! Flag: this patch will remain.
      real                               :: total_area   ! Area of removed patches
      real                               :: elim_area    ! Area of removed patches
      real                               :: new_area     ! Just to make sure area is 1.
      real                               :: area_scale   ! Scaling area factor.
      !------------------------------------------------------------------------------------!

      allocate (remain_table(csite%npatches))
      remain_table(:) = .true.

      !------------------------------------------------------------------------------------!
      !     Loop through all the patches in this site and determine which of these patches !
      ! is too small in area to be valid. Remove these patches via the mask function.      !
      ! Realocate a new site with only the valid patches, and normalize their areas and    !
      ! plant densities to reflect the area loss.                                          !
      !------------------------------------------------------------------------------------!
      elim_area  = 0.0
      total_area = 0.0
      do ipa = 1,csite%npatches
         if (csite%area(ipa) < min_patch_area) then
            elim_area = elim_area + csite%area(ipa)
            remain_table(ipa) = .false.
         end if
         total_area = total_area + csite%area(ipa)
      end do

      !----- Use the mask to resize the patch vectors in the current site. ----------------!
      allocate(tempsite)
      call allocate_sitetype(tempsite,count(remain_table))
      call copy_sitetype_mask(csite,tempsite,remain_table,size(remain_table)               &
                             ,count(remain_table))
      call deallocate_sitetype(csite)
      call allocate_sitetype(csite,count(remain_table))

      remain_table(:)                   = .false.
      remain_table(1:tempsite%npatches) = .true.
      call copy_sitetype_mask(tempsite,csite,remain_table(1:tempsite%npatches)             &
                             ,count(remain_table),count(remain_table))
      call deallocate_sitetype(tempsite)
      deallocate(tempsite)
      deallocate(remain_table)

      !------------------------------------------------------------------------------------!
      !    Renormalize the total area.                                                     !
      !------------------------------------------------------------------------------------!
      new_area   = 0.
      area_scale = 1.0 / (total_area - elim_area)
      do ipa = 1,csite%npatches
         csite%area(ipa) = csite%area(ipa) * area_scale
         new_area        = new_area + csite%area(ipa)
      end do

      if (abs(new_area-1.0) > 1.e-5) then
         write (unit=*,fmt='(a,1x,es12.5)') ' + ELIM_AREA:',elim_area
         write (unit=*,fmt='(a,1x,es12.5)') ' + NEW_AREA: ',new_area
         call fatal_error('New_area should be 1 but it isn''t!!!','terminate_patches'      &
                         ,'fuse_fiss_utils.f90')
      end if 
      
      return
   end subroutine terminate_patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will rescale the area of the patches.  This is almost the same as !
   ! the terminate_patches subroutine, except that no patch is removed.                    !
   !---------------------------------------------------------------------------------------!
   subroutine rescale_patches(csite)
      use ed_state_vars, only : polygontype        & ! Structure
                              , sitetype           & ! Structure
                              , patchtype          ! ! Structure
      use disturb_coms , only : min_patch_area     ! ! intent(in)
      use ed_misc_coms , only : iqoutput           & ! intent(in)
                              , imoutput           & ! intent(in)
                              , idoutput           ! ! intent(in)
      use allometry    , only : size2bl            ! ! function
      use ed_max_dims  , only : n_dist_types       & ! intent(in)
                              , n_pft              ! ! intent(in)
      
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target      :: csite        ! Current site
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer     :: cpatch       ! Pointer to current site
      type(sitetype)         , pointer     :: tempsite     ! Scratch site
      integer                              :: ipa          ! Counter
      integer                              :: ico          ! Counter
      integer                              :: ilu          ! Counter
      integer                              :: ipft         ! Counter
      logical                              :: onlyone      ! Is this a single patch?
      logical, dimension  (:), allocatable :: remain_table ! Flag: this patch shall remain.
      real   , dimension  (:), allocatable :: old_area     ! Area before rescaling
      real   , dimension  (:), allocatable :: elim_area    ! Area of removed patches
      real                                 :: laimax       ! Max. LAI for current cohort
      real   , dimension  (:), allocatable :: lu_area      ! Area of disturbance type
      real                                 :: n_scale      ! Scaling nplant factor.
      real                                 :: site_area    ! Total site area (must be 1).
      real   , dimension  (:), allocatable :: patch_laimax ! Total LAI_max per patch
      real   , dimension  (:), allocatable :: lu_laimax    ! total bleaf per land use type
      integer, dimension  (:), allocatable :: lu_npatch    ! # of patches per land use type
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     No need to re-scale patches if there is a single patch left.                   !
      !------------------------------------------------------------------------------------!
      onlyone = csite%npatches == 1
      if (onlyone .and. csite%area(1) == 1.) return
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop through all the patches in this site and determine which of these patches !
      ! is too small in area to be valid. These will be removed and the remaining patches  !
      ! will be rescaled later to account for this change and changes in biomass           !
      !------------------------------------------------------------------------------------!
      allocate (remain_table(csite%npatches))
      remain_table(:) = .true.

      allocate (elim_area(n_dist_types))
      elim_area (:) = 0.0

      do ipa = 1,csite%npatches
         if (csite%area(ipa) < min_patch_area .and. (.not. onlyone)) then
            ilu = csite%dist_type(ipa)
            elim_area(ilu)= elim_area (ilu) + csite%area(ipa)
            remain_table(ipa) = .false.
         end if
      end do
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Remove patches that are smaller than the minimum non-negligible area.          !
      !------------------------------------------------------------------------------------!
      if ( sum(elim_area) > 0.0 ) then
         !----- Use the mask to resize the patch vectors in the current site. -------------!
         allocate(tempsite)
         call allocate_sitetype(tempsite,count(remain_table))
         call copy_sitetype_mask(csite,tempsite,remain_table,size(remain_table)            &
                                ,count(remain_table))
         call deallocate_sitetype(csite)
         call allocate_sitetype(csite,count(remain_table))

         remain_table(:)                   = .false.
         remain_table(1:tempsite%npatches) = .true.
         call copy_sitetype_mask(tempsite,csite,remain_table(1:tempsite%npatches)          &
                                ,count(remain_table),count(remain_table))
         call deallocate_sitetype(tempsite)
         deallocate(tempsite)
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Allocate a temporary array that will contain the potential leaf biomass of     !
      ! each patch.  This is done so phenology doesn't impact the area.                    !
      !------------------------------------------------------------------------------------!
      allocate (patch_laimax(csite%npatches))
      allocate (old_area    (csite%npatches))
      patch_laimax (:) = 0.0
      old_area     (:) = 0.0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Loop through land use types, total area per disturbance type will remain       !
      !  unchanged.  Area of the patches (PFTs) within a disturbance type will be rescaled !
      !  based on their area.                                                              !
      !------------------------------------------------------------------------------------!
      allocate (lu_npatch(n_dist_types))
      allocate (lu_laimax(n_dist_types))
      allocate (lu_area  (n_dist_types))
      lu_npatch (:) = 0
      lu_area   (:) = 0.0
      lu_laimax (:) = 0.0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Add the leaf biomass of all sites.                                             !
      !------------------------------------------------------------------------------------!
      site_area = 0.0
      do ipa=1,csite%npatches
         !----- This patch. ---------------------------------------------------------------!
         cpatch => csite%patch(ipa)
         !---------------------------------------------------------------------------------!

         !----- Find the disturbance type. ------------------------------------------------!
         ilu = csite%dist_type(ipa)
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !      Determine the maximum attainable leaf area index for each PFT, which will  !
         ! be used to rescale the relative area of the patch.                              !
         !---------------------------------------------------------------------------------!
         do ico = 1,cpatch%ncohorts
            ipft              = cpatch%pft(ico)
            laimax            = cpatch%nplant(ico) * cpatch%sla(ico)                       &
                              * size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
            patch_laimax(ipa) = patch_laimax(ipa) + laimax * csite%area(ipa)
         end do

         lu_area  (ilu) = lu_area  (ilu) + csite%area  (ipa)
         lu_laimax(ilu) = lu_laimax(ilu) + patch_laimax(ipa)
         lu_npatch(ilu) = lu_npatch(ilu) + 1

         site_area      = site_area + csite%area(ipa) 
      end do
      !------------------------------------------------------------------------------------!



      !------ Normalise the area of each land use type so the total is going to be one. ---!
      lu_area(:) = lu_area(:) / site_area
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Renormalize the total area.                                                     !
      !------------------------------------------------------------------------------------!
      site_area = 0.0
      do ipa = 1,csite%npatches

         !----- Find the disturbance type. ------------------------------------------------!
         ilu = csite%dist_type(ipa)
         !---------------------------------------------------------------------------------!


         !----- Find the new area, based on the fraction of biomass. ----------------------!
         old_area(ipa)   = csite%area(ipa)
         csite%area(ipa) = patch_laimax(ipa) / lu_laimax(ilu) * lu_area(ilu)
         site_area       = site_area + csite%area(ipa)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     We must change the scale for nplant (and all "extensive" cohort-level       !
         ! variables) so the site-level maximum LAI and maximum biomass stay the same.     !
         !---------------------------------------------------------------------------------!
         n_scale = old_area(ipa) / csite%area(ipa)
         !---------------------------------------------------------------------------------!

         cpatch => csite%patch(ipa)
         call update_cohort_extensive_props(cpatch,1,cpatch%ncohorts,n_scale)
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!






      !------------------------------------------------------------------------------------!
      !     Sanity check: total new area must be 1.0.                                      !
      !------------------------------------------------------------------------------------!
      if (abs(site_area-1.0) > 1.e-5) then
         write (unit=*,fmt='(a)'          ) '---------------------------------------------'
         write (unit=*,fmt='(a)'          ) ' PATCH BIOMASS:'
         write (unit=*,fmt='(a)'          ) ' '
         write (unit=*,fmt='(7(1x,a))'    ) '       PATCH','     LU_TYPE','PATCH_LAIMAX'   &
                                           ,'   LU_LAIMAX','     LU_AREA','    OLD_AREA'   &
                                           ,'    NEW_AREA'
         do ipa=1,csite%npatches
            ilu = csite%dist_type(ipa)
            write(unit=*,fmt='(2(1x,i12),5(1x,es12.5))')                                   &
                    ipa,ilu,patch_laimax(ipa),lu_laimax(ilu),lu_area(ilu),old_area(ipa)    &
                   ,csite%area(ipa)
         end do
         write (unit=*,fmt='(a,1x,es12.5)') ' SITE LAIMAX  :',sum(lu_laimax)
         write (unit=*,fmt='(a,1x,es12.5)') ' SITE_AREA    :',site_area
         write (unit=*,fmt='(a)'          ) '---------------------------------------------'
         call fatal_error('New_area should be 1 but it isn''t!!!','rescale_patches'        &
                         ,'fuse_fiss_utils.f90')
      end if 
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Free memory before we leave the sub-routine.                                   !
      !------------------------------------------------------------------------------------!
      deallocate(patch_laimax)
      deallocate(remain_table)
      deallocate(old_area    )
      deallocate(lu_area     )
      deallocate(lu_laimax   )
      deallocate(lu_npatch   )
      deallocate(elim_area   )
      !------------------------------------------------------------------------------------!
      return
   end subroutine rescale_patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will perform cohort fusion based on various similarity criteria to  !
   ! determine whether they can be fused with no significant loss of information. The user !
   ! is welcome to set up a benchmark, but should be aware that no miracles will happen    !
   ! here. If there are more very distinct cohorts than maxcohort, then the user will need !
   ! to live with that and accept life is not always fair with those with limited          !
   ! computational resources.                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_cohorts(csite,ipa, green_leaf_factor, lsl, fuse_initial)

      use ed_state_vars       , only : sitetype            & ! Structure
                                     , patchtype           ! ! Structure
      use pft_coms            , only : rho                 & ! intent(in)
                                     , b1Ht                & ! intent(in)
                                     , hgt_max             & ! intent(in)
                                     , sla                 & ! intent(in)
                                     , is_grass            & ! intent(in)
                                     , hgt_ref             ! ! intent(in)
      use fusion_fission_coms , only : fusetol_h           & ! intent(in)
                                     , fusetol             & ! intent(in)
                                     , lai_fuse_tol        & ! intent(in)
                                     , fuse_relax          & ! intent(in)
                                     , coh_tolerance_max   ! ! intent(in)
      use ed_max_dims         , only : n_pft               ! ! intent(in)
      use mem_polygons        , only : maxcohort           ! ! intent(in)
      use canopy_layer_coms   , only : crown_mod           ! ! intent(in)
      use allometry           , only : dbh2h               & ! function
                                     , size2bl             ! ! function
      use ed_misc_coms        , only : igrass              ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target      :: csite             ! Current site
      integer                , intent(in)  :: ipa               ! Current patch ID
      real, dimension(n_pft) , intent(in)  :: green_leaf_factor ! 
      integer                , intent(in)  :: lsl               ! Lowest soil level
      logical                , intent(in)  :: fuse_initial      ! Initialisation step?
      !----- Local variables --------------------------------------------------------------!
      logical, dimension(:)  , allocatable :: fuse_table     ! Flag, remaining cohorts
      type(patchtype)        , pointer     :: cpatch         ! Current patch
      type(patchtype)        , pointer     :: temppatch      ! Scratch patch
      integer                              :: donc,recc,ico3 ! Counters
      logical                              :: fusion_test    ! Flag: proceed with fusion?
      real                                 :: newn           ! new nplants of merged coh.
      real                                 :: lai_max        ! Maximum LAI the fused 
                                                             !    cohort could have.
      real                                 :: total_size     ! Total size
      real                                 :: tolerance_mult ! Multiplication factor
      integer                              :: ncohorts_old   ! # of coh. before fusion test
      real                                 :: mean_dbh       ! Mean DBH           (???)
      real                                 :: mean_hite      ! Mean height        (???)
      real                                 :: new_size       ! New size
      integer                              :: ntall          ! # of tall cohorts  (???)
      integer                              :: nshort         ! # of short cohorts (???)
      logical                              :: any_fusion     ! Flag: was there any fusion?
      !------------------------------------------------------------------------------------!


      !----- Start with no factor ---------------------------------------------------------!
      tolerance_mult = 1.0

      cpatch => csite%patch(ipa)

      !------------------------------------------------------------------------------------!
      !     Return if maxcohort is 0 (flag for no cohort fusion), or if the patch is empty !
      ! or has a single cohort.                                                            !
      !------------------------------------------------------------------------------------!
      if (maxcohort == 0 .or. cpatch%ncohorts < 2) return

      !------------------------------------------------------------------------------------!
      !    Calculate mean DBH and HITE to help with the normalization of differences mean  !
      ! hite is not being used right now, but can be optioned in the future if it seems    !
      ! advantageous.                                                                      !
      !------------------------------------------------------------------------------------!
      mean_dbh  = 0.0
      mean_hite = 0.0
      nshort    = 0
      ntall     = 0
      do ico3 = 1,cpatch%ncohorts
         !---------------------------------------------------------------------------------!
         !    Get fusion height threshold.  Height is a good predictor when plants are     !
         ! growing in height, but it approaches the maximum height DBH becomes the only    !
         ! possible predictor because height saturates.                                    !
         !---------------------------------------------------------------------------------!
         if (cpatch%hite(ico3) < (0.95 * hgt_max(cpatch%pft(ico3))) ) then
            mean_hite = mean_hite + cpatch%hite(ico3)
            nshort    = nshort + 1
         else
            mean_dbh  = mean_dbh + cpatch%dbh(ico3)
            ntall     = ntall + 1
         end if
      end do 
      !------------------------------------------------------------------------------------!
      if (ntall  > 0) mean_dbh = mean_dbh   / real(ntall)
      if (nshort > 0) mean_hite= mean_hite  / real(nshort)

      !----- Initialize table. In principle, all cohorts stay. ----------------------------!
      allocate(fuse_table(cpatch%ncohorts))
      fuse_table(:) = .true.

      force_fusion: do
         
         ncohorts_old =  count(fuse_table) ! Save current number of cohorts ---------------!
         
         donloop:do donc = 1,cpatch%ncohorts-1
            if (.not. fuse_table(donc)) cycle donloop ! This one is gone, move to next.

            recloop: do recc = donc+1,cpatch%ncohorts
               if (.not. fuse_table(recc)) cycle recloop ! This one is gone, move to next.
                                                         ! Hope it never happens...

               !---------------------------------------------------------------------------!
               !     Test for similarity.  Again, we use height to assess similarity only  !
               ! when the cohort is not approaching the maximum height.  If this is the    !
               ! case, then we use DBH to test.                                            !
               !---------------------------------------------------------------------------!
               if (cpatch%hite(donc) >= (0.95 * hgt_max(cpatch%pft(donc))) ) then
                  mean_dbh=0.5*(cpatch%dbh(donc)+cpatch%dbh(recc))
                  fusion_test = ( abs(cpatch%dbh(donc) - cpatch%dbh(recc)))/mean_dbh       &
                              < fusetol * tolerance_mult
               elseif (fuse_relax) then
                  fusion_test = ( abs(cpatch%hite(donc) - cpatch%hite(recc))               &
                                     / (0.5*(cpatch%hite(donc) + cpatch%hite(recc)))  <    &
                                fusetol * tolerance_mult)  
               else
                  fusion_test = (abs(cpatch%hite(donc) - cpatch%hite(recc))  <             &
                                fusetol_h * tolerance_mult)
               end if

               if (fusion_test) then

                  !----- New cohort has the total number of plants ------------------------!
                  newn = cpatch%nplant(donc) + cpatch%nplant(recc)

                  !------------------------------------------------------------------------!
                  !     We now check the maximum LAI the fused cohorts could have.  We     !
                  ! don't want the cohort to have a very large LAI.  If both cohorts have  !
                  ! leaves fully flushed, this is the same as adding the individual LAIs,  !
                  ! but if they are not, we need to consider that LAI may grow...          !
                  !------------------------------------------------------------------------!
                  if (is_grass(cpatch%pft(donc)) .and. igrass==1) then
                      !--use actual bleaf for grass
                      lai_max = ( cpatch%nplant(recc) * cpatch%bleaf(recc)                 &
                                + cpatch%nplant(donc) * cpatch%bleaf(donc) )               &
                                * cpatch%sla(recc)
                  else
                      !--use dbh for trees
                      lai_max = ( cpatch%nplant(recc)                                      &
                                * size2bl(cpatch%dbh(recc),cpatch%hite(recc)               &
                                         ,cpatch%pft(recc))                                &
                                + cpatch%nplant(donc)                                      &
                                * size2bl(cpatch%dbh(donc),cpatch%hite(donc)               &
                                         ,cpatch%pft(donc)))                               &
                                * cpatch%sla(recc)
                  end if

                  !----- Checking the total size of this cohort before and after fusion. --!
                  total_size = cpatch%nplant(donc) * ( cpatch%balive(donc)                 &
                                                     + cpatch%bdead(donc)                  &
                                                     + cpatch%bstorage(donc) )             &
                             + cpatch%nplant(recc) * ( cpatch%balive(recc)                 &
                                                     + cpatch%bdead(recc)                  &
                                                     + cpatch%bstorage(recc) )

                  
                  
                  
                  !------------------------------------------------------------------------!
                  !    Six conditions must be met to allow two cohorts to be fused:        !
                  ! 1. Both cohorts must have the same PFT;                                !
                  ! 2. Combined LAI won't be too large.                                    !
                  ! 3. Both cohorts must have the same status with respect to the first    !
                  !    census.                                                             !
                  ! 4. Both cohorts must have the same recruit status with respect to the  !
                  !    first census.                                                       !
                  ! 5. Both cohorts must have the same recruitment status with respect to  !
                  !    the DBH.                                                            !
                  ! 6. Both cohorts must have the same recruitment status with respect to  !
                  !    the census.                                                         !
                  ! 7. Both cohorts must have the same phenology status.                   !
                  !------------------------------------------------------------------------!
                  if (     cpatch%pft(donc)              == cpatch%pft(recc)               &
                     .and. lai_max                        < lai_fuse_tol*tolerance_mult    &
                     .and. cpatch%first_census(donc)     == cpatch%first_census(recc)      &
                     .and. cpatch%new_recruit_flag(donc) == cpatch%new_recruit_flag(recc)  &
                     .and. cpatch%recruit_dbh     (donc) == cpatch%recruit_dbh(recc)       &
                     .and. cpatch%census_status   (donc) == cpatch%census_status(recc)     &
                     .and. cpatch%phenology_status(donc) == cpatch%phenology_status(recc)  &
                     ) then

                     !----- Proceed with fusion -------------------------------------------!
                     call fuse_2_cohorts(cpatch,donc,recc,newn                             &
                                        ,green_leaf_factor(cpatch%pft(donc))               &
                                        ,csite%can_prss(ipa),csite%can_shv(ipa),lsl        &
                                        ,fuse_initial)

                     !----- Flag donating cohort as gone, so it won't be checked again. ---!
                     fuse_table(donc) = .false.
                     
                     !----- Check whether total size and LAI are conserved. ---------------!
                     new_size = cpatch%nplant(recc) * ( cpatch%balive(recc)                &
                                                      + cpatch%bdead(recc)                 &
                                                      + cpatch%bstorage(recc) )
                     if (new_size < 0.99* total_size .or. new_size > 1.01* total_size )    &
                     then
                        write (unit=*,fmt='(a,1x,es14.7)') 'OLD SIZE: ',total_size
                        write (unit=*,fmt='(a,1x,es14.7)') 'NEW SIZE: ',new_size
                        call fatal_error('Cohort fusion didn''t conserve plant size!!!'    &
                                        &,'fuse_2_cohorts','fuse_fiss_utils.f90')
                     end if
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    Recalculate the means                                            !
                     !---------------------------------------------------------------------!
                     mean_dbh  = 0.0
                     mean_hite = 0.0
                     nshort    = 0
                     ntall     = 0
                     recalcloop: do ico3 = 1,cpatch%ncohorts
                        if (.not. fuse_table(ico3)) cycle recalcloop
                        !----- Get fusion height threshold --------------------------------!
                        if (cpatch%hite(ico3) < (0.95 * hgt_max(cpatch%pft(ico3))) ) then
                           mean_hite = mean_hite + cpatch%hite(ico3)
                           nshort = nshort+1
                        else
                           mean_dbh = mean_dbh + cpatch%dbh(ico3)
                           ntall=ntall+1
                        end if
                     end do recalcloop
                     !---------------------------------------------------------------------!
                     cycle donloop
                  end if
                  !------------------------------------------------------------------------!
               end if
            end do recloop
         end do donloop

         !------ If we are under maxcohort, no need to continue fusing. -------------------!
         if ( count(fuse_table) <= abs(maxcohort)) exit force_fusion
         !------ If no fusion happened and the tolerance exceeded the maximum, I give up. -!
         if ( count(fuse_table) == ncohorts_old .and. tolerance_mult > coh_tolerance_max ) &
            exit force_fusion

         tolerance_mult = tolerance_mult * 1.01
         ncohorts_old = count(fuse_table)
      end do force_fusion

      !----- If any fusion has happened, then we need to rearrange cohorts. ---------------!
      any_fusion = .not. all(fuse_table)
      if (any_fusion) then

         !---------------------------------------------------------------------------------!
         !     Now copy the merged patch to a temporary patch using the fuse_table as a    !
         ! mask.  Then allocate a temporary patch, copy the remaining cohorts there.       !
         !---------------------------------------------------------------------------------!
         nullify (temppatch)
         allocate(temppatch)
         call allocate_patchtype(temppatch,cpatch%ncohorts)
         call copy_patchtype_mask(cpatch,temppatch,fuse_table,size(fuse_table)             &
                                 ,count(fuse_table))

         !----- Now I reallocate the current patch with its new reduced size. -------------!
         call deallocate_patchtype(cpatch)  
         call allocate_patchtype(cpatch,count(fuse_table))
  
         !----- Make fuse_table true to all remaining cohorts. ----------------------------!
         fuse_table(:)                 = .false.
         fuse_table(1:cpatch%ncohorts) = .true.
         call copy_patchtype_mask(temppatch,cpatch,fuse_table,size(fuse_table)             &
                                 ,count(fuse_table))

         !----- Discard the scratch patch. ------------------------------------------------!
         call deallocate_patchtype(temppatch)
         deallocate(temppatch)  

         !----- Sort cohorts by size again, and update the cohort census for this patch. --!
         call sort_cohorts(cpatch)
         csite%cohort_count(ipa) = count(fuse_table)
      end if

      !----- Deallocate the aux. table ----------------------------------------------------!
      deallocate(fuse_table)
     
      return
   end subroutine fuse_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will split two cohorts if its LAI has become too large.  This is    !
   ! only necessary when we solve radiation cohort by cohort rather than layer by layer.   !
   !---------------------------------------------------------------------------------------!
   subroutine split_cohorts(cpatch, green_leaf_factor, lsl)

      use ed_state_vars        , only : patchtype              & ! structure
                                      , copy_patchtype         ! ! sub-routine
      use pft_coms             , only : q                      & ! intent(in), lookup table
                                      , qsw                    & ! intent(in), lookup table
                                      , is_grass               ! ! intent(in)
      use fusion_fission_coms  , only : lai_tol                ! ! intent(in)
      use ed_max_dims          , only : n_pft                  ! ! intent(in)
      use allometry            , only : dbh2h                  & ! function
                                      , bd2dbh                 & ! function
                                      , bl2dbh                 & ! function
                                      , bl2h                   & ! function
                                      , dbh2bd                 ! ! function
      use ed_misc_coms         , only : iqoutput               & ! intent(in)
                                      , imoutput               & ! intent(in)
                                      , idoutput               & ! intent(in)
                                      , igrass                 ! ! intent(in)
      use canopy_layer_coms    , only : crown_mod              ! ! intent(in)
      implicit none
      !----- Constants --------------------------------------------------------------------!
      real                   , parameter   :: epsilon=0.0001    ! Tweak factor...
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype)        , target      :: cpatch            ! Current patch
      real, dimension(n_pft) , intent(in)  :: green_leaf_factor !
      integer                , intent(in)  :: lsl               ! Lowest soil level
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer     :: temppatch         ! Temporary patch
      logical, dimension(:)  , allocatable :: split_mask        ! Flag: split this cohort
      integer                              :: ico               ! Counter
      integer                              :: inew              ! Counter
      integer                              :: ncohorts_new      ! New # of cohorts
      integer                              :: tobesplit         ! # of cohorts to be split
      integer                              :: ipft              ! PFT type
      real                                 :: stai              ! Potential TAI
      real                                 :: old_nplant        ! Old nplant
      real                                 :: new_nplant        ! New nplant
      real                                 :: old_size          ! Old size
      real                                 :: new_size          ! New size
      !------------------------------------------------------------------------------------!


      !----- Initialize the vector with splitting table -----------------------------------!
      allocate(split_mask(cpatch%ncohorts))
      split_mask(:) = .false.
      old_nplant = 0.
      old_size   = 0.
      !----- Loop through cohorts ---------------------------------------------------------!
      do ico = 1,cpatch%ncohorts
         ipft = cpatch%pft(ico)

         !---------------------------------------------------------------------------------! 
         !     STAI is the potential TAI that this cohort has when its leaves are fully    !
         ! flushed.                                                                        !
         !---------------------------------------------------------------------------------! 
         stai = cpatch%nplant(ico) * cpatch%balive(ico) * green_leaf_factor(ipft)          &
              * q(ipft) / ( 1.0 + q(ipft) + qsw(ipft) * cpatch%hite(ico) )                 &
              * cpatch%sla(ico) + cpatch%wai(ico)

         !----- If the resulting TAI is too large, split this cohort. ---------------------!
         split_mask(ico) = stai > lai_tol
         
         old_nplant = old_nplant + cpatch%nplant(ico)
         old_size   = old_size   + cpatch%nplant(ico) * ( cpatch%balive(ico)               &
                                                        + cpatch%bdead(ico)                &
                                                        + cpatch%bstorage(ico) )
      end do

      !----- Compute the new number of cohorts. -------------------------------------------!
      tobesplit    = count(split_mask)
      ncohorts_new = cpatch%ncohorts + tobesplit
      
      if (tobesplit > 0) then

         !----- Allocate the temppatch. ---------------------------------------------------!
         nullify(temppatch)
         allocate(temppatch)
         call allocate_patchtype(temppatch,cpatch%ncohorts)

         !----- Fill the temp space with the current patches. -----------------------------!
         call copy_patchtype(cpatch,temppatch,1,cpatch%ncohorts,1,cpatch%ncohorts)

         !----- Deallocate the current patch. ---------------------------------------------!
         call deallocate_patchtype(cpatch)

         !----- Re-allocate the current patch. --------------------------------------------!
         call allocate_patchtype(cpatch,ncohorts_new)

         !----- Transfer the temp values back in. -----------------------------------------!
         call copy_patchtype(temppatch,cpatch,1,temppatch%ncohorts,1,temppatch%ncohorts)

         !----- Remove the temporary patch. -----------------------------------------------!
         call deallocate_patchtype(temppatch)
         deallocate(temppatch)
     
         inew = size(split_mask)
         do ico = 1,size(split_mask)

            if (split_mask(ico)) then

               !---------------------------------------------------------------------------!
               !   Half the densities of the original cohort.  All "extensive" variables   !
               ! must be rescaled.                                                         !
               !---------------------------------------------------------------------------!
               call update_cohort_extensive_props(cpatch,ico,ico,0.5)
               !---------------------------------------------------------------------------!


               !----- Apply these values to the new cohort. -------------------------------!
               inew = inew+1
               call copy_patchtype(cpatch,cpatch,ico,ico,inew,inew)
               !---------------------------------------------------------------------------!

               !----- Tweaking bdead, to ensure carbon is conserved. ----------------------!
               if (is_grass(cpatch%pft(ico)) .and. igrass==1) then 
                   !-- use bleaf for grass
                   cpatch%bleaf(ico)  = cpatch%bleaf(ico) * (1.-epsilon)
                   cpatch%dbh  (ico)  = bl2dbh(cpatch%bleaf(ico), cpatch%pft(ico))
                   cpatch%hite (ico)  = bl2h(cpatch%bleaf(ico), cpatch%pft(ico))

                   cpatch%bleaf(inew)  = cpatch%bleaf(inew) * (1.+epsilon)
                   cpatch%dbh  (inew)  = bl2dbh(cpatch%bleaf(inew), cpatch%pft(inew))
                   cpatch%hite (inew)  = bl2h(cpatch%bleaf(inew), cpatch%pft(inew))               
               else
                   !-- use bdead for trees
                   cpatch%bdead(ico)  = cpatch%bdead(ico) * (1.-epsilon)
                   cpatch%dbh  (ico)  = bd2dbh(cpatch%pft(ico), cpatch%bdead(ico))
                   cpatch%hite (ico)  = dbh2h(cpatch%pft(ico), cpatch%dbh(ico))

                   cpatch%bdead(inew) = cpatch%bdead(inew) * (1.+epsilon)
                   cpatch%dbh  (inew) = bd2dbh(cpatch%pft(inew), cpatch%bdead(inew))
                   cpatch%hite (inew) = dbh2h(cpatch%pft(inew), cpatch%dbh(inew))
               end if
               !---------------------------------------------------------------------------!

            end if
         end do

         !----- After splitting, cohorts may need to be sorted again... -------------------!
         call sort_cohorts(cpatch)

         !----- Checking whether the total # of plants is conserved... --------------------!
         new_nplant = 0.
         new_size   = 0.
         do ico=1,cpatch%ncohorts
            new_nplant = new_nplant + cpatch%nplant(ico)
            new_size   = new_size   + cpatch%nplant(ico) * ( cpatch%balive(ico)            &
                                                           + cpatch%bdead(ico)             &
                                                           + cpatch%bstorage(ico) )
         end do
         if (new_nplant < 0.99 * old_nplant .or. new_nplant > 1.01 * old_nplant .or.       &
             new_size   < 0.99 * old_size   .or. new_size   > 1.01 * old_size) then
            write (unit=*,fmt='(a,1x,es14.7)') 'OLD NPLANT: ',old_nplant
            write (unit=*,fmt='(a,1x,es14.7)') 'NEW NPLANT: ',new_nplant
            write (unit=*,fmt='(a,1x,es14.7)') 'OLD SIZE:   ',old_size
            write (unit=*,fmt='(a,1x,es14.7)') 'NEW SIZE:   ',new_size
            call fatal_error('Cohort splitting didn''t conserve plants!!!'                 &
                                        &,'split_cohorts','fuse_fiss_utils.f90')
         end if
         
      end if
      deallocate(split_mask)
      return
   end subroutine split_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will merge two cohorts into 1. The donating cohort (donc) is the  !
   ! one that will be deallocated, while the receptor cohort (recc) will contain the       !
   !  information from both cohorts.                                                       !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_cohorts(cpatch,donc,recc, newn,green_leaf_factor,can_prss,can_shv,lsl &
                            ,fuse_initial)
      use ed_state_vars      , only : patchtype              ! ! Structure
      use pft_coms           , only : q                      & ! intent(in), lookup table
                                    , qsw                    & ! intent(in), lookup table
                                    , is_grass               ! ! intent(in)
      use therm_lib          , only : uextcm2tl              & ! subroutine
                                    , vpdefil                & ! subroutine
                                    , qslif                  ! ! function
      use allometry          , only : dbh2krdepth            & ! function
                                    , bd2dbh                 & ! function
                                    , bl2dbh                 & ! function
                                    , bl2h                   & ! function
                                    , dbh2h                  ! ! function
      use ed_max_dims        , only : n_mort                 ! ! intent(in)
      use ed_misc_coms       , only : writing_long           & ! intent(in)
                                    , writing_eorq           & ! intent(in)
                                    , writing_dcyc           & ! intent(in)
                                    , ndcycle                & ! intent(in)
                                    , igrass                 ! ! intent(in)
      use consts_coms        , only : lnexp_min              & ! intent(in)
                                    , lnexp_max              & ! intent(in)
                                    , tiny_num               ! ! intent(in)
      use fusion_fission_coms, only : corr_cohort
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype) , target     :: cpatch            ! Current patch
      integer                      :: donc              ! Donating cohort.
      integer                      :: recc              ! Receptor cohort.
      real            , intent(in) :: newn              ! New nplant
      real            , intent(in) :: green_leaf_factor ! Green leaf factor
      real            , intent(in) :: can_prss          ! Canopy air pressure
      real            , intent(in) :: can_shv           ! Canopy air specific humidity
      integer         , intent(in) :: lsl               ! Lowest soil level
      logical         , intent(in) :: fuse_initial      ! Called from initialisation
      !----- Local variables --------------------------------------------------------------!
      integer                      :: imon              ! Month for cb loop
      integer                      :: t                 ! Time of day for dcycle loop
      integer                      :: imty              ! Mortality type
      real                         :: newni             ! Inverse of new nplants
      real                         :: exp_mort_donc     ! Exp(mortality) donor
      real                         :: exp_mort_recc     ! Exp(mortality) receptor
      real                         :: rlai              ! LAI of receiver
      real                         :: dlai              ! LAI of donor
      real                         :: rwai              ! WAI of receiver
      real                         :: dwai              ! WAI of donor
      real                         :: rnplant           ! nplant of receiver
      real                         :: dnplant           ! nplant of donor
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the scaling factor for variables that are not "extensive".                 !
      !  - If the unit is X/plant, then we scale by nplant.                                !
      !  - If the unit is X/m2_leaf, then we scale by LAI.                                 !
      !  - If the unit is X/m2_wood, then we scale by WAI.                                 !
      !  - If the unit is X/m2_gnd, then we add, since they are "extensive".               !
      !------------------------------------------------------------------------------------!
      newni   = 1.0 / newn
      rnplant = cpatch%nplant(recc) * newni
      dnplant = 1.0 - rnplant
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    This is a fix for when two cohorts with very very low LAI are fused             !
      ! it prevents numerical errors (RGK 8-18-2014).                                      !
      ! (MLO 11-24-2014): turned rlai and dlai relative weights, so it works in all cases. !
      ! Also, applied the same idea to WAI-dependent variables.                            !
      !------------------------------------------------------------------------------------!
      if (cpatch%lai(recc) + cpatch%lai(donc) > 0 ) then
         rlai    = cpatch%lai(recc) / ( cpatch%lai(recc) + cpatch%lai(donc) )
         dlai    = 1.0 - rlai
      else
         rlai    = 0.5
         dlai    = 0.5
      end if
      if (cpatch%wai(recc) + cpatch%wai(donc) > 0 ) then
         rwai    = cpatch%wai(recc) / ( cpatch%wai(recc) + cpatch%wai(donc) )
         dwai    = 1.0 - rwai
      else
         rwai    = 0.5
         dwai    = 0.5
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Find DBH and height.  Make sure that carbon is conserved.                     !
      !------------------------------------------------------------------------------------!
      if (is_grass(cpatch%pft(donc)) .and. igrass == 1) then
          !----- New grass scheme, use bleaf then find DBH and height. --------------------!
          cpatch%bleaf(recc) = cpatch%bleaf(recc) * rnplant + cpatch%bleaf(donc) * dnplant
          cpatch%dbh(recc)   = bl2dbh(cpatch%bleaf(recc), cpatch%pft(recc))
          cpatch%hite(recc)  = bl2h  (cpatch%bleaf(recc), cpatch%pft(recc))
          !--------------------------------------------------------------------------------!
      else
          !----- Trees, or old grass scheme.  Use bdead then find DBH and height. ---------!
          cpatch%bdead(recc) = cpatch%bdead(recc) * rnplant + cpatch%bdead(donc) * dnplant
          cpatch%dbh(recc)   = bd2dbh(cpatch%pft(recc), cpatch%bdead(recc))
          cpatch%hite(recc)  = dbh2h(cpatch%pft(recc),  cpatch%dbh(recc))
          !--------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !----- Rooting depth. ---------------------------------------------------------------!
      cpatch%krdepth(recc) = dbh2krdepth(cpatch%hite(recc),cpatch%dbh(recc)                &
                                        ,cpatch%pft(recc),lsl)
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Conserving carbon to get balive, bleaf, and bstorage.                          !
      !------------------------------------------------------------------------------------!
      cpatch%balive           (recc) = cpatch%balive          (recc) * rnplant             &
                                     + cpatch%balive          (donc) * dnplant
      cpatch%broot            (recc) = cpatch%broot           (recc) * rnplant             &
                                     + cpatch%broot           (donc) * dnplant
      cpatch%bsapwooda        (recc) = cpatch%bsapwooda       (recc) * rnplant             &
                                     + cpatch%bsapwooda       (donc) * dnplant
      cpatch%bsapwoodb        (recc) = cpatch%bsapwoodb       (recc) * rnplant             &
                                     + cpatch%bsapwoodb       (donc) * dnplant
      cpatch%bstorage         (recc) = cpatch%bstorage        (recc) * rnplant             &
                                     + cpatch%bstorage        (donc) * dnplant
      cpatch%bseeds           (recc) = cpatch%bseeds          (recc) * rnplant             &
                                     + cpatch%bseeds          (donc) * dnplant
      cpatch%leaf_maintenance (recc) = cpatch%leaf_maintenance(recc) * rnplant             &
                                     + cpatch%leaf_maintenance(donc) * dnplant
      cpatch%root_maintenance (recc) = cpatch%root_maintenance(recc) * rnplant             &
                                     + cpatch%root_maintenance(donc) * dnplant
      cpatch%leaf_drop        (recc) = cpatch%leaf_drop       (recc) * rnplant             &
                                     + cpatch%leaf_drop       (donc) * dnplant
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Bleaf must be zero if phenology status is 2.  This is probably done correctly   !
      ! throughout the code, but being safe here.                                          !
      !------------------------------------------------------------------------------------!
      if (cpatch%phenology_status(recc) /= -2) then
         cpatch%bleaf(recc)  = rnplant * cpatch%bleaf(recc) + dnplant * cpatch%bleaf(donc)
      else
         cpatch%bleaf(recc)  = 0.
      end if
      !------------------------------------------------------------------------------------!



      !------ Energy, water, and heat capacity are extensive, add them. -------------------!
      cpatch%leaf_energy(recc) = cpatch%leaf_energy(recc) + cpatch%leaf_energy(donc)
      cpatch%leaf_water (recc) = cpatch%leaf_water (recc) + cpatch%leaf_water (donc)
      cpatch%leaf_hcap  (recc) = cpatch%leaf_hcap  (recc) + cpatch%leaf_hcap  (donc)
      cpatch%wood_energy(recc) = cpatch%wood_energy(recc) + cpatch%wood_energy(donc)
      cpatch%wood_water (recc) = cpatch%wood_water (recc) + cpatch%wood_water (donc)
      cpatch%wood_hcap  (recc) = cpatch%wood_hcap  (recc) + cpatch%wood_hcap  (donc)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      We update temperature and liquid water fraction.  We check whether the heat   !
      ! capacity is non-zero.  If it is a normal number, use the standard thermodynamic    !
      ! library, otherwise, average temperature, this is probably a blend of tiny cohorts  !
      ! that couldn't be solved, or the wood is not solved.                                !
      !------------------------------------------------------------------------------------!
      if ( cpatch%leaf_hcap(recc) > 0. ) then
         !----- Update temperature using the standard thermodynamics. ---------------------!
         call uextcm2tl(cpatch%leaf_energy(recc),cpatch%leaf_water(recc)                   &
                       ,cpatch%leaf_hcap(recc),cpatch%leaf_temp(recc)                      &
                       ,cpatch%leaf_fliq(recc))
         
         
      else 
         !----- Leaf temperature cannot be found using uextcm2tl, this is a singularity. --!
         cpatch%leaf_temp(recc)  = cpatch%leaf_temp(recc) * rnplant                        &
                                 + cpatch%leaf_temp(donc) * dnplant
         cpatch%leaf_fliq(recc)  = 0.0
      end if


      if ( cpatch%wood_hcap(recc) > 0. ) then
         !----- Update temperature using the standard thermodynamics. ---------------------!
         call uextcm2tl(cpatch%wood_energy(recc),cpatch%wood_water(recc)                   &
                       ,cpatch%wood_hcap(recc),cpatch%wood_temp(recc)                      &
                       ,cpatch%wood_fliq(recc))
      else 
         !----- Wood temperature cannot be found using uextcm2tl, this is a singularity. --!
         cpatch%wood_temp(recc)  = cpatch%wood_temp(recc) * rnplant                        &
                                 + cpatch%wood_temp(donc)  * dnplant
         cpatch%wood_fliq(recc)  = 0.0
      end if
      
      !----- Set time-steps temperatures as the current. ----------------------------------!
      cpatch%leaf_temp_pv(recc) = cpatch%leaf_temp(recc)
      cpatch%wood_temp_pv(recc) = cpatch%wood_temp(recc)
      !------------------------------------------------------------------------------------!

      !------ Find the intercellular value assuming saturation. ---------------------------!
      cpatch%lint_shv(recc) = qslif(can_prss,cpatch%leaf_temp(recc))
      !------------------------------------------------------------------------------------!

      !------ Find the vapour pressure deficit. -------------------------------------------!
      cpatch%leaf_vpdef(recc) = vpdefil(can_prss,cpatch%leaf_temp(recc),can_shv,.true.)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     CB and CB_Xmax are scaled by population, as they are in kgC/plant/yr.          !
      !------------------------------------------------------------------------------------!
      ! RK: I think the below comment is no longer true. Per gh-24 reverting again to      !
      ! calculate CBR from running means of CB and CB_Xmax. 
      ! "The relative carbon balance, however, is no longer derived from the annual values !
      ! of CB, CB_LightMax, and CB_MoistMax, but tracked independently as it used to be    !
      ! done in ED-1.0."                                                                   !
      !------------------------------------------------------------------------------------!
      do imon = 1,13
         cpatch%cb         (imon,recc) = cpatch%cb          (imon,recc) * rnplant          &
                                       + cpatch%cb          (imon,donc) * dnplant
         cpatch%cb_lightmax(imon,recc) = cpatch%cb_lightmax (imon,recc) * rnplant          &
                                       + cpatch%cb_lightmax (imon,donc) * dnplant
         cpatch%cb_moistmax(imon,recc) = cpatch%cb_moistmax (imon,recc) * rnplant          &
                                       + cpatch%cb_moistmax (imon,donc) * dnplant
         cpatch%cb_mlmax   (imon,recc) = cpatch%cb_mlmax    (imon,recc) * rnplant          &
                                       + cpatch%cb_mlmax    (imon,donc) * dnplant
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Relative carbon balance is also averaged between the cohorts, to avoid wild    !
      ! oscillations in mortality when cohorts are fused.  This is the original method     !
      ! used in ED-1.0.                                                                    !
      !------------------------------------------------------------------------------------!
      cpatch%cbr_bar(recc) = cpatch%cbr_bar(recc) * rnplant + cpatch%cbr_bar(donc) * dnplant
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Today variables are extensive.                                                  !
      !------------------------------------------------------------------------------------!
      cpatch%today_gpp          (recc) = cpatch%today_gpp          (recc)                  &
                                       + cpatch%today_gpp          (donc)
                                  
      cpatch%today_nppleaf      (recc) = cpatch%today_nppleaf      (recc)                  &
                                       + cpatch%today_nppleaf      (donc)
                                  
      cpatch%today_nppfroot     (recc) = cpatch%today_nppfroot     (recc)                  &
                                       + cpatch%today_nppfroot     (donc)
                                  
      cpatch%today_nppsapwood   (recc) = cpatch%today_nppsapwood   (recc)                  &
                                       + cpatch%today_nppsapwood   (donc)
                                  
      cpatch%today_nppcroot     (recc) = cpatch%today_nppcroot     (recc)                  &
                                       + cpatch%today_nppcroot     (donc)
                                  
      cpatch%today_nppseeds     (recc) = cpatch%today_nppseeds     (recc)                  &
                                       + cpatch%today_nppseeds     (donc)
                                  
      cpatch%today_nppwood      (recc) = cpatch%today_nppwood      (recc)                  &
                                       + cpatch%today_nppwood      (donc)
                                  
      cpatch%today_nppdaily     (recc) = cpatch%today_nppdaily     (recc)                  &
                                       + cpatch%today_nppdaily     (donc)
                                  
      cpatch%today_gpp_pot      (recc) = cpatch%today_gpp_pot      (recc)                  &
                                       + cpatch%today_gpp_pot      (donc)

      cpatch%today_gpp_lightmax (recc) = cpatch%today_gpp_lightmax (recc)                  &
                                       + cpatch%today_gpp_lightmax (donc)

      cpatch%today_gpp_moistmax (recc) = cpatch%today_gpp_moistmax (recc)                  &
                                       + cpatch%today_gpp_moistmax (donc)

      cpatch%today_gpp_mlmax    (recc) = cpatch%today_gpp_mlmax    (recc)                  &
                                       + cpatch%today_gpp_mlmax    (donc)

      cpatch%today_leaf_resp    (recc) = cpatch%today_leaf_resp    (recc)                  &
                                       + cpatch%today_leaf_resp    (donc)

      cpatch%today_root_resp    (recc) = cpatch%today_root_resp    (recc)                  &
                                       + cpatch%today_root_resp    (donc)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Fuse the leaf surface and intenal properties.  Since they are intensive         !
      ! properties, they are scaled by the number of plants.  These numbers are diagnostic !
      ! and this should be used for the output only.                                       !
      !------------------------------------------------------------------------------------!
      cpatch%lsfc_shv_open  (recc) = cpatch%lsfc_shv_open  (recc)  * rlai                  &
                                   + cpatch%lsfc_shv_open  (donc)  * dlai
      cpatch%lsfc_shv_closed(recc) = cpatch%lsfc_shv_closed(recc)  * rlai                  &
                                   + cpatch%lsfc_shv_closed(donc)  * dlai
      cpatch%lsfc_co2_open  (recc) = cpatch%lsfc_co2_open  (recc)  * rlai                  &
                                   + cpatch%lsfc_co2_open  (donc)  * dlai
      cpatch%lsfc_co2_closed(recc) = cpatch%lsfc_co2_closed(recc)  * rlai                  &
                                   + cpatch%lsfc_co2_closed(donc)  * dlai
      cpatch%lint_co2_open  (recc) = cpatch%lint_co2_open  (recc)  * rlai                  &
                                   + cpatch%lint_co2_open  (donc)  * dlai
      cpatch%lint_co2_closed(recc) = cpatch%lint_co2_closed(recc)  * rlai                  &
                                   + cpatch%lint_co2_closed(donc)  * dlai
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Fuse mortality rates.  By definition, mortality rate M for is given by:         !
      !                                                                                    !
      !                                    ln(A) - ln(N)                                   !
      !                               m = ---------------                                  !
      !                                          dt                                        !
      !                                                                                    !
      ! where A is the population that was previously alive, and N is the population that  !
      ! survived the mortality.  The cohorts represent a group of individuals with the     !
      ! same size and PFT, so they don't mix new recruits and old plants. Therefore, we    !
      ! can assume that N is actually nplant.  We don't know A, if the mortality rate is   !
      ! assumed constant during the interval dt, A = N * exp(m dt).                        !
      !                                                                                    !
      ! For fusion we don't really care about dt, so any number will do as long as it is   !
      ! the same for both cohorts.  With these assumptions, the mortality rate for the     !
      ! fused cohort mf is:                                                                !
      !                                                                                    !
      !  mf   =  ln (Ad+Ar) - ln(Nd+Nr) = ln[Nd*exp(md) + Nr*exp(mr)] - ln[Nd+Nr]          !
      !                                                                                    !
      !             / Nd*exp(md) + Nr*exp(mr) \                                            !
      !  mf   =  ln |-------------------------|                                            !
      !             \        Nd + Nr          /                                            !
      !------------------------------------------------------------------------------------!
      do imty=1,n_mort
         exp_mort_donc = exp(max(lnexp_min,min(lnexp_max,cpatch%mort_rate(imty,donc))))
         exp_mort_recc = exp(max(lnexp_min,min(lnexp_max,cpatch%mort_rate(imty,recc))))

         cpatch%mort_rate(imty,recc) = log( rnplant * exp_mort_recc                        &
                                          + dnplant * exp_mort_donc )
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Light level.  Using the intensive way of fusing.                                !
      !------------------------------------------------------------------------------------!
      cpatch%light_level     (recc) = cpatch%light_level      (recc) * rnplant             &
                                    + cpatch%light_level      (donc) * dnplant
      cpatch%light_level_beam(recc) = cpatch%light_level_beam (recc) * rnplant             &
                                    + cpatch%light_level_beam (donc) * dnplant
      cpatch%light_level_diff(recc) = cpatch%light_level_diff (recc) * rnplant             &
                                    + cpatch%light_level_diff (donc) * dnplant
      cpatch%par_level_beam  (recc) = cpatch%par_level_beam   (recc) * rnplant             &
                                    + cpatch%par_level_beam   (donc) * dnplant
      cpatch%par_level_diffd (recc) = cpatch%par_level_diffd  (recc) * rnplant             &
                                    + cpatch%par_level_diffd  (donc) * dnplant
      cpatch%par_level_diffu (recc) = cpatch%par_level_diffu  (recc) * rnplant             &
                                    + cpatch%par_level_diffu  (donc) * dnplant
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Long-trm respiration terms are in kgC/plant/yr, use nplant to fuse them.        !
      !------------------------------------------------------------------------------------!
      cpatch%leaf_growth_resp   (recc) = cpatch%leaf_growth_resp   (recc) * rnplant        &
                                       + cpatch%leaf_growth_resp   (donc) * dnplant
      cpatch%root_growth_resp   (recc) = cpatch%root_growth_resp   (recc) * rnplant        &
                                       + cpatch%root_growth_resp   (donc) * dnplant
      cpatch%sapa_growth_resp   (recc) = cpatch%sapa_growth_resp   (recc) * rnplant        &
                                       + cpatch%sapa_growth_resp   (donc) * dnplant
      cpatch%sapb_growth_resp   (recc) = cpatch%sapb_growth_resp   (recc) * rnplant        &
                                       + cpatch%sapb_growth_resp   (donc) * dnplant
      cpatch%storage_respiration(recc) = cpatch%storage_respiration(recc) * rnplant        &
                                       + cpatch%storage_respiration(donc) * dnplant
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Water demand is in kg/m2_leaf/s, so we scale them by LAI.  Water supply is in   !
      ! kg/m2_ground/s, so we just add them.                                               !
      !------------------------------------------------------------------------------------!
      cpatch%psi_open    (recc) = cpatch%psi_open  (recc) * rlai                           &
                                + cpatch%psi_open  (donc) * dlai
      cpatch%psi_closed  (recc) = cpatch%psi_closed(recc) * rlai                           &
                                + cpatch%psi_closed(donc) * dlai
      cpatch%water_supply(recc) = cpatch%water_supply(recc) + cpatch%water_supply(donc)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Carbon demand is in kg_C/m2_leaf/s, so we scale them by LAI.  FSW and FSN are   !
      ! really related to leaves, so we scale them by LAI.                                 !
      !------------------------------------------------------------------------------------!
      cpatch%A_open  (recc)     = cpatch%A_open  (recc) * rlai                             &
                                + cpatch%A_open  (donc) * dlai
      cpatch%A_closed(recc)     = cpatch%A_closed(recc) * rlai                             &
                                + cpatch%A_closed(donc) * dlai
      cpatch%A_light (recc)     = cpatch%A_light (donc) * rlai                             &
                                + cpatch%A_light (recc) * dlai
      cpatch%A_rubp  (recc)     = cpatch%A_rubp  (donc) * rlai                             &
                                + cpatch%A_rubp  (recc) * dlai
      cpatch%A_co2   (recc)     = cpatch%A_co2   (donc) * rlai                             &
                                + cpatch%A_co2   (recc) * dlai
      cpatch%fsw     (recc)     = cpatch%fsw     (recc) * rlai                             &
                                + cpatch%fsw     (donc) * dlai
      cpatch%fsn     (recc)     = cpatch%fsn     (recc) * rlai                             &
                                + cpatch%fsn     (donc) * dlai
      cpatch%fs_open (recc)     = cpatch%fsw(recc) * cpatch%fsn(recc)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Merge biomass and basal area.  Contrary to the patch/site/polygon levels,       !
      ! these variables are "intensive" (or per plant) at the cohort level, so we must     !
      ! average them.                                                                      !
      !------------------------------------------------------------------------------------!
      cpatch%agb      (recc) = cpatch%agb      (recc) * rnplant                            &
                             + cpatch%agb      (donc) * dnplant
      cpatch%basarea  (recc) = cpatch%basarea  (recc) * rnplant                            &
                             + cpatch%basarea  (donc) * dnplant
      cpatch%dagb_dt  (recc) = cpatch%dagb_dt  (recc) * rnplant                            &
                             + cpatch%dagb_dt  (donc) * dnplant
      cpatch%dba_dt   (recc) = cpatch%dba_dt   (recc) * rnplant                            &
                             + cpatch%dba_dt   (donc) * dnplant
      cpatch%dlndbh_dt(recc) = cpatch%dlndbh_dt(recc) * rnplant                            &
                             + cpatch%dlndbh_dt(donc) * dnplant
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Updating the tendency of plant density and the relative mortality rate.  The   !
      ! first is extensive (i.e. per unit of area), so it must be added, not scaled.  The  !
      ! second is intensive so it must be scaled.                                          !
      !------------------------------------------------------------------------------------!
      cpatch%monthly_dndt  (recc) = cpatch%monthly_dndt  (recc)                            &
                                  + cpatch%monthly_dndt  (donc)
      cpatch%monthly_dlnndt(recc) = cpatch%monthly_dlnndt(recc) * rnplant                  &
                                  + cpatch%monthly_dlnndt(donc) * dnplant
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update the carbon fluxes. They are fluxes per unit of area, so they should be  !
      ! added, not scaled.                                                                 !
      !------------------------------------------------------------------------------------!
      cpatch%gpp             (recc) = cpatch%gpp             (recc)                        &
                                    + cpatch%gpp             (donc)
      cpatch%leaf_respiration(recc) = cpatch%leaf_respiration(recc)                        &
                                    + cpatch%leaf_respiration(donc)
      cpatch%root_respiration(recc) = cpatch%root_respiration(recc)                        &
                                    + cpatch%root_respiration(donc)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Potential available water and elongation factor can be consider "intensive"    !
      ! variable water.                                                                    !
      !------------------------------------------------------------------------------------!
      cpatch%paw_avg(recc) = cpatch%paw_avg(recc) * rnplant + cpatch%paw_avg(donc) * dnplant
      cpatch%elongf (recc) = cpatch%elongf(recc)  * rnplant + cpatch%elongf (donc) * dnplant
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Light-phenology characteristics (MLO I am not sure if they should be scaled by  !
      ! nplant or LAI, it seems LAI would make more sense...).                             !
      !------------------------------------------------------------------------------------!
      cpatch%turnover_amp(recc) = cpatch%turnover_amp(recc) * rnplant                      &
                                + cpatch%turnover_amp(donc) * dnplant
      cpatch%llspan      (recc) = cpatch%llspan      (recc) * rnplant                      &
                                + cpatch%llspan      (donc) * dnplant
      cpatch%vm_bar      (recc) = cpatch%vm_bar      (recc) * rnplant                      &
                                + cpatch%vm_bar      (donc) * dnplant
      cpatch%sla         (recc) = cpatch%sla         (recc) * rnplant                      &
                                + cpatch%sla         (donc) * dnplant
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
      !    Fast averages.                                                                  !
      !------------------------------------------------------------------------------------!
      if (.not. fuse_initial) then
         !---------------------------------------------------------------------------------!
         !    Intensive variables, scaled by plant density.                                !
         !---------------------------------------------------------------------------------!
         cpatch%fmean_gpp             (recc) = cpatch%fmean_gpp             (recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_gpp             (donc)         &
                                             * dnplant
         cpatch%fmean_npp             (recc) = cpatch%fmean_npp             (recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_npp             (donc)         &
                                             * dnplant
         cpatch%fmean_leaf_resp       (recc) = cpatch%fmean_leaf_resp       (recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_leaf_resp       (donc)         &
                                             * dnplant
         cpatch%fmean_root_resp       (recc) = cpatch%fmean_root_resp       (recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_root_resp       (donc)         &
                                             * dnplant
         cpatch%fmean_leaf_growth_resp(recc) = cpatch%fmean_leaf_growth_resp(recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_leaf_growth_resp(donc)         &
                                             * dnplant
         cpatch%fmean_root_growth_resp(recc) = cpatch%fmean_root_growth_resp(recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_root_growth_resp(donc)         &
                                             * dnplant
         cpatch%fmean_sapa_growth_resp(recc) = cpatch%fmean_sapa_growth_resp(recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_sapa_growth_resp(donc)         &
                                             * dnplant
         cpatch%fmean_sapb_growth_resp(recc) = cpatch%fmean_sapb_growth_resp(recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_sapb_growth_resp(donc)         &
                                             * dnplant
         cpatch%fmean_storage_resp    (recc) = cpatch%fmean_storage_resp    (recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_storage_resp    (donc)         &
                                             * dnplant
         cpatch%fmean_plresp          (recc) = cpatch%fmean_plresp          (recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_plresp          (donc)         &
                                             * dnplant
         cpatch%fmean_light_level     (recc) = cpatch%fmean_light_level     (recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_light_level     (donc)         &
                                             * dnplant
         cpatch%fmean_light_level_beam(recc) = cpatch%fmean_light_level_beam(recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_light_level_beam(donc)         &
                                             * dnplant
         cpatch%fmean_light_level_diff(recc) = cpatch%fmean_light_level_diff(recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_light_level_diff(donc)         &
                                             * dnplant
         cpatch%fmean_par_level_beam  (recc) = cpatch%fmean_par_level_beam  (recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_par_level_beam  (donc)         &
                                             * dnplant
         cpatch%fmean_par_level_diffd (recc) = cpatch%fmean_par_level_diffd (recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_par_level_diffd (donc)         &
                                             * dnplant
         cpatch%fmean_par_level_diffu (recc) = cpatch%fmean_par_level_diffu (recc)         &
                                             * rnplant                                     &
                                             + cpatch%fmean_par_level_diffu (donc)         &
                                             * dnplant
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Intensive variables, scaled by LAI.                                          !
         !---------------------------------------------------------------------------------!
         cpatch%fmean_leaf_gsw        (recc) = cpatch%fmean_leaf_gsw        (recc) * rlai  &
                                             + cpatch%fmean_leaf_gsw        (donc) * dlai
         cpatch%fmean_leaf_gbw        (recc) = cpatch%fmean_leaf_gbw        (recc) * rlai  &
                                             + cpatch%fmean_leaf_gbw        (donc) * dlai
         cpatch%fmean_fs_open         (recc) = cpatch%fmean_fs_open         (recc) * rlai  &
                                             + cpatch%fmean_fs_open         (donc) * dlai
         cpatch%fmean_fsw             (recc) = cpatch%fmean_fsw             (recc) * rlai  &
                                             + cpatch%fmean_fsw             (donc) * dlai
         cpatch%fmean_fsn             (recc) = cpatch%fmean_fsn             (recc) * rlai  &
                                             + cpatch%fmean_fsn             (donc) * dlai
         cpatch%fmean_A_open          (recc) = cpatch%fmean_A_open          (recc) * rlai  &
                                             + cpatch%fmean_A_open          (donc) * dlai
         cpatch%fmean_A_closed        (recc) = cpatch%fmean_A_closed        (recc) * rlai  &
                                             + cpatch%fmean_A_closed        (donc) * dlai
         cpatch%fmean_A_net           (recc) = cpatch%fmean_A_net           (recc) * rlai  &
                                             + cpatch%fmean_A_net           (donc) * dlai
         cpatch%fmean_A_light         (recc) = cpatch%fmean_A_light         (recc) * rlai  &
                                             + cpatch%fmean_A_light         (donc) * dlai
         cpatch%fmean_A_rubp          (recc) = cpatch%fmean_A_rubp          (recc) * rlai  &
                                             + cpatch%fmean_A_rubp          (donc) * dlai
         cpatch%fmean_A_co2           (recc) = cpatch%fmean_A_co2           (recc) * rlai  &
                                             + cpatch%fmean_A_co2           (donc) * dlai
         cpatch%fmean_psi_open        (recc) = cpatch%fmean_psi_open        (recc) * rlai  &
                                             + cpatch%fmean_psi_open        (donc) * dlai
         cpatch%fmean_psi_closed      (recc) = cpatch%fmean_psi_closed      (recc) * rlai  &
                                             + cpatch%fmean_psi_closed      (donc) * dlai
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Intensive variables, scaled by WAI.                                          !
         !---------------------------------------------------------------------------------!
         cpatch%fmean_wood_gbw        (recc) = cpatch%fmean_wood_gbw        (recc) * rwai  &
                                             + cpatch%fmean_wood_gbw        (donc) * dwai
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Extensive variables.                                                         !
         !---------------------------------------------------------------------------------!
         cpatch%fmean_leaf_energy     (recc) = cpatch%fmean_leaf_energy     (recc)         &
                                             + cpatch%fmean_leaf_energy     (donc)
         cpatch%fmean_leaf_water      (recc) = cpatch%fmean_leaf_water      (recc)         &
                                             + cpatch%fmean_leaf_water      (donc)
         cpatch%fmean_leaf_hcap       (recc) = cpatch%fmean_leaf_hcap       (recc)         &
                                             + cpatch%fmean_leaf_hcap       (donc)
         cpatch%fmean_wood_energy     (recc) = cpatch%fmean_wood_energy     (recc)         &
                                             + cpatch%fmean_wood_energy     (donc)
         cpatch%fmean_wood_water      (recc) = cpatch%fmean_wood_water      (recc)         &
                                             + cpatch%fmean_wood_water      (donc)
         cpatch%fmean_wood_hcap       (recc) = cpatch%fmean_wood_hcap       (recc)         &
                                             + cpatch%fmean_wood_hcap       (donc)
         cpatch%fmean_water_supply    (recc) = cpatch%fmean_water_supply    (recc)         &
                                             + cpatch%fmean_water_supply    (donc)
         cpatch%fmean_par_l           (recc) = cpatch%fmean_par_l           (recc)         &
                                             + cpatch%fmean_par_l           (donc)
         cpatch%fmean_par_l_beam      (recc) = cpatch%fmean_par_l_beam      (recc)         &
                                             + cpatch%fmean_par_l_beam      (donc)
         cpatch%fmean_par_l_diff      (recc) = cpatch%fmean_par_l_diff      (recc)         &
                                             + cpatch%fmean_par_l_diff      (donc)
         cpatch%fmean_rshort_l        (recc) = cpatch%fmean_rshort_l        (recc)         &
                                             + cpatch%fmean_rshort_l        (donc)
         cpatch%fmean_rlong_l         (recc) = cpatch%fmean_rlong_l         (recc)         &
                                             + cpatch%fmean_rlong_l         (donc)
         cpatch%fmean_sensible_lc     (recc) = cpatch%fmean_sensible_lc     (recc)         &
                                             + cpatch%fmean_sensible_lc     (donc)
         cpatch%fmean_vapor_lc        (recc) = cpatch%fmean_vapor_lc        (recc)         &
                                             + cpatch%fmean_vapor_lc        (donc)
         cpatch%fmean_transp          (recc) = cpatch%fmean_transp          (recc)         &
                                             + cpatch%fmean_transp          (donc)
         cpatch%fmean_intercepted_al  (recc) = cpatch%fmean_intercepted_al  (recc)         &
                                             + cpatch%fmean_intercepted_al  (donc)
         cpatch%fmean_wshed_lg        (recc) = cpatch%fmean_wshed_lg        (recc)         &
                                             + cpatch%fmean_wshed_lg        (donc)
         cpatch%fmean_rshort_w        (recc) = cpatch%fmean_rshort_w        (recc)         &
                                             + cpatch%fmean_rshort_w        (donc)
         cpatch%fmean_rlong_w         (recc) = cpatch%fmean_rlong_w         (recc)         &
                                             + cpatch%fmean_rlong_w         (donc)
         cpatch%fmean_rad_profile   (:,recc) = cpatch%fmean_rad_profile   (:,recc)         &
                                             + cpatch%fmean_rad_profile   (:,donc)
         cpatch%fmean_sensible_wc     (recc) = cpatch%fmean_sensible_wc     (recc)         &
                                             + cpatch%fmean_sensible_wc     (donc)
         cpatch%fmean_vapor_wc        (recc) = cpatch%fmean_vapor_wc        (recc)         &
                                             + cpatch%fmean_vapor_wc        (donc)
         cpatch%fmean_intercepted_aw  (recc) = cpatch%fmean_intercepted_aw  (recc)         &
                                             + cpatch%fmean_intercepted_aw  (donc)
         cpatch%fmean_wshed_wg        (recc) = cpatch%fmean_wshed_wg        (recc)         &
                                             + cpatch%fmean_wshed_wg        (donc)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      We update temperature and liquid water fraction.  We check whether the     !
         ! heat capacity is non-zero.  If it is a normal number, use the standard thermo-  !
         ! dynamic library, otherwise, average temperature, this is probably a blend of    !
         ! tiny cohorts that couldn't be solved, or the wood is not solved.                !
         !---------------------------------------------------------------------------------!
         !------ Leaf. --------------------------------------------------------------------!
         if ( cpatch%fmean_leaf_hcap(recc) > 0. ) then
            !----- Update temperature and liquid fraction using standard thermodynamics. --!
            call uextcm2tl(cpatch%fmean_leaf_energy(recc),cpatch%fmean_leaf_water(recc)    &
                          ,cpatch%fmean_leaf_hcap  (recc),cpatch%fmean_leaf_temp (recc)    &
                          ,cpatch%fmean_leaf_fliq  (recc))
            !------------------------------------------------------------------------------!


            !----- Scale vapour pressure deficit using LAI. -------------------------------!
            cpatch%fmean_leaf_vpdef   (recc) = cpatch%fmean_leaf_vpdef      (recc) * rlai  &
                                             + cpatch%fmean_leaf_vpdef      (donc) * dlai
            !------------------------------------------------------------------------------!
         else
            !----- None of the cohorts has leaf biomass use nplant to scale them. ---------!
            cpatch%fmean_leaf_temp (recc) = cpatch%fmean_leaf_temp (recc) * rnplant        &
                                          + cpatch%fmean_leaf_temp (donc) * dnplant
            cpatch%fmean_leaf_fliq (recc) = 0.0
            cpatch%fmean_leaf_vpdef(recc) = cpatch%fmean_leaf_vpdef(recc) * rnplant        &
                                          + cpatch%fmean_leaf_vpdef(donc) * dnplant
            !------------------------------------------------------------------------------!
         end if
         !------ Wood. --------------------------------------------------------------------!
         if ( cpatch%fmean_wood_hcap(recc) > 0. ) then
            !----- Update temperature using the standard thermodynamics. ------------------!
            call uextcm2tl(cpatch%fmean_wood_energy(recc),cpatch%fmean_wood_water(recc)    &
                          ,cpatch%fmean_wood_hcap  (recc),cpatch%fmean_wood_temp (recc)    &
                          ,cpatch%fmean_wood_fliq  (recc))
            !------------------------------------------------------------------------------!
         else                                                                              
            !----- Wood temperature can't be found using uextcm2tl (singularity). ---------!
            cpatch%fmean_wood_temp(recc) = cpatch%fmean_wood_temp(recc) * rnplant          &
                                         + cpatch%fmean_wood_temp(donc) * dnplant
            cpatch%fmean_wood_fliq(recc) = 0.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Daily means.                                                                    !
      !------------------------------------------------------------------------------------!
      if (writing_long .and. (.not. fuse_initial)) then

         !---------------------------------------------------------------------------------!
         !    Intensive variables, scaled by plant density.                                !
         !---------------------------------------------------------------------------------!
         cpatch%dmean_nppleaf         (recc) = cpatch%dmean_nppleaf         (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_nppleaf         (donc)         &
                                             * dnplant
         cpatch%dmean_nppfroot        (recc) = cpatch%dmean_nppfroot        (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_nppfroot        (donc)         &
                                             * dnplant
         cpatch%dmean_nppsapwood      (recc) = cpatch%dmean_nppsapwood      (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_nppsapwood      (donc)         &
                                             * dnplant
         cpatch%dmean_nppcroot        (recc) = cpatch%dmean_nppcroot        (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_nppcroot        (donc)         &
                                             * dnplant
         cpatch%dmean_nppseeds        (recc) = cpatch%dmean_nppseeds        (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_nppseeds        (donc)         &
                                             * dnplant
         cpatch%dmean_nppwood         (recc) = cpatch%dmean_nppwood         (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_nppwood         (donc)         &
                                             * dnplant
         cpatch%dmean_nppdaily        (recc) = cpatch%dmean_nppdaily        (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_nppdaily        (donc)         &
                                             * dnplant
         cpatch%dmean_gpp             (recc) = cpatch%dmean_gpp             (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_gpp             (donc)         &
                                             * dnplant
         cpatch%dmean_npp             (recc) = cpatch%dmean_npp             (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_npp             (donc)         &
                                             * dnplant
         cpatch%dmean_leaf_resp       (recc) = cpatch%dmean_leaf_resp       (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_leaf_resp       (donc)         &
                                             * dnplant
         cpatch%dmean_root_resp       (recc) = cpatch%dmean_root_resp       (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_root_resp       (donc)         &
                                             * dnplant
         cpatch%dmean_leaf_growth_resp(recc) = cpatch%dmean_leaf_growth_resp(recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_leaf_growth_resp(donc)         &
                                             * dnplant
         cpatch%dmean_root_growth_resp(recc) = cpatch%dmean_root_growth_resp(recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_root_growth_resp(donc)         &
                                             * dnplant
         cpatch%dmean_sapa_growth_resp(recc) = cpatch%dmean_sapa_growth_resp(recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_sapa_growth_resp(donc)         &
                                             * dnplant
         cpatch%dmean_sapb_growth_resp(recc) = cpatch%dmean_sapb_growth_resp(recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_sapb_growth_resp(donc)         &
                                             * dnplant
         cpatch%dmean_storage_resp    (recc) = cpatch%dmean_storage_resp    (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_storage_resp    (donc)         &
                                             * dnplant
         cpatch%dmean_plresp          (recc) = cpatch%dmean_plresp          (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_plresp          (donc)         &
                                             * dnplant
         cpatch%dmean_light_level     (recc) = cpatch%dmean_light_level     (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_light_level     (donc)         &
                                             * dnplant
         cpatch%dmean_light_level_beam(recc) = cpatch%dmean_light_level_beam(recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_light_level_beam(donc)         &
                                             * dnplant
         cpatch%dmean_light_level_diff(recc) = cpatch%dmean_light_level_diff(recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_light_level_diff(donc)         &
                                             * dnplant
         cpatch%dmean_par_level_beam  (recc) = cpatch%dmean_par_level_beam  (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_par_level_beam  (donc)         &
                                             * dnplant
         cpatch%dmean_par_level_diffd (recc) = cpatch%dmean_par_level_diffd (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_par_level_diffd (donc)         &
                                             * dnplant
         cpatch%dmean_par_level_diffu (recc) = cpatch%dmean_par_level_diffu (recc)         &
                                             * rnplant                                     &
                                             + cpatch%dmean_par_level_diffu (donc)         &
                                             * dnplant
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Intensive variables, scaled by LAI.                                          !
         !---------------------------------------------------------------------------------!
         cpatch%dmean_leaf_gsw        (recc) = cpatch%dmean_leaf_gsw        (recc) * rlai  &
                                             + cpatch%dmean_leaf_gsw        (donc) * dlai
         cpatch%dmean_leaf_gbw        (recc) = cpatch%dmean_leaf_gbw        (recc) * rlai  &
                                             + cpatch%dmean_leaf_gbw        (donc) * dlai
         cpatch%dmean_fs_open         (recc) = cpatch%dmean_fs_open         (recc) * rlai  &
                                             + cpatch%dmean_fs_open         (donc) * dlai
         cpatch%dmean_fsw             (recc) = cpatch%dmean_fsw             (recc) * rlai  &
                                             + cpatch%dmean_fsw             (donc) * dlai
         cpatch%dmean_fsn             (recc) = cpatch%dmean_fsn             (recc) * rlai  &
                                             + cpatch%dmean_fsn             (donc) * dlai
         cpatch%dmean_A_open          (recc) = cpatch%dmean_A_open          (recc) * rlai  &
                                             + cpatch%dmean_A_open          (donc) * dlai
         cpatch%dmean_A_closed        (recc) = cpatch%dmean_A_closed        (recc) * rlai  &
                                             + cpatch%dmean_A_closed        (donc) * dlai
         cpatch%dmean_A_net           (recc) = cpatch%dmean_A_net           (recc) * rlai  &
                                             + cpatch%dmean_A_net           (donc) * dlai
         cpatch%dmean_A_light         (recc) = cpatch%dmean_A_light         (recc) * rlai  &
                                             + cpatch%dmean_A_light         (donc) * dlai
         cpatch%dmean_A_rubp          (recc) = cpatch%dmean_A_rubp          (recc) * rlai  &
                                             + cpatch%dmean_A_rubp          (donc) * dlai
         cpatch%dmean_A_co2           (recc) = cpatch%dmean_A_co2           (recc) * rlai  &
                                             + cpatch%dmean_A_co2           (donc) * dlai
         cpatch%dmean_psi_open        (recc) = cpatch%dmean_psi_open        (recc) * rlai  &
                                             + cpatch%dmean_psi_open        (donc) * dlai
         cpatch%dmean_psi_closed      (recc) = cpatch%dmean_psi_closed      (recc) * rlai  &
                                             + cpatch%dmean_psi_closed      (donc) * dlai
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Intensive variables, scaled by WAI.                                          !
         !---------------------------------------------------------------------------------!
         cpatch%dmean_wood_gbw        (recc) = cpatch%dmean_wood_gbw        (recc) * rwai  &
                                             + cpatch%dmean_wood_gbw        (donc) * dwai
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Extensive variables.                                                         !
         !---------------------------------------------------------------------------------!
         cpatch%dmean_leaf_energy     (recc) = cpatch%dmean_leaf_energy     (recc)         &
                                             + cpatch%dmean_leaf_energy     (donc)
         cpatch%dmean_leaf_water      (recc) = cpatch%dmean_leaf_water      (recc)         &
                                             + cpatch%dmean_leaf_water      (donc)
         cpatch%dmean_leaf_hcap       (recc) = cpatch%dmean_leaf_hcap       (recc)         &
                                             + cpatch%dmean_leaf_hcap       (donc)
         cpatch%dmean_wood_energy     (recc) = cpatch%dmean_wood_energy     (recc)         &
                                             + cpatch%dmean_wood_energy     (donc)
         cpatch%dmean_wood_water      (recc) = cpatch%dmean_wood_water      (recc)         &
                                             + cpatch%dmean_wood_water      (donc)
         cpatch%dmean_wood_hcap       (recc) = cpatch%dmean_wood_hcap       (recc)         &
                                             + cpatch%dmean_wood_hcap       (donc)
         cpatch%dmean_water_supply    (recc) = cpatch%dmean_water_supply    (recc)         &
                                             + cpatch%dmean_water_supply    (donc)
         cpatch%dmean_par_l           (recc) = cpatch%dmean_par_l           (recc)         &
                                             + cpatch%dmean_par_l           (donc)
         cpatch%dmean_par_l_beam      (recc) = cpatch%dmean_par_l_beam      (recc)         &
                                             + cpatch%dmean_par_l_beam      (donc)
         cpatch%dmean_par_l_diff      (recc) = cpatch%dmean_par_l_diff      (recc)         &
                                             + cpatch%dmean_par_l_diff      (donc)
         cpatch%dmean_rshort_l        (recc) = cpatch%dmean_rshort_l        (recc)         &
                                             + cpatch%dmean_rshort_l        (donc)
         cpatch%dmean_rlong_l         (recc) = cpatch%dmean_rlong_l         (recc)         &
                                             + cpatch%dmean_rlong_l         (donc)
         cpatch%dmean_sensible_lc     (recc) = cpatch%dmean_sensible_lc     (recc)         &
                                             + cpatch%dmean_sensible_lc     (donc)
         cpatch%dmean_vapor_lc        (recc) = cpatch%dmean_vapor_lc        (recc)         &
                                             + cpatch%dmean_vapor_lc        (donc)
         cpatch%dmean_transp          (recc) = cpatch%dmean_transp          (recc)         &
                                             + cpatch%dmean_transp          (donc)
         cpatch%dmean_intercepted_al  (recc) = cpatch%dmean_intercepted_al  (recc)         &
                                             + cpatch%dmean_intercepted_al  (donc)
         cpatch%dmean_wshed_lg        (recc) = cpatch%dmean_wshed_lg        (recc)         &
                                             + cpatch%dmean_wshed_lg        (donc)
         cpatch%dmean_rshort_w        (recc) = cpatch%dmean_rshort_w        (recc)         &
                                             + cpatch%dmean_rshort_w        (donc)
         cpatch%dmean_rlong_w         (recc) = cpatch%dmean_rlong_w         (recc)         &
                                             + cpatch%dmean_rlong_w         (donc)
         cpatch%dmean_rad_profile   (:,recc) = cpatch%dmean_rad_profile   (:,recc)         &
                                             + cpatch%dmean_rad_profile   (:,donc)
         cpatch%dmean_sensible_wc     (recc) = cpatch%dmean_sensible_wc     (recc)         &
                                             + cpatch%dmean_sensible_wc     (donc)
         cpatch%dmean_vapor_wc        (recc) = cpatch%dmean_vapor_wc        (recc)         &
                                             + cpatch%dmean_vapor_wc        (donc)
         cpatch%dmean_intercepted_aw  (recc) = cpatch%dmean_intercepted_aw  (recc)         &
                                             + cpatch%dmean_intercepted_aw  (donc)
         cpatch%dmean_wshed_wg        (recc) = cpatch%dmean_wshed_wg        (recc)         &
                                             + cpatch%dmean_wshed_wg        (donc)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      We update temperature and liquid water fraction.  We check whether the     !
         ! heat capacity is non-zero.  If it is a normal number, use the standard thermo-  !
         ! dynamic library, otherwise, average temperature, this is probably a blend of    !
         ! tiny cohorts that couldn't be solved, or the wood is not solved.                !
         !---------------------------------------------------------------------------------!
         !------ Leaf. --------------------------------------------------------------------!
         if ( cpatch%dmean_leaf_hcap(recc) > 0. ) then
            !----- Update temperature and liquid fraction using standard thermodynamics. --!
            call uextcm2tl(cpatch%dmean_leaf_energy(recc),cpatch%dmean_leaf_water(recc)    &
                          ,cpatch%dmean_leaf_hcap  (recc),cpatch%dmean_leaf_temp (recc)    &
                          ,cpatch%dmean_leaf_fliq  (recc))                                                  
            !------------------------------------------------------------------------------!


            !----- Scale vapour pressure deficit using LAI. -------------------------------!
            cpatch%dmean_leaf_vpdef   (recc) = cpatch%dmean_leaf_vpdef      (recc) * rlai  &
                                             + cpatch%dmean_leaf_vpdef      (donc) * dlai
            !------------------------------------------------------------------------------!
         else
            !----- None of the cohorts has leaf biomass use nplant to scale them. ---------!
            cpatch%dmean_leaf_temp (recc) = cpatch%dmean_leaf_temp (recc) * rnplant        &
                                          + cpatch%dmean_leaf_temp (donc) * dnplant
            cpatch%dmean_leaf_fliq (recc) = 0.0
            cpatch%dmean_leaf_vpdef(recc) = cpatch%dmean_leaf_vpdef(recc) * rnplant        &
                                          + cpatch%dmean_leaf_vpdef(donc) * dnplant
            !------------------------------------------------------------------------------!
         end if                                                                            
         !------ Wood. --------------------------------------------------------------------!
         if ( cpatch%dmean_wood_hcap(recc) > 0. ) then                                     
            !----- Update temperature using the standard thermodynamics. ------------------!
            call uextcm2tl(cpatch%dmean_wood_energy(recc),cpatch%dmean_wood_water(recc)    &
                          ,cpatch%dmean_wood_hcap  (recc),cpatch%dmean_wood_temp (recc)    &
                          ,cpatch%dmean_wood_fliq  (recc))
            !------------------------------------------------------------------------------!
         else                                                                              
            !----- Wood temperature can't be found using uextcm2tl (singularity). ---------!
            cpatch%dmean_wood_temp(recc) = cpatch%dmean_wood_temp(recc) * rnplant          &
                                         + cpatch%dmean_wood_temp(donc) * dnplant
            cpatch%dmean_wood_fliq(recc) = 0.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Monthly means.                                                                  !
      !------------------------------------------------------------------------------------!
      if (writing_eorq .and. (.not. fuse_initial)) then
         !---------------------------------------------------------------------------------!
         !     First we merge the squares, as they require the means.                      !
         !---------------------------------------------------------------------------------!
         !----- Intensive. ----------------------------------------------------------------!
         cpatch%mmsqu_gpp        (recc) = fuse_msqu( cpatch%mmean_gpp        (recc)        &
                                                   , cpatch%mmsqu_gpp        (recc)        &
                                                   , rnplant                               &
                                                   , cpatch%mmean_gpp        (donc)        &
                                                   , cpatch%mmsqu_gpp        (donc)        &
                                                   , cpatch%nplant           (donc)        &
                                                   , corr_cohort, .false.)
         cpatch%mmsqu_npp        (recc) = fuse_msqu( cpatch%mmean_npp        (recc)        &
                                                   , cpatch%mmsqu_npp        (recc)        &
                                                   , rnplant                               &
                                                   , cpatch%mmean_npp        (donc)        &
                                                   , cpatch%mmsqu_npp        (donc)        &
                                                   , cpatch%nplant           (donc)        &
                                                   , corr_cohort, .false.)
         cpatch%mmsqu_plresp     (recc) = fuse_msqu( cpatch%mmean_plresp     (recc)        &
                                                   , cpatch%mmsqu_plresp     (recc)        &
                                                   , rnplant                               &
                                                   , cpatch%mmean_plresp     (donc)        &
                                                   , cpatch%mmsqu_plresp     (donc)        &
                                                   , dnplant                               &
                                                   , corr_cohort, .false.)
         !----- Extensive variables. ------------------------------------------------------!
         cpatch%mmsqu_sensible_lc(recc) = fuse_msqu( cpatch%mmean_sensible_lc(recc)        &
                                                   , cpatch%mmsqu_sensible_lc(recc)        &
                                                   , rnplant                               &
                                                   , cpatch%mmean_sensible_lc(donc)        &
                                                   , cpatch%mmsqu_sensible_lc(donc)        &
                                                   , dnplant                               &
                                                   , corr_cohort, .true. )
         cpatch%mmsqu_vapor_lc   (recc) = fuse_msqu( cpatch%mmean_vapor_lc   (recc)        &
                                                   , cpatch%mmsqu_vapor_lc   (recc)        &
                                                   , rnplant                               &
                                                   , cpatch%mmean_vapor_lc   (donc)        &
                                                   , cpatch%mmsqu_vapor_lc   (donc)        &
                                                   , dnplant                               &
                                                   , corr_cohort, .true. )
         cpatch%mmsqu_transp     (recc) = fuse_msqu( cpatch%mmean_transp     (recc)        &
                                                   , cpatch%mmsqu_transp     (recc)        &
                                                   , rnplant                               &
                                                   , cpatch%mmean_transp     (donc)        &
                                                   , cpatch%mmsqu_transp     (donc)        &
                                                   , dnplant                               &
                                                   , corr_cohort, .true. )
         cpatch%mmsqu_sensible_wc(recc) = fuse_msqu( cpatch%mmean_sensible_wc(recc)        &
                                                   , cpatch%mmsqu_sensible_wc(recc)        &
                                                   , rnplant                               &
                                                   , cpatch%mmean_sensible_wc(donc)        &
                                                   , cpatch%mmsqu_sensible_wc(donc)        &
                                                   , dnplant                               &
                                                   , corr_cohort, .true. )
         cpatch%mmsqu_vapor_wc   (recc) = fuse_msqu( cpatch%mmean_vapor_wc   (recc)        &
                                                   , cpatch%mmsqu_vapor_wc   (recc)        &
                                                   , rnplant                               &
                                                   , cpatch%mmean_vapor_wc   (donc)        &
                                                   , cpatch%mmsqu_vapor_wc   (donc)        &
                                                   , dnplant                               &
                                                   , corr_cohort, .true. )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Intensive variables, scaled by plant density.                                !
         !---------------------------------------------------------------------------------!
         cpatch%mmean_nppleaf         (recc) = cpatch%mmean_nppleaf         (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_nppleaf         (donc)         &
                                             * dnplant
         cpatch%mmean_nppfroot        (recc) = cpatch%mmean_nppfroot        (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_nppfroot        (donc)         &
                                             * dnplant
         cpatch%mmean_nppsapwood      (recc) = cpatch%mmean_nppsapwood      (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_nppsapwood      (donc)         &
                                             * dnplant
         cpatch%mmean_nppcroot        (recc) = cpatch%mmean_nppcroot        (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_nppcroot        (donc)         &
                                             * dnplant
         cpatch%mmean_nppseeds        (recc) = cpatch%mmean_nppseeds        (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_nppseeds        (donc)         &
                                             * dnplant
         cpatch%mmean_nppwood         (recc) = cpatch%mmean_nppwood         (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_nppwood         (donc)         &
                                             * dnplant
         cpatch%mmean_nppdaily        (recc) = cpatch%mmean_nppdaily        (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_nppdaily        (donc)         &
                                             * dnplant
         cpatch%mmean_gpp             (recc) = cpatch%mmean_gpp             (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_gpp             (donc)         &
                                             * dnplant
         cpatch%mmean_npp             (recc) = cpatch%mmean_npp             (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_npp             (donc)         &
                                             * dnplant
         cpatch%mmean_leaf_resp       (recc) = cpatch%mmean_leaf_resp       (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_leaf_resp       (donc)         &
                                             * dnplant
         cpatch%mmean_root_resp       (recc) = cpatch%mmean_root_resp       (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_root_resp       (donc)         &
                                             * dnplant
         cpatch%mmean_leaf_growth_resp(recc) = cpatch%mmean_leaf_growth_resp(recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_leaf_growth_resp(donc)         &
                                             * dnplant
         cpatch%mmean_root_growth_resp(recc) = cpatch%mmean_root_growth_resp(recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_root_growth_resp(donc)         &
                                             * dnplant
         cpatch%mmean_sapa_growth_resp(recc) = cpatch%mmean_sapa_growth_resp(recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_sapa_growth_resp(donc)         &
                                             * dnplant
         cpatch%mmean_sapb_growth_resp(recc) = cpatch%mmean_sapb_growth_resp(recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_sapb_growth_resp(donc)         &
                                             * dnplant
         cpatch%mmean_storage_resp    (recc) = cpatch%mmean_storage_resp    (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_storage_resp    (donc)         &
                                             * dnplant
         cpatch%mmean_plresp          (recc) = cpatch%mmean_plresp          (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_plresp          (donc)         &
                                             * dnplant
         cpatch%mmean_light_level     (recc) = cpatch%mmean_light_level     (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_light_level     (donc)         &
                                             * dnplant
         cpatch%mmean_light_level_beam(recc) = cpatch%mmean_light_level_beam(recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_light_level_beam(donc)         &
                                             * dnplant
         cpatch%mmean_light_level_diff(recc) = cpatch%mmean_light_level_diff(recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_light_level_diff(donc)         &
                                             * dnplant
         cpatch%mmean_bleaf           (recc) = cpatch%mmean_bleaf           (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_bleaf           (donc)         &
                                             * dnplant
         cpatch%mmean_broot           (recc) = cpatch%mmean_broot           (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_broot           (donc)         &
                                             * dnplant
         cpatch%mmean_bstorage        (recc) = cpatch%mmean_bstorage        (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_bstorage        (donc)         &
                                             * dnplant
         cpatch%mmean_leaf_maintenance(recc) = cpatch%mmean_leaf_maintenance(recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_leaf_maintenance(donc)         &
                                             * dnplant
         cpatch%mmean_root_maintenance(recc) = cpatch%mmean_root_maintenance(recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_root_maintenance(donc)         &
                                             * dnplant
         cpatch%mmean_leaf_drop       (recc) = cpatch%mmean_leaf_drop       (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_leaf_drop       (donc)         &
                                             * dnplant
         cpatch%mmean_cb              (recc) = cpatch%mmean_cb              (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_cb              (donc)         &
                                             * dnplant
         cpatch%mmean_par_level_beam  (recc) = cpatch%mmean_par_level_beam  (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_par_level_beam  (donc)         &
                                             * dnplant
         cpatch%mmean_par_level_diffd (recc) = cpatch%mmean_par_level_diffd (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_par_level_diffd (donc)         &
                                             * dnplant
         cpatch%mmean_par_level_diffu (recc) = cpatch%mmean_par_level_diffu (recc)         &
                                             * rnplant                                     &
                                             + cpatch%mmean_par_level_diffu (donc)         &
                                             * dnplant
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Intensive variables, scaled by LAI.                                          !
         !---------------------------------------------------------------------------------!
         cpatch%mmean_leaf_gsw        (recc) = cpatch%mmean_leaf_gsw        (recc) * rlai  &
                                             + cpatch%mmean_leaf_gsw        (donc) * dlai
         cpatch%mmean_leaf_gbw        (recc) = cpatch%mmean_leaf_gbw        (recc) * rlai  &
                                             + cpatch%mmean_leaf_gbw        (donc) * dlai
         cpatch%mmean_fs_open         (recc) = cpatch%mmean_fs_open         (recc) * rlai  &
                                             + cpatch%mmean_fs_open         (donc) * dlai
         cpatch%mmean_fsw             (recc) = cpatch%mmean_fsw             (recc) * rlai  &
                                             + cpatch%mmean_fsw             (donc) * dlai
         cpatch%mmean_fsn             (recc) = cpatch%mmean_fsn             (recc) * rlai  &
                                             + cpatch%mmean_fsn             (donc) * dlai
         cpatch%mmean_A_open          (recc) = cpatch%mmean_A_open          (recc) * rlai  &
                                             + cpatch%mmean_A_open          (donc) * dlai
         cpatch%mmean_A_closed        (recc) = cpatch%mmean_A_closed        (recc) * rlai  &
                                             + cpatch%mmean_A_closed        (donc) * dlai
         cpatch%mmean_A_net           (recc) = cpatch%mmean_A_net           (recc) * rlai  &
                                             + cpatch%mmean_A_net           (donc) * dlai
         cpatch%mmean_A_light         (recc) = cpatch%mmean_A_light         (recc) * rlai  &
                                             + cpatch%mmean_A_light         (donc) * dlai
         cpatch%mmean_A_rubp          (recc) = cpatch%mmean_A_rubp          (recc) * rlai  &
                                             + cpatch%mmean_A_rubp          (donc) * dlai
         cpatch%mmean_A_co2           (recc) = cpatch%mmean_A_co2           (recc) * rlai  &
                                             + cpatch%mmean_A_co2           (donc) * dlai
         cpatch%mmean_psi_open        (recc) = cpatch%mmean_psi_open        (recc) * rlai  &
                                             + cpatch%mmean_psi_open        (donc) * dlai
         cpatch%mmean_psi_closed      (recc) = cpatch%mmean_psi_closed      (recc) * rlai  &
                                             + cpatch%mmean_psi_closed      (donc) * dlai
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Intensive variables, scaled by WAI.                                          !
         !---------------------------------------------------------------------------------!
         cpatch%mmean_wood_gbw        (recc) = cpatch%mmean_wood_gbw        (recc) * rwai  &
                                             + cpatch%mmean_wood_gbw        (donc) * dwai
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Extensive variables.                                                         !
         !---------------------------------------------------------------------------------!
         cpatch%mmean_lai             (recc) = cpatch%mmean_lai             (recc)         &
                                             + cpatch%mmean_lai             (donc)
         cpatch%mmean_leaf_energy     (recc) = cpatch%mmean_leaf_energy     (recc)         &
                                             + cpatch%mmean_leaf_energy     (donc)
         cpatch%mmean_leaf_water      (recc) = cpatch%mmean_leaf_water      (recc)         &
                                             + cpatch%mmean_leaf_water      (donc)
         cpatch%mmean_leaf_hcap       (recc) = cpatch%mmean_leaf_hcap       (recc)         &
                                             + cpatch%mmean_leaf_hcap       (donc)
         cpatch%mmean_wood_energy     (recc) = cpatch%mmean_wood_energy     (recc)         &
                                             + cpatch%mmean_wood_energy     (donc)
         cpatch%mmean_wood_water      (recc) = cpatch%mmean_wood_water      (recc)         &
                                             + cpatch%mmean_wood_water      (donc)
         cpatch%mmean_wood_hcap       (recc) = cpatch%mmean_wood_hcap       (recc)         &
                                             + cpatch%mmean_wood_hcap       (donc)
         cpatch%mmean_water_supply    (recc) = cpatch%mmean_water_supply    (recc)         &
                                             + cpatch%mmean_water_supply    (donc)
         cpatch%mmean_par_l           (recc) = cpatch%mmean_par_l           (recc)         &
                                             + cpatch%mmean_par_l           (donc)
         cpatch%mmean_par_l_beam      (recc) = cpatch%mmean_par_l_beam      (recc)         &
                                             + cpatch%mmean_par_l_beam      (donc)
         cpatch%mmean_par_l_diff      (recc) = cpatch%mmean_par_l_diff      (recc)         &
                                             + cpatch%mmean_par_l_diff      (donc)

         cpatch%mmean_rshort_l        (recc) = cpatch%mmean_rshort_l        (recc)         &
                                             + cpatch%mmean_rshort_l        (donc)
         cpatch%mmean_rlong_l         (recc) = cpatch%mmean_rlong_l         (recc)         &
                                             + cpatch%mmean_rlong_l         (donc)
         cpatch%mmean_sensible_lc     (recc) = cpatch%mmean_sensible_lc     (recc)         &
                                             + cpatch%mmean_sensible_lc     (donc)
         cpatch%mmean_vapor_lc        (recc) = cpatch%mmean_vapor_lc        (recc)         &
                                             + cpatch%mmean_vapor_lc        (donc)
         cpatch%mmean_transp          (recc) = cpatch%mmean_transp          (recc)         &
                                             + cpatch%mmean_transp          (donc)
         cpatch%mmean_intercepted_al  (recc) = cpatch%mmean_intercepted_al  (recc)         &
                                             + cpatch%mmean_intercepted_al  (donc)
         cpatch%mmean_wshed_lg        (recc) = cpatch%mmean_wshed_lg        (recc)         &
                                             + cpatch%mmean_wshed_lg        (donc)
         cpatch%mmean_rshort_w        (recc) = cpatch%mmean_rshort_w        (recc)         &
                                             + cpatch%mmean_rshort_w        (donc)
         cpatch%mmean_rlong_w         (recc) = cpatch%mmean_rlong_w         (recc)         &
                                             + cpatch%mmean_rlong_w         (donc)
         cpatch%mmean_rad_profile   (:,recc) = cpatch%mmean_rad_profile   (:,recc)         &
                                             + cpatch%mmean_rad_profile   (:,donc)
         cpatch%mmean_sensible_wc     (recc) = cpatch%mmean_sensible_wc     (recc)         &
                                             + cpatch%mmean_sensible_wc     (donc)
         cpatch%mmean_vapor_wc        (recc) = cpatch%mmean_vapor_wc        (recc)         &
                                             + cpatch%mmean_vapor_wc        (donc)
         cpatch%mmean_intercepted_aw  (recc) = cpatch%mmean_intercepted_aw  (recc)         &
                                             + cpatch%mmean_intercepted_aw  (donc)
         cpatch%mmean_wshed_wg        (recc) = cpatch%mmean_wshed_wg        (recc)         &
                                             + cpatch%mmean_wshed_wg        (donc)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      We update temperature and liquid water fraction.  We check whether the     !
         ! heat capacity is non-zero.  If it is a normal number, use the standard thermo-  !
         ! dynamic library, otherwise, average temperature, this is probably a blend of    !
         ! tiny cohorts that couldn't be solved, or the wood is not solved.                !
         !---------------------------------------------------------------------------------!
         !------ Leaf. --------------------------------------------------------------------!
         if ( cpatch%mmean_leaf_hcap(recc) > 0. ) then
            !----- Update temperature and liquid fraction using standard thermodynamics. --!
            call uextcm2tl(cpatch%mmean_leaf_energy(recc),cpatch%mmean_leaf_water(recc)    &
                          ,cpatch%mmean_leaf_hcap  (recc),cpatch%mmean_leaf_temp (recc)    &
                          ,cpatch%mmean_leaf_fliq  (recc))                                                  
            !------------------------------------------------------------------------------!


            !----- Scale vapour pressure deficit using LAI. -------------------------------!
            cpatch%mmean_leaf_vpdef   (recc) = cpatch%mmean_leaf_vpdef      (recc) * rlai  &
                                             + cpatch%mmean_leaf_vpdef      (donc) * dlai
            !------------------------------------------------------------------------------!
         else
            !----- None of the cohorts has leaf biomass use nplant to scale them. ---------!
            cpatch%mmean_leaf_temp (recc) = cpatch%mmean_leaf_temp (recc) * rnplant        &
                                          + cpatch%mmean_leaf_temp (donc) * dnplant
            cpatch%mmean_leaf_fliq (recc) = 0.0
            cpatch%mmean_leaf_vpdef(recc) = cpatch%mmean_leaf_vpdef(recc) * rnplant        &
                                          + cpatch%mmean_leaf_vpdef(donc) * dnplant
            !------------------------------------------------------------------------------!
         end if                                                                            
         !------ Wood. --------------------------------------------------------------------!
         if ( cpatch%mmean_wood_hcap(recc) > 0. ) then                                     
            !----- Update temperature using the standard thermodynamics. ------------------!
            call uextcm2tl(cpatch%mmean_wood_energy(recc),cpatch%mmean_wood_water(recc)    &
                          ,cpatch%mmean_wood_hcap  (recc),cpatch%mmean_wood_temp (recc)    &
                          ,cpatch%mmean_wood_fliq  (recc))
            !------------------------------------------------------------------------------!
         else                                                                              
            !----- Wood temperature can't be found using uextcm2tl (singularity). ---------!
            cpatch%mmean_wood_temp(recc) = cpatch%mmean_wood_temp(recc) * rnplant          &
                                         + cpatch%mmean_wood_temp(donc) * dnplant
            cpatch%mmean_wood_fliq(recc) = 0.0
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Fuse mortality rates.  By definition, mortality rate M for is given by:      !
         !                                                                                 !
         !                                    ln(A) - ln(N)                                !
         !                               m = ---------------                               !
         !                                          dt                                     !
         !                                                                                 !
         ! where A is the population that was previously alive, and N is the population    !
         ! that survived the mortality.  The cohorts represent a group of individuals with !
         ! the same size and PFT, so they don't mix new recruits and old plants. There-    !
         ! fore, we can assume that N is actually nplant.  We don't know A, if the         !
         ! mortality rate is assumed constant during the interval dt, A = N * exp(m dt).   !
         !                                                                                 !
         ! For fusion we don't really care about dt, so any number will do as long as it   !
         ! is the same for both cohorts.  With these assumptions, the mortality rate for   !
         ! the fused cohort mf is:                                                         !
         !                                                                                 !
         !  mf   =  ln (Ad+Ar) - ln(Nd+Nr) = ln[Nd*exp(md) + Nr*exp(mr)] - ln[Nd+Nr]       !
         !                                                                                 !
         !             / Nd*exp(md) + Nr*exp(mr) \                                         !
         !  mf   =  ln |-------------------------|                                         !
         !             \        Nd + Nr          /                                         !
         !---------------------------------------------------------------------------------!
         do imty=1,n_mort
            exp_mort_donc = exp(max(lnexp_min,min(lnexp_max                                &
                                                 ,cpatch%mmean_mort_rate(imty,donc))))
            exp_mort_recc = exp(max(lnexp_min,min(lnexp_max                                &
                                                 ,cpatch%mmean_mort_rate(imty,recc))))
            cpatch%mmean_mort_rate(imty,recc) = log( exp_mort_recc * rnplant               &
                                                   + exp_mort_donc * dnplant )
         end do
         !------------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Mean diel.                                                                      !
      !------------------------------------------------------------------------------------!
      if (writing_dcyc .and. (.not. fuse_initial)) then
         !---------------------------------------------------------------------------------!
         !     First we merge the squares, as they require the means.                      !
         !---------------------------------------------------------------------------------!
         do t=1,ndcycle
            !----- Intensive. -------------------------------------------------------------!
            cpatch%qmsqu_gpp        (t,recc) = fuse_msqu(cpatch%qmean_gpp        (t,recc)  &
                                                        ,cpatch%qmsqu_gpp        (t,recc)  &
                                                        ,rnplant                           &
                                                        ,cpatch%qmean_gpp        (t,donc)  &
                                                        ,cpatch%qmsqu_gpp        (t,donc)  &
                                                        ,dnplant                           &
                                                        ,corr_cohort,.false.)
            cpatch%qmsqu_npp        (t,recc) = fuse_msqu(cpatch%qmean_npp        (t,recc)  &
                                                        ,cpatch%qmsqu_npp        (t,recc)  &
                                                        ,rnplant                           &
                                                        ,cpatch%qmean_npp        (t,donc)  &
                                                        ,cpatch%qmsqu_npp        (t,donc)  &
                                                        ,dnplant                           &
                                                        ,corr_cohort,.false.)
            cpatch%qmsqu_plresp     (t,recc) = fuse_msqu(cpatch%qmean_plresp     (t,recc)  &
                                                        ,cpatch%qmsqu_plresp     (t,recc)  &
                                                        ,rnplant                           &
                                                        ,cpatch%qmean_plresp     (t,donc)  &
                                                        ,cpatch%qmsqu_plresp     (t,donc)  &
                                                        ,dnplant                           &
                                                        ,corr_cohort,.false.)
            !----- Extensive variables. ---------------------------------------------------!
            cpatch%qmsqu_sensible_lc(t,recc) = fuse_msqu(cpatch%qmean_sensible_lc(t,recc)  &
                                                        ,cpatch%qmsqu_sensible_lc(t,recc)  &
                                                        ,rnplant                           &
                                                        ,cpatch%qmean_sensible_lc(t,donc)  &
                                                        ,cpatch%qmsqu_sensible_lc(t,donc)  &
                                                        ,dnplant                           &
                                                        ,corr_cohort,.true. )
            cpatch%qmsqu_vapor_lc   (t,recc) = fuse_msqu(cpatch%qmean_vapor_lc   (t,recc)  &
                                                        ,cpatch%qmsqu_vapor_lc   (t,recc)  &
                                                        ,rnplant                           &
                                                        ,cpatch%qmean_vapor_lc   (t,donc)  &
                                                        ,cpatch%qmsqu_vapor_lc   (t,donc)  &
                                                        ,dnplant                           &
                                                        ,corr_cohort,.true. )
            cpatch%qmsqu_transp     (t,recc) = fuse_msqu(cpatch%qmean_transp     (t,recc)  &
                                                        ,cpatch%qmsqu_transp     (t,recc)  &
                                                        ,rnplant                           &
                                                        ,cpatch%qmean_transp     (t,donc)  &
                                                        ,cpatch%qmsqu_transp     (t,donc)  &
                                                        ,dnplant                           &
                                                        ,corr_cohort,.true. )
            cpatch%qmsqu_sensible_wc(t,recc) = fuse_msqu(cpatch%qmean_sensible_wc(t,recc)  &
                                                        ,cpatch%qmsqu_sensible_wc(t,recc)  &
                                                        ,rnplant                           &
                                                        ,cpatch%qmean_sensible_wc(t,donc)  &
                                                        ,cpatch%qmsqu_sensible_wc(t,donc)  &
                                                        ,dnplant                           &
                                                        ,corr_cohort,.true. )
            cpatch%qmsqu_vapor_wc   (t,recc) = fuse_msqu(cpatch%qmean_vapor_wc   (t,recc)  &
                                                        ,cpatch%qmsqu_vapor_wc   (t,recc)  &
                                                        ,rnplant                           &
                                                        ,cpatch%qmean_vapor_wc   (t,donc)  &
                                                        ,cpatch%qmsqu_vapor_wc   (t,donc)  &
                                                        ,dnplant                           &
                                                        ,corr_cohort,.true. )
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Intensive variables, scaled by plant density.                                !
         !---------------------------------------------------------------------------------!
         cpatch%qmean_gpp             (:,recc) = cpatch%qmean_gpp             (:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_gpp             (:,donc)     &
                                               * dnplant
         cpatch%qmean_npp             (:,recc) = cpatch%qmean_npp             (:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_npp             (:,donc)     &
                                               * dnplant
         cpatch%qmean_leaf_resp       (:,recc) = cpatch%qmean_leaf_resp       (:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_leaf_resp       (:,donc)     &
                                               * dnplant
         cpatch%qmean_root_resp       (:,recc) = cpatch%qmean_root_resp       (:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_root_resp       (:,donc)     &
                                               * dnplant
         cpatch%qmean_leaf_growth_resp(:,recc) = cpatch%qmean_leaf_growth_resp(:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_leaf_growth_resp(:,donc)     &
                                               * dnplant
         cpatch%qmean_root_growth_resp(:,recc) = cpatch%qmean_root_growth_resp(:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_root_growth_resp(:,donc)     &
                                               * dnplant
         cpatch%qmean_sapa_growth_resp(:,recc) = cpatch%qmean_sapa_growth_resp(:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_sapa_growth_resp(:,donc)     &
                                               * dnplant
         cpatch%qmean_sapb_growth_resp(:,recc) = cpatch%qmean_sapb_growth_resp(:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_sapb_growth_resp(:,donc)     &
                                               * dnplant
         cpatch%qmean_storage_resp    (:,recc) = cpatch%qmean_storage_resp    (:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_storage_resp    (:,donc)     &
                                               * dnplant
         cpatch%qmean_plresp          (:,recc) = cpatch%qmean_plresp          (:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_plresp          (:,donc)     &
                                               * dnplant
         cpatch%qmean_light_level     (:,recc) = cpatch%qmean_light_level     (:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_light_level     (:,donc)     &
                                               * dnplant
         cpatch%qmean_light_level_beam(:,recc) = cpatch%qmean_light_level_beam(:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_light_level_beam(:,donc)     &
                                               * dnplant
         cpatch%qmean_light_level_diff(:,recc) = cpatch%qmean_light_level_diff(:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_light_level_diff(:,donc)     &
                                               * dnplant
         cpatch%qmean_par_level_beam  (:,recc) = cpatch%qmean_par_level_beam  (:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_par_level_beam  (:,donc)     &
                                               * dnplant
         cpatch%qmean_par_level_diffd (:,recc) = cpatch%qmean_par_level_diffd (:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_par_level_diffd (:,donc)     &
                                               * dnplant
         cpatch%qmean_par_level_diffu (:,recc) = cpatch%qmean_par_level_diffu (:,recc)     &
                                               * rnplant                                   &
                                               + cpatch%qmean_par_level_diffu (:,donc)     &
                                               * dnplant
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Intensive variables, scaled by LAI.                                          !
         !---------------------------------------------------------------------------------!
         cpatch%qmean_leaf_gsw      (:,recc) = cpatch%qmean_leaf_gsw      (:,recc) * rlai  &
                                             + cpatch%qmean_leaf_gsw      (:,donc) * dlai
         cpatch%qmean_leaf_gbw      (:,recc) = cpatch%qmean_leaf_gbw      (:,recc) * rlai  &
                                             + cpatch%qmean_leaf_gbw      (:,donc) * dlai
         cpatch%qmean_fs_open       (:,recc) = cpatch%qmean_fs_open       (:,recc) * rlai  &
                                             + cpatch%qmean_fs_open       (:,donc) * dlai
         cpatch%qmean_fsw           (:,recc) = cpatch%qmean_fsw           (:,recc) * rlai  &
                                             + cpatch%qmean_fsw           (:,donc) * dlai
         cpatch%qmean_fsn           (:,recc) = cpatch%qmean_fsn           (:,recc) * rlai  &
                                             + cpatch%qmean_fsn           (:,donc) * dlai
         cpatch%qmean_A_open        (:,recc) = cpatch%qmean_A_open        (:,recc) * rlai  &
                                             + cpatch%qmean_A_open        (:,donc) * dlai
         cpatch%qmean_A_closed      (:,recc) = cpatch%qmean_A_closed      (:,recc) * rlai  &
                                             + cpatch%qmean_A_closed      (:,donc) * dlai
         cpatch%qmean_A_net         (:,recc) = cpatch%qmean_A_net         (:,recc) * rlai  &
                                             + cpatch%qmean_A_net         (:,donc) * dlai
         cpatch%qmean_A_light       (:,recc) = cpatch%qmean_A_light       (:,recc) * rlai  &
                                             + cpatch%qmean_A_light       (:,donc) * dlai
         cpatch%qmean_A_rubp        (:,recc) = cpatch%qmean_A_rubp        (:,recc) * rlai  &
                                             + cpatch%qmean_A_rubp        (:,donc) * dlai
         cpatch%qmean_A_co2         (:,recc) = cpatch%qmean_A_co2         (:,recc) * rlai  &
                                             + cpatch%qmean_A_co2         (:,donc) * dlai
         cpatch%qmean_psi_open      (:,recc) = cpatch%qmean_psi_open      (:,recc) * rlai  &
                                             + cpatch%qmean_psi_open      (:,donc) * dlai
         cpatch%qmean_psi_closed    (:,recc) = cpatch%qmean_psi_closed    (:,recc) * rlai  &
                                             + cpatch%qmean_psi_closed    (:,donc) * dlai
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Intensive variables, scaled by WAI.                                          !
         !---------------------------------------------------------------------------------!
         cpatch%qmean_wood_gbw      (:,recc) = cpatch%qmean_wood_gbw      (:,recc) * rwai  &
                                             + cpatch%qmean_wood_gbw      (:,donc) * dwai
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Extensive variables.                                                         !
         !---------------------------------------------------------------------------------!
         cpatch%qmean_leaf_energy   (:,recc) = cpatch%qmean_leaf_energy   (:,recc)         &
                                             + cpatch%qmean_leaf_energy   (:,donc)
         cpatch%qmean_leaf_water    (:,recc) = cpatch%qmean_leaf_water    (:,recc)         &
                                             + cpatch%qmean_leaf_water    (:,donc)
         cpatch%qmean_leaf_hcap     (:,recc) = cpatch%qmean_leaf_hcap     (:,recc)         &
                                             + cpatch%qmean_leaf_hcap     (:,donc)
         cpatch%qmean_wood_energy   (:,recc) = cpatch%qmean_wood_energy   (:,recc)         &
                                             + cpatch%qmean_wood_energy   (:,donc)
         cpatch%qmean_wood_water    (:,recc) = cpatch%qmean_wood_water    (:,recc)         &
                                             + cpatch%qmean_wood_water    (:,donc)
         cpatch%qmean_wood_hcap     (:,recc) = cpatch%qmean_wood_hcap     (:,recc)         &
                                             + cpatch%qmean_wood_hcap     (:,donc)
         cpatch%qmean_water_supply  (:,recc) = cpatch%qmean_water_supply  (:,recc)         &
                                             + cpatch%qmean_water_supply  (:,donc)
         cpatch%qmean_par_l         (:,recc) = cpatch%qmean_par_l         (:,recc)         &
                                             + cpatch%qmean_par_l         (:,donc)
         cpatch%qmean_par_l_beam    (:,recc) = cpatch%qmean_par_l_beam    (:,recc)         &
                                             + cpatch%qmean_par_l_beam    (:,donc)
         cpatch%qmean_par_l_diff    (:,recc) = cpatch%qmean_par_l_diff    (:,recc)         &
                                             + cpatch%qmean_par_l_diff    (:,donc)

 

         cpatch%qmean_rshort_l      (:,recc) = cpatch%qmean_rshort_l      (:,recc)         &
                                             + cpatch%qmean_rshort_l      (:,donc)
         cpatch%qmean_rlong_l       (:,recc) = cpatch%qmean_rlong_l       (:,recc)         &
                                             + cpatch%qmean_rlong_l       (:,donc)
         cpatch%qmean_sensible_lc   (:,recc) = cpatch%qmean_sensible_lc   (:,recc)         &
                                             + cpatch%qmean_sensible_lc   (:,donc)
         cpatch%qmean_vapor_lc      (:,recc) = cpatch%qmean_vapor_lc      (:,recc)         &
                                             + cpatch%qmean_vapor_lc      (:,donc)
         cpatch%qmean_transp        (:,recc) = cpatch%qmean_transp        (:,recc)         &
                                             + cpatch%qmean_transp        (:,donc)
         cpatch%qmean_intercepted_al(:,recc) = cpatch%qmean_intercepted_al(:,recc)         &
                                             + cpatch%qmean_intercepted_al(:,donc)
         cpatch%qmean_wshed_lg      (:,recc) = cpatch%qmean_wshed_lg      (:,recc)         &
                                             + cpatch%qmean_wshed_lg      (:,donc)
         cpatch%qmean_rshort_w      (:,recc) = cpatch%qmean_rshort_w      (:,recc)         &
                                             + cpatch%qmean_rshort_w      (:,donc)
         cpatch%qmean_rlong_w       (:,recc) = cpatch%qmean_rlong_w       (:,recc)         &
                                             + cpatch%qmean_rlong_w       (:,donc)
         cpatch%qmean_rad_profile (:,:,recc) = cpatch%qmean_rad_profile (:,:,recc)         &
                                             + cpatch%qmean_rad_profile (:,:,donc)
         cpatch%qmean_sensible_wc   (:,recc) = cpatch%qmean_sensible_wc   (:,recc)         &
                                             + cpatch%qmean_sensible_wc   (:,donc)
         cpatch%qmean_vapor_wc      (:,recc) = cpatch%qmean_vapor_wc      (:,recc)         &
                                             + cpatch%qmean_vapor_wc      (:,donc)
         cpatch%qmean_intercepted_aw(:,recc) = cpatch%qmean_intercepted_aw(:,recc)         &
                                             + cpatch%qmean_intercepted_aw(:,donc)
         cpatch%qmean_wshed_wg      (:,recc) = cpatch%qmean_wshed_wg      (:,recc)         &
                                             + cpatch%qmean_wshed_wg      (:,donc)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      We update temperature and liquid water fraction.  We check whether the     !
         ! heat capacity is non-zero.  If it is a normal number, use the standard thermo-  !
         ! dynamic library, otherwise, average temperature, this is probably a blend of    !
         ! tiny cohorts that couldn't be solved, or the wood is not solved.                !
         !---------------------------------------------------------------------------------!
         do t=1,ndcycle
            !------ Leaf. -----------------------------------------------------------------!
            if ( cpatch%qmean_leaf_hcap(t,recc) > 0. ) then
               !----- Update temperature and liquid fraction. -----------------------------!
               call uextcm2tl(cpatch%qmean_leaf_energy(t,recc)                             &
                             ,cpatch%qmean_leaf_water (t,recc)                             &
                             ,cpatch%qmean_leaf_hcap  (t,recc)                             &
                             ,cpatch%qmean_leaf_temp  (t,recc)                             &
                             ,cpatch%qmean_leaf_fliq  (t,recc))
               !---------------------------------------------------------------------------!


               !----- Scale vapour pressure deficit using LAI. ----------------------------!
               cpatch%qmean_leaf_vpdef (t,recc) = cpatch%qmean_leaf_vpdef  (t,recc) * rlai &
                                                + cpatch%qmean_leaf_vpdef  (t,donc) * dlai
               !---------------------------------------------------------------------------!
            else
               !----- None of the cohorts has leaf biomass use nplant to scale them. ------!
               cpatch%qmean_leaf_temp (t,recc) = cpatch%qmean_leaf_temp (t,recc) * rnplant &
                                               + cpatch%qmean_leaf_temp (t,donc) * dnplant
               cpatch%qmean_leaf_fliq (t,recc) = 0.0
               cpatch%qmean_leaf_vpdef(t,recc) = cpatch%qmean_leaf_vpdef(t,recc) * rnplant &
                                               + cpatch%qmean_leaf_vpdef(t,donc) * dnplant
               !---------------------------------------------------------------------------!
            end if
            !------ Wood. -----------------------------------------------------------------!
            if ( cpatch%qmean_wood_hcap(t,recc) > 0. ) then
               !----- Update temperature and liquid fraction. -----------------------------!
               call uextcm2tl(cpatch%qmean_wood_energy(t,recc)                             &
                             ,cpatch%qmean_wood_water (t,recc)                             &
                             ,cpatch%qmean_wood_hcap  (t,recc)                             &
                             ,cpatch%qmean_wood_temp  (t,recc)                             &
                             ,cpatch%qmean_wood_fliq  (t,recc))
               !---------------------------------------------------------------------------!
            else
               !----- None of the cohorts has leaf biomass use nplant to scale them. ------!
               cpatch%qmean_wood_temp (t,recc) = cpatch%qmean_wood_temp (t,recc) * rnplant &
                                               + cpatch%qmean_wood_temp (t,donc) * dnplant
               cpatch%qmean_wood_fliq (t,recc) = 0.0
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Lastly, we update nplant, LAI, WAI, and crown area.  They may be used as       !
      ! weighting factors so we leave them to the end.                                     !
      !------------------------------------------------------------------------------------!
      cpatch%nplant     (recc) = cpatch%nplant(recc) + cpatch%nplant(donc)
      cpatch%lai        (recc) = cpatch%lai   (recc) + cpatch%lai   (donc)
      cpatch%wai        (recc) = cpatch%wai   (recc) + cpatch%wai   (donc)
      !----- Make sure that crown area is bounded. ----------------------------------------!
      cpatch%crown_area (recc) = min(1.,cpatch%crown_area(recc)  + cpatch%crown_area(donc))
      !------------------------------------------------------------------------------------!

      return
   end subroutine fuse_2_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will sort the patches by age (1st = oldest, last = youngest.)       !
   !---------------------------------------------------------------------------------------!
   subroutine sort_patches(csite)

      use ed_state_vars, only  :  sitetype   ! ! Structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype), target   :: csite      ! Current site, that will have patches sorted.
      !----- Local variables --------------------------------------------------------------!
      type(sitetype), pointer  :: tempsite   ! Structure to temporarily host the sorted site
      integer                  :: ipa        ! Counters
      integer                  :: oldpa      ! Index of oldest patch
      logical                  :: sorted     ! Flag: the site is already sorted
      !------------------------------------------------------------------------------------!
      
      !----- No need to sort a site with a single patch. ----------------------------------!
      if (csite%npatches < 2) return

      !------------------------------------------------------------------------------------!
      !     Check whether this site is already sorted.   We don't want to do the entire    !
      ! deallocating/copying/allocating thing if it's not needed as this takes up too much !
      ! time.                                                                              !
      !------------------------------------------------------------------------------------!
      sorted = .true.
      sortcheck: do ipa=1,csite%npatches-1
         sorted = csite%age(ipa) >= csite%age(ipa+1)
         if (.not. sorted) exit sortcheck
      end do sortcheck
      if (sorted) return
      !------------------------------------------------------------------------------------!



      !----- Assign a scratch patch. ------------------------------------------------------!
      nullify (tempsite)
      allocate(tempsite)
      call allocate_sitetype(tempsite,csite%npatches)
      
      ipa = 0
      !---- Loop until all patches were sorted. -------------------------------------------!
      do while (ipa < csite%npatches)
         ipa = ipa + 1
      
         !----- Find the oldest site. -----------------------------------------------------!
         oldpa = maxloc(csite%age,dim=1)
         
         !----- Copy to patch the scratch structure. --------------------------------------!
         call copy_sitetype(csite,tempsite,oldpa,oldpa,ipa,ipa)
         
         !----- Put a non-sense age so this patch will never "win" again. -----------------!
         csite%age(oldpa) = -huge(1.)
      end do

      !------ Reset the actual patch, and re-allocate it. ---------------------------------!
      call deallocate_sitetype(csite)
      call allocate_sitetype  (csite,tempsite%npatches)

      !------ Copy the scratch patch to the regular one. ----------------------------------!
      call copy_sitetype(tempsite,csite,1,tempsite%npatches,1,tempsite%npatches)
      call deallocate_sitetype(tempsite)
      deallocate(tempsite)

      return

   end subroutine sort_patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will perform patch fusion based on some similarity criteria to      !
   ! determine whether they can be fused with no significant loss of information. The user !
   ! is welcome to set up a benchmark, but they should be aware that no miracles will      !
   ! happen here. If there are more very distinct patches than maxpatch, then the user     !
   ! will need to live with that and accept life is not always fair with those with        !
   ! limited computational resources.                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_patches(cgrid,ifm,fuse_initial)
      use ed_state_vars       , only : edtype              & ! structure
                                     , polygontype         & ! structure
                                     , sitetype            & ! structure
                                     , patchtype           ! ! structure
      use fusion_fission_coms , only : ff_nhgt             & ! intent(in)
                                     , niter_patfus        & ! intent(in)
                                     , dark_cumlai_min     & ! intent(in)
                                     , dark_cumlai_max     & ! intent(in)
                                     , dark_cumlai_mult    & ! intent(in)
                                     , sunny_cumlai_min    & ! intent(in)
                                     , sunny_cumlai_max    & ! intent(in)
                                     , sunny_cumlai_mult   & ! intent(in)
                                     , print_fuse_details  & ! intent(in)
                                     , light_toler_min     & ! intent(in)
                                     , light_toler_max     & ! intent(in)
                                     , light_toler_mult    & ! intent(in)
                                     , min_oldgrowth       & ! intent(in)
                                     , fuse_prefix         ! ! intent(in)
      use ed_max_dims         , only : n_pft               & ! intent(in)
                                     , str_len             ! ! intent(in)
      use mem_polygons        , only : maxpatch            & ! intent(in)
                                     , maxcohort           ! ! intent(in)
      use ed_node_coms        , only : mynum               ! ! intent(in)
      use ed_misc_coms        , only : current_time        ! ! intent(in)
      use grid_coms           , only : nzg                 & ! intent(in)
                                     , nzs                 ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)          , target      :: cgrid           ! Current grid
      integer               , intent(in)  :: ifm             ! Current grid index
      logical               , intent(in)  :: fuse_initial    ! Called from initialisation?
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)     , pointer     :: cpoly           ! Current polygon
      type(polygontype)     , pointer     :: jpoly           ! Current polygon
      type(sitetype)        , pointer     :: csite           ! Current site
      type(patchtype)       , pointer     :: cpatch          ! Current patch
      type(patchtype)       , pointer     :: donpatch        ! Donor patch
      type(patchtype)       , pointer     :: recpatch        ! Receptor patch
      type(sitetype)        , pointer     :: tempsite        ! Temporary site
      logical, dimension(:) , allocatable :: fuse_table      ! Flag: this will remain.
      character(len=str_len)              :: fuse_fout       ! Filename for detailed output
      integer                             :: ipy             ! Counters
      integer                             :: isi             ! Counters
      integer                             :: jpy             ! Counters
      integer                             :: jsi             ! Counters
      integer                             :: ipa             ! Counters
      integer                             :: ico             ! Counters
      integer                             :: donp            ! Counters
      integer                             :: recp            ! Counters
      integer                             :: rec_lu          ! Land use of receptor patch
      integer                             :: don_lu          ! Land use of donor patch
      integer                             :: ihgt            ! Counters
      integer                             :: ifus            ! Counters
      integer                             :: npatches_new    ! New # of patches
      integer                             :: npatches_old    ! Old # of patches
      integer                             :: npatches_orig   ! Original # of patches
      logical                             :: fuse_flag       ! Flag: fusion will happen
      logical                             :: recp_found      ! Found a receptor candidate
      logical                             :: sunny_donp      ! Donor patch bin too sunny
      logical                             :: sunny_recp      ! Receptor patch bin too sunny
      logical                             :: dark_donp       ! Donor patch bin too small
      logical                             :: dark_recp       ! Receptor patch bin too small
      logical                             :: same_age        ! Patches with same age
      logical                             :: old_or_same_lu  ! Old patches or the same LU.
      real                                :: diff            ! Absolute difference in prof.
      real                                :: refv            ! Reference value of bin
      real                                :: norm            ! Normalised difference
      real                                :: llevel_donp     ! Light level of donor patch
      real                                :: llevel_recp     ! Light level of rec.  patch
      real                                :: sunny_toler     ! Light layer tolerance.
      real                                :: dark_lai80      ! Minimum dark layer.
      real                                :: dark_toler      ! Dark layer tolerance.
      real                                :: light_toler     ! Light level Relative toler.
      real                                :: old_area        ! For area conservation check
      real                                :: new_area        ! For area conservation check
      real                                :: now_area        ! Area for LAI80
      real                                :: old_lai_tot     ! Old total LAI
      real                                :: old_nplant_tot  ! Old total nplant
      real                                :: new_lai_tot     ! New total LAI
      real                                :: new_nplant_tot  ! New total nplant
      real                                :: elim_nplant     ! Elim. nplant during 1 fusion
      real                                :: elim_lai        ! Elim. LAI during 1 fusion
      real                                :: elim_nplant_tot ! Total eliminated nplant
      real                                :: elim_lai_tot    ! Elim. eliminated LAI
      real                                :: cumlai_recp     ! Cumulative LAI (receptor)
      real                                :: cumlai_donp     ! Cumulative LAI (donor)
      integer                             :: tot_npolygons   ! Total # of polygons
      integer                             :: tot_nsites      ! Total # of sites
      integer                             :: tot_npatches    ! Total # of patches
      integer                             :: tot_ncohorts    ! Total # of cohorts
      !----- Locally saved variables. --------------------------------------------------------!
      logical                   , save    :: first_time = .true.
      logical                   , save    :: dont_force_fuse
      logical                   , save    :: force_fuse
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     First time here.  Delete all files.                                            !
      !------------------------------------------------------------------------------------!
      if (first_time) then
         force_fuse      = abs(maxpatch) == 1
         dont_force_fuse = .not. force_fuse
         if (print_fuse_details) then
            do jpy = 1, cgrid%npolygons
               jpoly => cgrid%polygon(jpy)
               do jsi = 1, jpoly%nsites
                  write (fuse_fout,fmt='(a,2(a,i4.4),a)')                                  &
                        trim(fuse_prefix),'polygon_',jpy,'_site_',jsi,'.txt'
                  open (unit=72,file=trim(fuse_fout),status='replace',action='write')
                  write(unit=72,fmt='(a)')       '----------------------------------------'
                  write(unit=72,fmt='(a)')       ' Patch Fusion log for: '
                  write(unit=72,fmt='(a,1x,i5)') ' POLYGON: ',jpy 
                  write(unit=72,fmt='(a,1x,i5)') ' SITE:    ',jsi 
                  write(unit=72,fmt='(a)')       '----------------------------------------'
                  write(unit=72,fmt='(a)')       ' '
                  close(unit=72,status='keep')
               end do
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Return if maxpatch is 0, this is a flag for no patch fusion.                   !
      !------------------------------------------------------------------------------------!
      if (maxpatch == 0) return
      !------------------------------------------------------------------------------------!




      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         
         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            write(fuse_fout,fmt='(a,2(a,i4.4),a)')                                         &
                     trim(fuse_prefix),'polygon_',ipy,'_site_',isi,'.txt'

            if (print_fuse_details) then
               open (unit=72,file=trim(fuse_fout),status='old',action='write'              &
                                                 ,position='append')

               write(unit=72,fmt='(2(a,i2.2),a,i4.4)')    ' - Date: ',current_time%month   &
                                                                 ,'/',current_time%date    &
                                                                 ,'/',current_time%year
               write(unit=72,fmt='(a,1x,i6)') '   + Initial number of patches: '           &
                                             ,csite%npatches
               write(unit=72,fmt='(a)')       '   + Looking for empty patches: '
               close(unit=72,status='keep')
            end if

            !----- Skip this site if it contains only one patch... ------------------------!
            if (csite%npatches < 2) cycle siteloop

            !----- Save original number of patches. ---------------------------------------!
            npatches_orig = csite%npatches

            !----- Allocate the swapper patches in the site type. -------------------------!
            nullify(tempsite)
            allocate(tempsite)
            call allocate_sitetype(tempsite, csite%npatches)
            allocate(fuse_table(csite%npatches))

            !------------------------------------------------------------------------------!
            !     Allocate the fusion flag vector, and set all elements to .true., which   !
            ! means that every patch can be fused.  As soon as the patch is fused, we will !
            ! switch the flag to false.                                                    !
            !------------------------------------------------------------------------------!
            fuse_table(:) = .true.
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find the original number of plants, LAI, and area, which will be used    !
            ! for sanity check.  We also compute the pft and size profile for all patches, !
            ! which will be used for the fusion criterion.                                 !
            !------------------------------------------------------------------------------!
            old_nplant_tot = 0.
            old_lai_tot    = 0.
            old_area       = 0.
            do ipa = 1,csite%npatches
               call patch_pft_size_profile(csite,ipa)

               old_area  = old_area + csite%area(ipa)
               cpatch => csite%patch(ipa)
               do ico = 1, cpatch%ncohorts
                  old_nplant_tot = old_nplant_tot + cpatch%nplant(ico) * csite%area(ipa)
                  old_lai_tot    = old_lai_tot    + cpatch%lai(ico)    * csite%area(ipa)
               end do
            end do
            !------------------------------------------------------------------------------!



            !----- Initialise the total eliminated nplant and LAI to zero. ----------------!
            elim_nplant_tot = 0.
            elim_lai_tot    = 0.
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     This loop will check whether there is at least two patches with exact    !
            ! same age.  This is common in initialisation with site-level measurements.    !
            ! In this case we want to merge any patch, regardless on their "age position", !
            ! which doesn't mean anything in this case.                                    !
            !------------------------------------------------------------------------------!
            same_age=.false.
            donloop_check: do donp=csite%npatches,2,-1
               recloop_check: do recp=donp-1,1,-1
                  same_age = csite%age(recp)       == csite%age(donp) .and.                &
                             csite%dist_type(recp) == csite%dist_type(donp)
                  if (print_fuse_details) then
                        open (unit=72,file=trim(fuse_fout),status='old',action='write'     &
                                              ,position='append')
                        write(unit=72,fmt='(a,1x,l1)') '     * same_age is ',same_age
                        close(unit=72,status='keep')
                  end if
                  !----- At least two patches have the same age. --------------------------!
                  if (same_age) exit donloop_check
              end do recloop_check
            end do donloop_check
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            donloope: do donp=csite%npatches,2,-1
               donpatch => csite%patch(donp)
               don_lu = csite%dist_type(donp)
               
               !----- If patch is not empty, or has already been fused, move on. ----------!
               if ( (.not. fuse_table(donp)) .or.                                          &
                    ( dont_force_fuse .and.  donpatch%ncohorts > 0) ) then
                  cycle donloope
               end if
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !     If we reach this point, it means that the donor patch is empty and    !
               ! hasn't been fused yet: look for an older empty patch and merge them.      !
               !---------------------------------------------------------------------------!
               if (print_fuse_details) then
                  open (unit=72,file=trim(fuse_fout),status='old',action='write'           &
                                                    ,position='append')
                  write(unit=72,fmt='(a,i6,a)') '     * Patch ',donp,' is empty.'
                  close(unit=72,status='keep')
               end if
               recloope: do recp=donp-1,1,-1
                  recpatch => csite%patch(recp)
                  rec_lu = csite%dist_type(recp)

                  !------------------------------------------------------------------------!
                  !     Set this flag that checks whether the patches have the same        !
                  ! disturbance type or are too old so we don't need to distinguish them.  !
                  !------------------------------------------------------------------------!
                  old_or_same_lu =   don_lu == rec_lu                         .or.         &
                                   ( csite%age(donp) >= min_oldgrowth(don_lu) .and.        &
                                     csite%age(recp) >= min_oldgrowth(rec_lu) )
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !     Skip the patch if it isn't empty, or it has already been fused, or !
                  ! if the donor and receptor have different disturbance types.            !
                  !------------------------------------------------------------------------!
                  if ( (.not. fuse_table(recp))                       .or.                 &
                       ( dont_force_fuse                              .and.                &
                         ( recpatch%ncohorts > 0 .or. (.not. old_or_same_lu) ) ) ) then
                     cycle recloope
                  end if
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !     Skip the patch if they don't have the same disturbance type and    !
                  ! are not too old.                                                       !
                  !------------------------------------------------------------------------!
                  if (.not. old_or_same_lu) cycle recloope
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !     Take an average of the patch properties of donpatch and recpatch,  !
                  ! and assign the average recpatch.                                       !
                  !------------------------------------------------------------------------!
                  call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)%prss          &
                                     ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)               &
                                     ,cpoly%green_leaf_factor(:,isi),fuse_initial          &
                                     ,elim_nplant,elim_lai)
                  !------------------------------------------------------------------------!


                  !----- Record the fusion if requested by the user. ----------------------!
                  if (print_fuse_details) then
                     open (unit=72,file=trim(fuse_fout),status='old',action='write'        &
                                                        ,position='append')
                     write(unit=72,fmt='(2(a,i6),a)') '     * Patches ',donp,' and ',recp  &
                                                     ,' were fused.'
                     close(unit=72,status='keep')
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Update total eliminated nplant and LAI.  This is actually not       !
                  ! necessary in this loop as both patches are empty, but we do it anyway  !
                  ! just to be consistent.                                                 !
                  !------------------------------------------------------------------------!
                  elim_nplant_tot = elim_nplant_tot + elim_nplant
                  elim_lai_tot    = elim_lai_tot    + elim_lai
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Recalculate the pft size profile for the averaged patch at recp.   !
                  ! Again, this is not really necessary as the receptor patch is empty,    !
                  ! but just to be consistent...                                           !
                  !------------------------------------------------------------------------!
                  call patch_pft_size_profile(csite,recp)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     The patch at index donp will be eliminated and should not be       !
                  ! checked for fusion again; we switch the fuse_table flag to .false. so  !
                  ! next time we reach this patch we will skip it.                         !
                  !------------------------------------------------------------------------!
                  fuse_table(donp) = .false.
                  !------------------------------------------------------------------------!

                  !------ We are done with donp, so we quit the recp loop. ----------------!
                  exit recloope

               end do recloope
            end do donloope
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Second loop.  Here we will fuse all patches with the same age and        !
            ! disturbance type.                                                            !
            !------------------------------------------------------------------------------!

            if (same_age) then
               !----- Start with no multiplication factor. --------------------------------!
               dark_toler   = dark_cumlai_max
               sunny_toler  = sunny_cumlai_min
               light_toler  = light_toler_min

               mainfuseloopa: do ifus=0,niter_patfus

                  npatches_old = count(fuse_table)
                  npatches_new = npatches_old

                  !------------------------------------------------------------------------!
                  !    Inform that the upcoming fusions are going to be with populated     !
                  ! patches, and record the tolerance used.                                !
                  !------------------------------------------------------------------------!
                  if (print_fuse_details) then
                     open (unit=72,file=trim(fuse_fout),status='old',action='write'        &
                                                        ,position='append')
                     write(unit=72,fmt='(a,1x,3(a,1x,es9.2,1x))')                          &
                                  '   + Looking for similar populated patches of same age' &
                                 ,' - Sunny Tolerance =',sunny_toler                       &
                                 ,' - Dark Tolerance  =',dark_toler                        &
                                 ,' - Rel. Tolerance  =',light_toler
                     close(unit=72,status='keep')
                  end if
                  !------------------------------------------------------------------------!

                  donloopa: do donp=csite%npatches,2,-1
                     donpatch => csite%patch(donp)
                     don_lu = csite%dist_type(donp)
                     
                     !----- If patch is not empty, or has already been fused, move on. ----!
                     if ( (.not. fuse_table(donp)) .or.                                    &
                          ( dont_force_fuse .and. donpatch%ncohorts == 0 ) ) then
                        cycle donloopa
                     end if
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     If we reach this point, it means that the donor patch is        !
                     ! populated and hasn't been fused yet: look for other patches with    !
                     ! the same age and merge them.                                        !
                     !---------------------------------------------------------------------!
                     if (print_fuse_details) then
                        open (unit=72,file=trim(fuse_fout),status='old',action='write'     &
                                                          ,position='append')
                        write(unit=72,fmt='(a,i6,a)') '     * Patch ',donp,' is populated.'
                        close(unit=72,status='keep')
                     end if
                     recloopa: do recp=donp-1,1,-1
                        recpatch => csite%patch(recp)
                        rec_lu = csite%dist_type(recp)

                        !------------------------------------------------------------------!
                        !     Skip the patch if it isn't empty, or it has already been     !
                        ! fused, or if the donor and receptor have different disturbance   !
                        ! types.                                                           !
                        !------------------------------------------------------------------!
                        if ( (.not. fuse_table(recp))                           .or.       &
                             ( dont_force_fuse                                 .and.       &
                               ( recpatch%ncohorts == 0                         .or.       &
                                 csite%dist_type(donp) /= csite%dist_type(recp) .or.       &
                                 csite%age(donp)       /=  csite%age(recp)           ) ) ) &
                        then
                           cycle recloopa
                        end if
                        !------------------------------------------------------------------!

                        !------------------------------------------------------------------!
                        !     Find the LAI that corresponds to 80% of the maximum LAI, to  !
                        ! avoid relaxing too much for forests.                             !
                        !------------------------------------------------------------------!
                        dark_lai80 = 0.40 * ( sum(csite%cumlai_profile(:,1,recp))          &
                                            + sum(csite%cumlai_profile(:,1,donp)) )


                        !------------------------------------------------------------------!
                        !     Compare the size profile for each PFT.  Here we compare the  !
                        ! maximum LAI for each PFT and height bin.  We switched the        !
                        ! classes from DBH to height because different PFTs may have       !
                        ! different heights for a given DBH, so we want to make sure the   !
                        ! light profile is right.                                          !
                        !------------------------------------------------------------------!
                        hgtloopa: do ihgt=1,ff_nhgt
                           cumlai_recp = sum(csite%cumlai_profile(:,ihgt,recp))
                           cumlai_donp = sum(csite%cumlai_profile(:,ihgt,donp))

                           !---------------------------------------------------------------!
                           !    Check whether these bins contain some LAI.  We don't need  !
                           ! to check the cohorts if the understorey is too dark, so once  !
                           ! both patches becomes very dark (very high LAI), we stop       !
                           ! checking the profiles.                                        !
                           !---------------------------------------------------------------!
                           dark_donp = cumlai_donp > dark_toler
                           dark_recp = cumlai_recp > dark_toler
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           if (dark_recp .and. dark_donp) then

                              if (print_fuse_details) then
                                 open  (unit=72,file=trim(fuse_fout),status='old'          &
                                               ,action='write',position='append')
                                 write (unit=72,fmt='(1(a,1x,i6,1x),4(a,1x,es9.2,1x)'//    &
                                                    ',2(a,1x,l1,1x))')                     &
                                    '       * IHGT =',ihgt                                 &
                                            ,'CUMLAI_RECP =',cumlai_recp                   &
                                            ,'CUMLAI_DONP =',cumlai_donp                   &
                                            ,'DARK_TOLER =',dark_toler                     &
                                            ,'DARK_LAI80 =',dark_lai80                     &
                                            ,'DARK_RECP =',dark_recp                       &
                                            ,'DARK_DONP =',dark_donp
                                 close (unit=72,status='keep')
                              end if
                              cycle hgtloopa
                           end if
                           !---------------------------------------------------------------!



                           
                           !---------------------------------------------------------------!
                           !    Check whether these bins contain some LAI.  Bins that have !
                           ! tiny cumulative LAI may differ by a lot in the relative       !
                           ! scale, but the actual value is so small that we don't really  !
                           ! care whether they are relatively different.                   !
                           !---------------------------------------------------------------!
                           sunny_donp = cumlai_donp <= sunny_toler
                           sunny_recp = cumlai_recp <= sunny_toler
                           !---------------------------------------------------------------!





                           !---------------------------------------------------------------!
                           !    If both patches have little or no biomass in this bin,     !
                           ! don't even bother checking the difference.                    !
                           !---------------------------------------------------------------!
                           if (sunny_donp .and. sunny_recp) then
                              if (print_fuse_details) then
                                 open  (unit=72,file=trim(fuse_fout),status='old'          &
                                               ,action='write',position='append')
                                 write (unit=72,fmt='(1(a,1x,i6,1x),3(a,1x,es9.2,1x)'//    &
                                                    ',2(a,1x,l1,1x))')                     &
                                    '       * IHGT=',ihgt                                  &
                                            ,'CUMLAI_RECP =',cumlai_recp                   &
                                            ,'CUMLAI_DONP =',cumlai_donp                   &
                                            ,'SUNNY_TOLER =',sunny_toler                   &
                                            ,'SUNNY_RECP =',sunny_recp                     &
                                            ,'SUNNY_RECP =',sunny_donp
                                 close (unit=72,status='keep')
                              end if
                              cycle hgtloopa
                           end if
                           !---------------------------------------------------------------!


                           !---------------------------------------------------------------!
                           !    Find the normalised difference in the density of this PFT  !
                           ! and size.  If one of the patches is missing any member of the !
                           ! profile the norm will be set to 2.0, which is the highest     !
                           ! value that the norm can be.                                   !
                           !---------------------------------------------------------------!
                           llevel_donp = exp(- 0.5 * cumlai_donp)
                           llevel_recp = exp(- 0.5 * cumlai_recp)
                           
                           diff = abs(llevel_donp - llevel_recp )
                           refv =    (llevel_donp + llevel_recp ) * 0.5
                           norm = diff / refv
                           fuse_flag = norm <= light_toler
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           if (print_fuse_details) then
                              open  (unit=72,file=trim(fuse_fout),status='old'             &
                                            ,action='write',position='append')
                              write (unit=72,fmt='(1(a,1x,i6,1x),7(a,1x,es9.2,1x)'//       &
                                                 ',1(a,1x,l1,1x))')                        &
                                 '       * IHGT=',ihgt                                     &
                                ,'CLAI_RECP =',cumlai_recp,'CLAI_DONP =',cumlai_donp       &
                                ,'LL_RECP =',llevel_recp,'LL_DONP =',llevel_donp           &
                                ,'DIFF =',diff,'REFV =',refv,'NORM =',norm                 &
                                ,'FUSE_FLAG =',fuse_flag
                              close (unit=72,status='keep')
                           end if
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           !     If fuse_flag is false, the patches aren't similar, move   !
                           ! to the next donor patch.                                      !
                           !---------------------------------------------------------------!
                           if (dont_force_fuse .and. (.not. fuse_flag)) cycle recloopa
                           !---------------------------------------------------------------!
                        end do hgtloopa
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Take an average of the patch properties of donpatch and      !
                        ! recpatch, and assign the average recpatch.                       !
                        !------------------------------------------------------------------!
                        call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)%prss    &
                                           ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)         &
                                           ,cpoly%green_leaf_factor(:,isi),fuse_initial    &
                                           ,elim_nplant,elim_lai)
                        !------------------------------------------------------------------!


                        !----- Record the fusion if requested by the user. ----------------!
                        if (print_fuse_details) then
                           open (unit=72,file=trim(fuse_fout),status='old',action='write'  &
                                                              ,position='append')
                           write(unit=72,fmt='(2(a,i6),a)') '     * Patches ',donp,' and ' &
                                                                   ,recp,' were fused.'
                           close(unit=72,status='keep')
                        end if
                        !------------------------------------------------------------------!


                        !------------------------------------------------------------------!
                        !    Update total eliminated nplant and LAI.  This is actually not !
                        ! necessary in this loop as both patches are empty, but we do it   !
                        ! anyway just to be consistent.                                    !
                        !------------------------------------------------------------------!
                        elim_nplant_tot = elim_nplant_tot + elim_nplant
                        elim_lai_tot    = elim_lai_tot    + elim_lai
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Recalculate the pft size profile for the averaged patch at   !
                        ! recp.  Again, this is not really necessary as the receptor patch !
                        ! is empty, but just to be consistent...                           !
                        !------------------------------------------------------------------!
                        call patch_pft_size_profile(csite,recp)
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     The patch at index donp will be eliminated and should not be !
                        ! checked for fusion again; we switch the fuse_table flag to       !
                        ! .false. so next time we reach this patch we will skip it.        !
                        !------------------------------------------------------------------!
                        fuse_table(donp) = .false.
                        !------------------------------------------------------------------!
                        !------------------------------------------------------------------!
                        !     Update the number of valid patches.                          !
                        !------------------------------------------------------------------!
                        npatches_new = npatches_new - 1
                        !------------------------------------------------------------------!

                        !------ We are done with donp, so we quit the recp loop. ----------!
                        exit recloopa
                     end do recloopa
                  end do donloopa

                  !------------------------------------------------------------------------!
                  !      Check how many patches are valid.  If the total number of patches !
                  ! is less than the target, or if we have reached the maximum tolerance   !
                  ! and the patch fusion still can't find similar patches, we quit the     !
                  ! fusion loop.                                                           !
                  !------------------------------------------------------------------------!
                  if (npatches_new <= abs(maxpatch)) exit mainfuseloopa
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Find the LAI that corresponds to 80% of the maximum LAI, to        !
                  ! avoid relaxing too much for forests.                                   !
                  !------------------------------------------------------------------------!
                  dark_lai80 = 0.
                  now_area   = 0.
                  do ipa=1,csite%npatches
                     if (fuse_table(ipa)) then
                        now_area   = now_area + csite%area(ipa)
                        dark_lai80 = dark_lai80                                            &
                                   + 0.80 * sum(csite%cumlai_profile(:,1,ipa))             &
                                   * csite%area(ipa)
                     end if
                  end do
                  if (now_area > 0.) dark_lai80 = dark_lai80 / now_area
                  !------------------------------------------------------------------------!



                  !----- Increment tolerance ----------------------------------------------!
                  sunny_toler =     sunny_toler * sunny_cumlai_mult
                  dark_toler  = max( dark_toler * dark_cumlai_mult , dark_lai80 )
                  light_toler =     light_toler * light_toler_mult
                  !------------------------------------------------------------------------!
               end do mainfuseloopa
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Third patch loop. Now that all empty patches have been fused, we will    !
            ! look for populated patches that have similar size and PFT structure, using   !
            ! the following algorithm:                                                     !
            !                                                                              !
            ! 1. Loop from the youngest to oldest patch;                                   !
            ! 2. Find next older patch with same dist_type;                                !
            ! 3. Check whether fusion criterion is met, and If so, fuse them.              !
            ! 4. After all the fusion, check how many patches we still have:               !
            !    If it is less than maxpatch, we quit, otherwise, we relax the tolerance a !
            !    little and try fusing more.  Notice that we will always try fusing        !
            !    patches at least once, even when the original number is less than         !
            !    maxpatch.                                                                 !
            !------------------------------------------------------------------------------!
            !----- Start with no multiplication factor. -----------------------------------!
            dark_toler   = dark_cumlai_max
            sunny_toler  = sunny_cumlai_min
            light_toler  = light_toler_min

            mainfuseloop: do ifus=0,niter_patfus

               npatches_old = count(fuse_table)
               npatches_new = npatches_old

               !---------------------------------------------------------------------------!
               !    Inform that the upcoming fusions are going to be with populated        !
               ! patches, and record the tolerance used.                                   !
               !---------------------------------------------------------------------------!
               if (print_fuse_details) then
                  open (unit=72,file=trim(fuse_fout),status='old',action='write'           &
                                                     ,position='append')
                  write(unit=72,fmt='(a,1x,3(a,1x,es9.2,1x))')                             &
                                              '   + Looking for similar populated patches' &
                                             ,' - Sunny Tolerance =',sunny_toler           &
                                             ,' - Dark Tolerance  =',dark_toler            &
                                             ,' - Rel. Tolerance  =',light_toler
                  close(unit=72,status='keep')
               end if
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !     Loop from youngest to the second oldest patch.                        !
               !---------------------------------------------------------------------------!
               donloopp: do donp = csite%npatches,2,-1
                  donpatch => csite%patch(donp)
                  don_lu = csite%dist_type(donp)

                  !------------------------------------------------------------------------!
                  !     If this is an empty patch, or has already been merged, we skip it. !
                  !------------------------------------------------------------------------!
                  if ((.not. fuse_table(donp))) then
                     if (print_fuse_details) then
                        open  (unit=72,file=trim(fuse_fout),status='old',action='write'    &
                                                           ,position='append')
                        write (unit=72,fmt='(a,1x,i6,1x,a)') '     - DONP:',donp           &
                                                            ,'has been already fused...'
                        close (unit=72,status='keep')
                     end if
                     cycle donloopp
                  end if
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !      If we have reached this place, the donor patch can be fused.  Now !
                  ! look for the next oldest patch that has the same disturbance type.  In !
                  ! case we can't find such patch, we will move to the next donor          !
                  ! candidate.                                                             !
                  !------------------------------------------------------------------------!
                  recp_found = .false.
                  recloopp: do recp=donp-1,1,-1

                     rec_lu = csite%dist_type(recp)

                     old_or_same_lu =   don_lu == rec_lu                         .or.      &
                                      ( csite%age(donp) >= min_oldgrowth(don_lu) .and.     &
                                        csite%age(recp) >= min_oldgrowth(rec_lu) )

                     recp_found = old_or_same_lu .and. fuse_table(recp) .and.              &
                                  (csite%dist_type(recp) == 1 .or. csite%age(recp) > 3.)

                     if (recp_found) then
                        recpatch => csite%patch(recp)
                        exit recloopp
                     end if
                  end do recloopp
                 
                  if (.not. recp_found) then
                     if (print_fuse_details) then
                        open  (unit=72,file=trim(fuse_fout),status='old',action='write'    &
                                                           ,position='append')
                        write (unit=72,fmt='(a)') '     - No receptor patch found. '
                        close (unit=72,status='keep')
                     end if
                     cycle donloopp
                  end if
                  !------------------------------------------------------------------------!

                  if (print_fuse_details) then
                     open  (unit=72,file=trim(fuse_fout),status='old',action='write'       &
                                                        ,position='append')
                     write (unit=72,fmt='(2(a,1x,i6,1x))') '     - DONP =',donp            &
                                                                 ,'RECP =',recp
                     close (unit=72,status='keep')
                  end if

                  !------------------------------------------------------------------------!
                  !     This should never happen because we have already fused all empty   !
                  ! patches, but, just in case... If both patches are empty they cannot be !
                  ! fused in this loop.                                                    !
                  !------------------------------------------------------------------------!
                  if ( dont_force_fuse        .and.                                        &
                       donpatch%ncohorts == 0 .and. recpatch%ncohorts == 0) then
                     cycle donloopp
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Find the LAI that corresponds to 80% of the maximum LAI, to avoid  !
                  ! relaxing too much for forests.                                         !
                  !------------------------------------------------------------------------!
                  dark_lai80 = 0.40 * ( sum(csite%cumlai_profile(:,1,recp))                &
                                      + sum(csite%cumlai_profile(:,1,donp)) )

                  !------------------------------------------------------------------------!
                  !     Compare the size profile for each PFT.  Here we compare the        !
                  ! maximum LAI for each PFT and height bin.  We switched the classes from !
                  ! DBH to height because different PFTs may have different heights for a  !
                  ! given DBH, so we want to make sure the light profile is right.         !
                  !------------------------------------------------------------------------!
                  hgtloop: do ihgt=1,ff_nhgt

                     cumlai_recp = sum(csite%cumlai_profile(:,ihgt,recp))
                     cumlai_donp = sum(csite%cumlai_profile(:,ihgt,donp))

                     !---------------------------------------------------------------------!
                     !    Check whether these bins contain some LAI.  We don't need to     !
                     ! check the cohorts if the understorey is too dark, so once both      !
                     ! patches becomes very dark (very high LAI), we stop checking the     !
                     ! profiles.                                                           !
                     !---------------------------------------------------------------------!
                     dark_donp = cumlai_donp > dark_toler
                     dark_recp = cumlai_recp > dark_toler
                     !---------------------------------------------------------------------!

                     !---------------------------------------------------------------------!
                     if (dark_recp .and. dark_donp) then

                        if (print_fuse_details) then
                           open  (unit=72,file=trim(fuse_fout),status='old',action='write' &
                                                              ,position='append')
                           write (unit=72,fmt='(1(a,1x,i6,1x),4(a,1x,es9.2,1x)'//          &
                                              ',2(a,1x,l1,1x))')                           &
                              '       * IHGT=',ihgt                                        &
                             ,'CUMLAI_RECP =',cumlai_recp,'CUMLAI_DONP =',cumlai_donp      &
                             ,'DARK_TOLER =',dark_toler,'DARK_LAI80 =',dark_lai80          &
                             ,'DARK_RECP =',dark_recp,'DARK_DONP =',dark_donp
                           close (unit=72,status='keep')
                        end if
                        cycle hgtloop
                     end if
                     !---------------------------------------------------------------------!

                     
                     !---------------------------------------------------------------------!
                     !    Check whether these bins contain some LAI.  Bins that have       !
                     ! tiny cumulative LAI may differ by a lot in the relative scale,      !
                     ! but the actual value is so small that we don't really care          !
                     ! whether they are relatively different.                              !
                     !---------------------------------------------------------------------!
                     sunny_donp = cumlai_donp <= sunny_toler
                     sunny_recp = cumlai_recp <= sunny_toler
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    If both patches have little or no biomass in this bin, don't     !
                     ! even bother checking the difference.                                !
                     !---------------------------------------------------------------------!
                     if (sunny_donp .and. sunny_recp) then

                        if (print_fuse_details) then
                           open  (unit=72,file=trim(fuse_fout),status='old',action='write' &
                                                              ,position='append')
                           write (unit=72,fmt='(1(a,1x,i6,1x),3(a,1x,es9.2,1x)'//          &
                                              ',2(a,1x,l1,1x))')                           &
                              '       * IHGT=',ihgt                                        &
                             ,'CUMLAI_RECP =',cumlai_recp,'CUMLAI_DONP =',cumlai_donp      &
                             ,'SUNNY_TOLER =',sunny_toler,'SUNNY_RECP =',sunny_recp        &
                             ,'SUNNY_RECP =',sunny_donp
                           close (unit=72,status='keep')
                        end if

                        cycle hgtloop
                     end if
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    Find the normalised difference in the density of this PFT and    !
                     ! size.  If one of the patches is missing any member of the           !
                     ! profile the norm will be set to 2.0, which is the highest value     !
                     ! that the norm can be.                                               !
                     !---------------------------------------------------------------------!
                     llevel_donp = exp(- 0.5 * cumlai_donp)
                     llevel_recp = exp(- 0.5 * cumlai_recp)
                     
                     diff = abs(llevel_donp - llevel_recp )
                     refv =    (llevel_donp + llevel_recp ) * 0.5
                     norm = diff / refv
                     fuse_flag = norm <= light_toler
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     if (print_fuse_details) then
                        open  (unit=72,file=trim(fuse_fout),status='old',action='write'    &
                                                           ,position='append')
                        write (unit=72,fmt='(1(a,1x,i6,1x),7(a,1x,es9.2,1x)'//             &
                                           ',1(a,1x,l1,1x))')                              &
                           '       * IHGT=',ihgt                                           &
                          ,'CLAI_RECP =',cumlai_recp,'CLAI_DONP =',cumlai_donp             &
                          ,'LL_RECP =',llevel_recp,'LL_DONP =',llevel_donp                 &
                          ,'DIFF =',diff,'REFV =',refv,'NORM =',norm                       &
                          ,'FUSE_FLAG =',fuse_flag
                        close (unit=72,status='keep')
                     end if
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !     If fuse_flag is false, the patches aren't similar, move to      !
                     ! the next donor patch.                                               !
                     !---------------------------------------------------------------------!
                     if (dont_force_fuse .and. (.not. fuse_flag)) cycle donloopp
                     !---------------------------------------------------------------------!
                  end do hgtloop
                  !------------------------------------------------------------------------!

                 

                  !------------------------------------------------------------------------!
                  !      Reaching this point means that the patches are sufficiently       !
                  ! similar so they will be fused.   We take the average of the patch      !
                  ! properties of donpatch and recpatch, and leave the averaged values at  !
                  ! recpatch.                                                              !
                  !------------------------------------------------------------------------!
                  call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)%prss          &
                                     ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)               &
                                     ,cpoly%green_leaf_factor(:,isi),fuse_initial          &
                                     ,elim_nplant,elim_lai)
                  !------------------------------------------------------------------------!

                  

                  !----- Record the fusion if requested by the user. ----------------------!
                  if (print_fuse_details) then
                     open (unit=72,file=trim(fuse_fout),status='old',action='write'        &
                                                        ,position='append')
                     write(unit=72,fmt='(2(a,i6),a)') '     * Patches ',donp,' and ',recp  &
                                                     ,' were fused.'
                     close(unit=72,status='keep')
                  end if
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !    Some cohorts may have been eliminated during the fusion process,    !
                  ! because they were way too small.  Add the eliminated plant density and !
                  ! LAI because we want to make sure that the fusion routine conserves the !
                  ! total plant density and LAI that remained in the polygon.              !
                  !------------------------------------------------------------------------!
                  elim_nplant_tot = elim_nplant_tot + elim_nplant
                  elim_lai_tot    = elim_lai_tot    + elim_lai
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Recalculate the pft size profile for the updated receptor patch.   !
                  !------------------------------------------------------------------------!
                  call patch_pft_size_profile(csite,recp)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     From now on donpatch should not be used: set fuse_table flag as    !
                  ! .false. so we won't check it again.                                    !
                  !------------------------------------------------------------------------!
                  fuse_table(donp) = .false.
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Update the number of valid patches.                                !
                  !------------------------------------------------------------------------!
                  npatches_new = npatches_new - 1
                  !------------------------------------------------------------------------!
               end do donloopp         ! do donp = csite%npatches,2,-1
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !      Check how many patches are valid.  If the total number of patches is !
               ! less than the target, or if we have reached the maximum tolerance and the !
               ! patch fusion still can't find similar patches, we quit the fusion loop.   !
               !---------------------------------------------------------------------------!
               if (npatches_new <= abs(maxpatch))  exit mainfuseloop
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Find the LAI that corresponds to 80% of the maximum LAI, to           !
               ! avoid relaxing too much for forests.                                      !
               !---------------------------------------------------------------------------!
               dark_lai80 = 0.
               now_area   = 0.
               do ipa=1,csite%npatches
                  if (fuse_table(ipa)) then
                     now_area   = now_area + csite%area(ipa)
                     dark_lai80 = dark_lai80                                               &
                                + 0.80 * sum(csite%cumlai_profile(:,1,ipa))                &
                                * csite%area(ipa)
                  end if
               end do
               if (now_area > 0.) dark_lai80 = dark_lai80 / now_area
               !---------------------------------------------------------------------------!

               !----- Increment tolerance -------------------------------------------------!
               sunny_toler =     sunny_toler * sunny_cumlai_mult
               dark_toler  = max( dark_toler * dark_cumlai_mult , dark_lai80 )
               light_toler =     light_toler * light_toler_mult
               !---------------------------------------------------------------------------!
            end do mainfuseloop
            !------------------------------------------------------------------------------!

            !----- Set the number of patches in the site to "npatches_new" ----------------!
            tempsite%npatches = npatches_new
            !------------------------------------------------------------------------------!



            !----- If there was any patch fusion, need to shrink csite --------------------!
            if (npatches_new < csite%npatches) then
               !---------------------------------------------------------------------------!
               !    Copy the selected data into the temporary space, args 1 and 3 must be  !
               ! dimension of arg 4. Argument 2 must be the dimension of the sum of the    !
               ! 3rd argument.                                                             !
               !---------------------------------------------------------------------------!
               call copy_sitetype_mask(csite,tempsite,fuse_table,size(fuse_table)          &
                                      ,npatches_new)
               call deallocate_sitetype(csite)

               !----- Reallocate the current site. ----------------------------------------!
               call allocate_sitetype(csite,npatches_new)

               !----- Copy the selected temporary data into the orignal site vectors. -----!
               call copy_sitetype(tempsite,csite,1,npatches_new,1,npatches_new)

               !---------------------------------------------------------------------------!
               !     The new and fused csite is now complete, clean up the temporary       !
               ! data. Deallocate it afterwards.                                           !
               !---------------------------------------------------------------------------!
               call deallocate_sitetype(tempsite)
            end if
            !------------------------------------------------------------------------------!



            !----- Deallocation should happen outside the "if" statement ------------------!
            deallocate(tempsite)
            deallocate(fuse_table)
            !------------------------------------------------------------------------------!


            !----- Make sure that patches are sorted from oldest to youngest. -------------!
            call sort_patches(csite)

            if (print_fuse_details) then
               open (unit=72,file=trim(fuse_fout),status='old',action='write'              &
                                                  ,position='append')
               write(unit=72,fmt='(a)')             '   + Patches were sorted. '
               write(unit=72,fmt='(2(a,1x,i6,1x))')                                        &
                                       '   + Number of patches.  Original =',npatches_orig &
                                                                ,'Current =',csite%npatches
               write(unit=72,fmt='(a)')       ' '
               close(unit=72,status='keep')
            end if


            !----- This is for mass conservation check ------------------------------------!
            new_nplant_tot = 0.
            new_lai_tot    = 0.
            new_area       = 0.
            do ipa = 1,csite%npatches
               new_area = new_area + csite%area(ipa)
               cpatch => csite%patch(ipa)
               do ico = 1, cpatch%ncohorts
                  new_nplant_tot = new_nplant_tot + cpatch%nplant(ico)*csite%area(ipa)
                  new_lai_tot    = new_lai_tot    + cpatch%lai(ico)*csite%area(ipa)
               end do
            end do
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Sanity check.  Except for the cohorts that were eliminated because they  !
            ! have become too small after fusion, the total plant count and LAI should be  !
            ! preserved.  In case something went wrong, we stop the run, it is likely to   !
            ! be a bug.                                                                    !
            !------------------------------------------------------------------------------!
            if (new_area       < 0.99 * old_area       .or.                                &
                new_area       > 1.01 * old_area            ) then
               write (unit=*,fmt='(a,1x,es12.5)') 'OLD_AREA:       ',old_area
               write (unit=*,fmt='(a,1x,es12.5)') 'NEW_AREA:       ',new_area
               write (unit=*,fmt='(a,1x,es12.5)') 'NEW_LAI_TOT:    ',new_lai_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'OLD_LAI_TOT:    ',old_lai_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'ELIM_LAI_TOT:   ',elim_lai_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'NEW_NPLANT_TOT: ',new_nplant_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'OLD_NPLANT_TOT: ',old_nplant_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'ELIM_NPLANT_TOT:',elim_nplant_tot
               call fatal_error('Conservation failed while fusing patches'                 &
                              &,'fuse_patches','fuse_fiss_utils.f90')
            end if
            !------------------------------------------------------------------------------!
            
         end do siteloop
      end do polyloop

      !------------------------------------------------------------------------------------!
      !     Print a banner to inform the user how many patches and cohorts exist.          !
      !------------------------------------------------------------------------------------!
      tot_npolygons = cgrid%npolygons
      tot_ncohorts  = 0
      tot_npatches  = 0
      tot_nsites    = 0
      do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         tot_nsites = tot_nsites + cpoly%nsites 
         do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            tot_npatches = tot_npatches + csite%npatches
            do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)
               tot_ncohorts = tot_ncohorts + cpatch%ncohorts
            end do
         end do
      end do
      write (unit=*,fmt='(6(a,1x,i8,1x))')                                                 &
        'Total count in node',mynum,'for grid',ifm,': POLYGONS=',tot_npolygons             &
       ,'SITES=',tot_nsites,'PATCHES=',tot_npatches,'COHORTS=',tot_ncohorts
      !------------------------------------------------------------------------------------!

      return
   end subroutine fuse_patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will merge two patches into 1.                                      !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_patches(csite,donp,recp,mzg,mzs,prss,lsl,ntext_soil,green_leaf_factor &
                            ,fuse_initial,elim_nplant,elim_lai)
      use ed_state_vars      , only : sitetype              & ! Structure 
                                    , patchtype             ! ! Structure
      use soil_coms          , only : soil                  & ! intent(in), lookup table
                                    , tiny_sfcwater_mass    & ! intent(in)
                                    , matric_potential      ! ! intent(in)
      use ed_max_dims        , only : n_pft                 & ! intent(in)
                                    , n_dbh                 ! ! intent(in)
      use mem_polygons       , only : maxcohort             ! ! intent(in)
      use therm_lib          , only : uextcm2tl             & ! subroutine
                                    , uint2tl               & ! subroutine
                                    , idealdenssh           & ! function
                                    , press2exner           & ! function
                                    , extheta2temp          ! ! function
      use ed_misc_coms       , only : writing_long          & ! intent(in)
                                    , writing_eorq          & ! intent(in)
                                    , writing_dcyc          & ! intent(in)
                                    , ndcycle               ! ! intent(in)
      use budget_utils       , only : update_budget         ! ! intent(in)
      use consts_coms        , only : wdns                  ! ! intent(in)
      use fusion_fission_coms, only : corr_patch            ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target      :: csite             ! Current site
      integer                , intent(in)  :: donp              ! Donating patch
      integer                , intent(in)  :: recp              ! Receptor patch
      integer                , intent(in)  :: lsl               ! Lowest soil level
      integer                , intent(in)  :: mzg               ! # of soil layers
      integer                , intent(in)  :: mzs               ! # of sfc. water layers
      integer, dimension(mzg), intent(in)  :: ntext_soil        ! Soil type
      real, dimension(n_pft) , intent(in)  :: green_leaf_factor ! Green leaf factor...
      real                   , intent(in)  :: prss              ! Sfc. air density
      logical                , intent(in)  :: fuse_initial      ! Initialisation?
      real                   , intent(out) :: elim_nplant       ! Eliminated nplant 
      real                   , intent(out) :: elim_lai          ! Eliminated lai
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer     :: cpatch            ! Current patch
      type(patchtype)        , pointer     :: temppatch         ! Temporary patch
      integer                              :: iii               ! Counters
      integer                              :: nsoil             ! Alias for soil texture
      integer                              :: t                 ! Counter for time of day
      integer                              :: ndc               ! # of cohorts - donp patch
      integer                              :: nrc               ! # of cohorts - recp patch
      real                                 :: can_exner         ! Exner function - CAS
      real                                 :: newarea           ! new patch area
      real                                 :: newareai          ! 1./(new patch area)
      real                                 :: area_scale        ! Cohort rescaling factor.
      !------------------------------------------------------------------------------------!
     
      !------------------------------------------------------------------------------------!
      !     This function fuses the two patches specified in the argument. It fuses the    !
      ! first patch in the argument (the "donor" = donp ) into the second patch in the     !
      ! argument (the "recipient" = recp ), and frees the memory associated with the donor !
      ! patch.                                                                             !
      !------------------------------------------------------------------------------------!
    
      !----- The new area is simply the sum of each patch area. ---------------------------!
      newarea  = csite%area(donp) + csite%area(recp)
      newareai = 1.0/newarea
      !------------------------------------------------------------------------------------!

      !----- Assign eliminated LAI and nplant to zero (everything stays) ------------------!
      elim_nplant = 0.
      elim_lai    = 0.
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     In case of old-growth stands, the disturbance type flag may be different.  In  !
      ! this case, we keep the type for the largest patch.                                 !
      !------------------------------------------------------------------------------------!
      if (csite%area(donp) > csite%area(recp)) then
         csite%dist_type(recp) = csite%dist_type(donp)
      end if
      !------------------------------------------------------------------------------------!



      !----- We now take the weighted average, scale by the individual patch area. --------!
      csite%age(recp)                = newareai *                                          &
                                     ( csite%age(donp)                * csite%area(donp)   &
                                     + csite%age(recp)                * csite%area(recp) )

      csite%fast_soil_C(recp)        = newareai *                                          &
                                     ( csite%fast_soil_C(donp)        * csite%area(donp)   &
                                     + csite%fast_soil_C(recp)        * csite%area(recp) )

      csite%slow_soil_C(recp)        = newareai *                                          &
                                     ( csite%slow_soil_C(donp)        * csite%area(donp)   &
                                     + csite%slow_soil_C(recp)        * csite%area(recp) )

      csite%structural_soil_C(recp)  = newareai *                                          &
                                     ( csite%structural_soil_C(donp)  * csite%area(donp)   &
                                     + csite%structural_soil_C(recp)  * csite%area(recp) )
                                     

      csite%structural_soil_L(recp)  = newareai *                                          &
                                     ( csite%structural_soil_L(donp)  * csite%area(donp)   &
                                     + csite%structural_soil_L(recp)  * csite%area(recp) )
                                     

      csite%mineralized_soil_N(recp) = newareai *                                          &
                                     ( csite%mineralized_soil_N(donp) * csite%area(donp)   &
                                     + csite%mineralized_soil_N(recp) * csite%area(recp) )

      csite%fast_soil_N(recp)        = newareai *                                          &
                                     ( csite%fast_soil_N(donp)        * csite%area(donp)   &
                                     + csite%fast_soil_N(recp)        * csite%area(recp) )

      csite%sum_dgd(recp)            = newareai *                                          &
                                     ( csite%sum_dgd(donp)            * csite%area(donp)   &
                                     + csite%sum_dgd(recp)            * csite%area(recp) )

      csite%sum_chd(recp)            = newareai *                                          &
                                     ( csite%sum_chd(donp)            * csite%area(donp)   &
                                     + csite%sum_chd(recp)            * csite%area(recp) )

      csite%can_co2(recp)            = newareai *                                          &
                                     ( csite%can_co2(donp)            * csite%area(donp)   &
                                     + csite%can_co2(recp)            * csite%area(recp) )

      csite%can_theta(recp)          = newareai *                                          &
                                     ( csite%can_theta(donp)          * csite%area(donp)   &
                                     + csite%can_theta(recp)          * csite%area(recp) )

      csite%can_temp_pv(recp)        = newareai *                                          &
                                     ( csite%can_temp_pv(donp)        * csite%area(donp)   &
                                     + csite%can_temp_pv(recp)        * csite%area(recp) )

      csite%can_theiv(recp)          = newareai *                                          &
                                     ( csite%can_theiv(donp)          * csite%area(donp)   &
                                     + csite%can_theiv(recp)          * csite%area(recp) )

      csite%can_vpdef(recp)          = newareai *                                          &
                                     ( csite%can_vpdef(donp)          * csite%area(donp)   &
                                     + csite%can_vpdef(recp)          * csite%area(recp) )

      csite%can_prss(recp)           = newareai *                                          &
                                     ( csite%can_prss(donp)           * csite%area(donp)   &
                                     + csite%can_prss(recp)           * csite%area(recp) )

      csite%can_shv(recp)            = newareai *                                          &
                                     ( csite%can_shv(donp)            * csite%area(donp)   &
                                     + csite%can_shv(recp)            * csite%area(recp) )

      csite%can_depth(recp)          = newareai *                                          &
                                     ( csite%can_depth(donp)          * csite%area(donp)   &
                                     + csite%can_depth(recp)          * csite%area(recp) )

      csite%ggbare(recp)             = newareai *                                          &
                                     ( csite%ggbare(donp)             * csite%area(donp)   &
                                     + csite%ggbare(recp)             * csite%area(recp) )

      csite%ggnet(recp)              = newareai *                                          &
                                     ( csite%ggnet(donp)              * csite%area(donp)   &
                                     + csite%ggnet(recp)              * csite%area(recp) )

      csite%ggsoil(recp)             = newareai *                                          &
                                     ( csite%ggsoil(donp)             * csite%area(donp)   &
                                     + csite%ggsoil(recp)             * csite%area(recp) )

      
      !------------------------------------------------------------------------------------!
      !    There is no guarantee that there will be a minimum amount of mass in the tempo- !
      ! rary layer, nor is there any reason for both patches to have the same number of    !
      ! layers. In order to be safe, the fusion must happen in 5 stages.                   !
      !------------------------------------------------------------------------------------!
      !----- 1. Find the "extensive" sfcwater_energy (convert from J/kg to J/m2); ---------!
      do iii=1,csite%nlev_sfcwater(recp)
         csite%sfcwater_energy(iii,recp) = csite%sfcwater_energy(iii,recp)                 &
                                         * csite%sfcwater_mass  (iii,recp)
      end do
      do iii=1,csite%nlev_sfcwater(donp)
         csite%sfcwater_energy(iii,donp) = csite%sfcwater_energy(iii,donp)                 &
                                         * csite%sfcwater_mass  (iii,donp)
      end do
      !------------------------------------------------------------------------------------!
      ! 2. Squeeze all layers into one.  If needed, the layer will be split again next     !
      !    time the Runge-Kutta integrator is called.  After adding the value to the first !
      !    layer, discard the value.                                                       !
      !------------------------------------------------------------------------------------!
      do iii=2,csite%nlev_sfcwater(recp)
         csite%sfcwater_energy(1,recp) = csite%sfcwater_energy(1,recp)                     &
                                       + csite%sfcwater_energy(iii,recp)
         csite%sfcwater_depth(1,recp)  = csite%sfcwater_depth(1,recp)                      &
                                       + csite%sfcwater_depth(iii,recp)
         csite%sfcwater_mass(1,recp)   = csite%sfcwater_mass(1,recp)                       &
                                       + csite%sfcwater_mass(iii,recp)
         csite%sfcwater_energy(iii,recp) = 0.
         csite%sfcwater_depth(iii,recp)  = 0.
         csite%sfcwater_mass(iii,recp)   = 0.
      end do
      do iii=2,csite%nlev_sfcwater(donp)
         csite%sfcwater_energy(1,donp) = csite%sfcwater_energy(1,donp)                     &
                                       + csite%sfcwater_energy(iii,donp)
         csite%sfcwater_depth(1,donp)  = csite%sfcwater_depth(1,donp)                      &
                                       + csite%sfcwater_depth(iii,donp)
         csite%sfcwater_mass(1,donp)   = csite%sfcwater_mass(1,donp)                       &
                                       + csite%sfcwater_mass(iii,donp)
         csite%sfcwater_energy(iii,donp) = 0.
         csite%sfcwater_depth(iii,donp)  = 0.
         csite%sfcwater_mass(iii,donp)   = 0.
      end do
      !----- 3. Merge the patches; --------------------------------------------------------!
      if (csite%nlev_sfcwater(recp) > 0 .or. csite%nlev_sfcwater(donp) > 0) then
         csite%sfcwater_mass(1,recp)   = newareai *                                        &
                                         (csite%sfcwater_mass(1,recp)  * csite%area(recp)  &
                                         +csite%sfcwater_mass(1,donp)  * csite%area(donp)  )
         csite%sfcwater_depth(1,recp)  = newareai *                                        &
                                         (csite%sfcwater_depth(1,recp) * csite%area(recp)  &
                                         +csite%sfcwater_depth(1,donp) * csite%area(donp)  )
         csite%sfcwater_energy(1,recp) = newareai *                                        &
                                         (csite%sfcwater_energy(1,recp) * csite%area(recp) &
                                         +csite%sfcwater_energy(1,donp) * csite%area(donp) )
      else
         csite%sfcwater_mass(1,recp)   = 0.
         csite%sfcwater_depth(1,recp)  = 0.
         csite%sfcwater_energy(1,recp) = 0.
      end if
      !------------------------------------------------------------------------------------!
      ! 4. Converting energy back to J/kg;                                                 !
      ! 5. Finding temperature and liquid water fraction;                                  !
      !    (Both are done in new_patch_sfc_props).                                         !
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
     

      !----- Merge soil energy and water. -------------------------------------------------!
      do iii=1,mzg
         csite%soil_energy(iii,recp)     = newareai *                                      &
                                         ( csite%soil_energy(iii,donp) * csite%area(donp)  &
                                         + csite%soil_energy(iii,recp) * csite%area(recp))

         csite%soil_water(iii,recp)      = newareai *                                      &
                                         ( csite%soil_water(iii,recp)  * csite%area(recp)  &
                                         + csite%soil_water(iii,donp)  * csite%area(donp))
      end do

      !------------------------------------------------------------------------------------!
      !     These variables shouldn't matter because they are reset every day/every month, !
      ! but just in case...                                                                !
      !------------------------------------------------------------------------------------!
      csite%avg_daily_temp       (recp) = newareai                                         &
                                        * ( csite%avg_daily_temp       (donp)              &
                                          * csite%area                 (donp)              &
                                          + csite%avg_daily_temp       (recp)              &
                                          * csite%area                 (recp) )

      csite%avg_monthly_gndwater (recp) = newareai                                         &
                                        * ( csite%avg_monthly_gndwater (donp)              &
                                          * csite%area                 (donp)              &
                                          + csite%avg_monthly_gndwater (recp)              &
                                          * csite%area                 (recp) )

      csite%avg_monthly_waterdef (recp) = newareai                                         &
                                        * ( csite%avg_monthly_waterdef (donp)              &
                                          * csite%area                 (donp)              &
                                          + csite%avg_monthly_waterdef (recp)              &
                                          * csite%area                 (recp) )
      !------------------------------------------------------------------------------------!


      

      !------------------------------------------------------------------------------------!
      !    This subroutine takes care of filling:                                          !
      !                                                                                    !
      ! + csite%ground_shv(recp)                                                           !
      ! + csite%ground_ssh(recp)                                                           !
      ! + csite%ground_temp(recp)                                                          !
      ! + csite%ground_fliq(recp)                                                          !
      ! + csite%soil_tempk(k,recp)                                                         !
      ! + csite%soil_fracliq(k,recp)                                                       !
      ! + csite%nlev_sfcwater(recp)                                                        !
      ! + csite%sfcwater_energy(k,recp) (Just converting back to J/kg)                     !
      ! + csite%csite%sfcwater_tempk(k,recp)                                               !
      ! + csite%sfcwater_fracliq(k,recp)                                                   !
      !------------------------------------------------------------------------------------!
      call new_patch_sfc_props(csite,recp,mzg,mzs,ntext_soil)
      !------------------------------------------------------------------------------------!

      csite%today_A_decomp(recp)         = newareai *                                      &
                                         ( csite%today_A_decomp(donp) * csite%area(donp)   &
                                         + csite%today_A_decomp(recp) * csite%area(recp) )

      csite%today_Af_decomp(recp)        = newareai *                                      &
                                         ( csite%today_Af_decomp(donp)* csite%area(donp)   &
                                         + csite%today_Af_decomp(recp)* csite%area(recp) )

      do iii = 1,n_pft
         csite%repro(iii,recp)           = newareai *                                      &
                                         ( csite%repro(iii,donp)      * csite%area(donp)   &
                                         + csite%repro(iii,recp)      * csite%area(recp) )
      end do
      !------------------------------------------------------------------------------------!


 


      !------------------------------------------------------------------------------------!
      !     Budget variables.                                                              !
      !------------------------------------------------------------------------------------!
      csite%co2budget_residual    (recp) = ( csite%co2budget_residual    (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%co2budget_residual    (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%co2budget_loss2atm    (recp) = ( csite%co2budget_loss2atm    (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%co2budget_loss2atm    (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%co2budget_denseffect  (recp) = ( csite%co2budget_denseffect  (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%co2budget_denseffect  (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%co2budget_gpp         (recp) = ( csite%co2budget_gpp         (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%co2budget_gpp         (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%co2budget_plresp      (recp) = ( csite%co2budget_plresp      (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%co2budget_plresp      (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%co2budget_rh          (recp) = ( csite%co2budget_rh          (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%co2budget_rh          (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%ebudget_residual      (recp) = ( csite%ebudget_residual      (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%ebudget_residual      (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%ebudget_netrad        (recp) = ( csite%ebudget_netrad        (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%ebudget_netrad        (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%ebudget_loss2atm      (recp) = ( csite%ebudget_loss2atm      (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%ebudget_loss2atm      (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%ebudget_denseffect    (recp) = ( csite%ebudget_denseffect    (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%ebudget_denseffect    (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%ebudget_prsseffect    (recp) = ( csite%ebudget_prsseffect    (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%ebudget_prsseffect    (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%ebudget_loss2runoff   (recp) = ( csite%ebudget_loss2runoff   (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%ebudget_loss2runoff   (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%ebudget_loss2drainage (recp) = ( csite%ebudget_loss2drainage (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%ebudget_loss2drainage (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%ebudget_precipgain    (recp) = ( csite%ebudget_precipgain    (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%ebudget_precipgain    (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%wbudget_residual      (recp) = ( csite%wbudget_residual      (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%wbudget_residual      (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%wbudget_loss2atm      (recp) = ( csite%wbudget_loss2atm      (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%wbudget_loss2atm      (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%wbudget_denseffect    (recp) = ( csite%wbudget_denseffect    (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%wbudget_denseffect    (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%wbudget_loss2runoff   (recp) = ( csite%wbudget_loss2runoff   (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%wbudget_loss2runoff   (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%wbudget_loss2drainage (recp) = ( csite%wbudget_loss2drainage (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%wbudget_loss2drainage (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      csite%wbudget_precipgain    (recp) = ( csite%wbudget_precipgain    (recp)            &
                                           * csite%area                  (recp)            &
                                           + csite%wbudget_precipgain    (donp)            &
                                           * csite%area                  (donp) )          &
                                         * newareai
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Sub-daily means.  This should be skipped during initialisation because all      !
      ! averages will be zero.                                                             !
      !------------------------------------------------------------------------------------!
      if (.not. fuse_initial) then
         csite%fmean_rh                   (recp) = ( csite%fmean_rh                (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_rh                (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_cwd_rh               (recp) = ( csite%fmean_cwd_rh            (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_cwd_rh            (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_nep                  (recp) = ( csite%fmean_nep               (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_nep               (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_rk4step              (recp) = ( csite%fmean_rk4step           (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_rk4step           (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_available_water      (recp) = ( csite%fmean_available_water   (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_available_water   (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_can_theiv            (recp) = ( csite%fmean_can_theiv         (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_can_theiv         (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_can_theta            (recp) = ( csite%fmean_can_theta         (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_can_theta         (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_can_vpdef            (recp) = ( csite%fmean_can_vpdef         (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_can_vpdef         (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_can_shv              (recp) = ( csite%fmean_can_shv           (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_can_shv           (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_can_co2              (recp) = ( csite%fmean_can_co2           (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_can_co2           (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_can_prss             (recp) = ( csite%fmean_can_prss          (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_can_prss          (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_gnd_temp             (recp) = ( csite%fmean_gnd_temp          (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_gnd_temp          (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_gnd_shv              (recp) = ( csite%fmean_gnd_shv           (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_gnd_shv           (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_can_ggnd             (recp) = ( csite%fmean_can_ggnd          (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_can_ggnd          (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_rshort_gnd           (recp) = ( csite%fmean_rshort_gnd        (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_rshort_gnd        (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_par_gnd              (recp) = ( csite%fmean_par_gnd           (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_par_gnd           (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_rlong_gnd            (recp) = ( csite%fmean_rlong_gnd         (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_rlong_gnd         (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_rlongup              (recp) = ( csite%fmean_rlongup           (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_rlongup           (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_parup                (recp) = ( csite%fmean_parup             (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_parup             (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_nirup                (recp) = ( csite%fmean_nirup             (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_nirup             (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_rshortup             (recp) = ( csite%fmean_rshortup          (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_rshortup          (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_rnet                 (recp) = ( csite%fmean_rnet              (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_rnet              (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_albedo               (recp) = ( csite%fmean_albedo            (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_albedo            (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_albedo_par           (recp) = ( csite%fmean_albedo_par        (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_albedo_par        (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_albedo_nir           (recp) = ( csite%fmean_albedo_nir        (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_albedo_nir        (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_rlong_albedo         (recp) = ( csite%fmean_rlong_albedo      (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_rlong_albedo      (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_ustar                (recp) = ( csite%fmean_ustar             (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_ustar             (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_tstar                (recp) = ( csite%fmean_tstar             (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_tstar             (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_qstar                (recp) = ( csite%fmean_qstar             (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_qstar             (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_cstar                (recp) = ( csite%fmean_cstar             (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_cstar             (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_carbon_ac            (recp) = ( csite%fmean_carbon_ac         (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_carbon_ac         (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_carbon_st            (recp) = ( csite%fmean_carbon_st         (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_carbon_st         (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_vapor_gc             (recp) = ( csite%fmean_vapor_gc          (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_vapor_gc          (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_vapor_ac             (recp) = ( csite%fmean_vapor_ac          (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_vapor_ac          (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_throughfall          (recp) = ( csite%fmean_throughfall       (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_throughfall       (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_runoff               (recp) = ( csite%fmean_runoff            (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_runoff            (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_drainage             (recp) = ( csite%fmean_drainage          (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_drainage          (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_sensible_gc          (recp) = ( csite%fmean_sensible_gc       (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_sensible_gc       (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_sensible_ac          (recp) = ( csite%fmean_sensible_ac       (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_sensible_ac       (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_qthroughfall         (recp) = ( csite%fmean_qthroughfall      (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_qthroughfall      (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_qrunoff              (recp) = ( csite%fmean_qrunoff           (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_qrunoff           (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_qdrainage            (recp) = ( csite%fmean_qdrainage         (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_qdrainage         (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_soil_energy        (:,recp) = ( csite%fmean_soil_energy     (:,recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_soil_energy     (:,donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_soil_water         (:,recp) = ( csite%fmean_soil_water      (:,recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_soil_water      (:,donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_smoist_gg          (:,recp) = ( csite%fmean_smoist_gg       (:,recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_smoist_gg       (:,donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_transloss          (:,recp) = ( csite%fmean_transloss       (:,recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_transloss       (:,donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_sensible_gg        (:,recp) = ( csite%fmean_sensible_gg     (:,recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_sensible_gg     (:,donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         !---------------------------------------------------------------------------------!

        
         !---------------------------------------------------------------------------------!
         !      Now we find the derived properties for the canopy air space.               !
         !---------------------------------------------------------------------------------!

         can_exner                   = press2exner (csite%fmean_can_prss(recp))
         csite%fmean_can_temp (recp) = extheta2temp(can_exner,csite%fmean_can_theta (recp))
         csite%fmean_can_rhos (recp) = idealdenssh ( csite%fmean_can_prss  (recp)          &
                                                   , csite%fmean_can_temp  (recp)          &
                                                   , csite%fmean_can_shv   (recp)          )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Find the soil mean temperature, liquid water fraction, and matric          !
         ! potential.                                                                      !
         !---------------------------------------------------------------------------------!
         do iii=lsl,mzg
            nsoil = ntext_soil(iii)
            call uextcm2tl( csite%fmean_soil_energy(iii,recp)                              &
                          , csite%fmean_soil_water (iii,recp) * wdns                       &
                          , soil(nsoil)%slcpd                                              &
                          , csite%fmean_soil_temp  (iii,recp)                              &
                          , csite%fmean_soil_fliq  (iii,recp))

            csite%fmean_soil_mstpot(iii,recp)  =                                           &
                                  matric_potential(nsoil,csite%fmean_soil_water (iii,recp))
         end do
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Find the temporary surface water properties.  They may not be available at  !
         ! all times, so we must check.                                                    !
         !---------------------------------------------------------------------------------!
         !----- Temporarily make energy extensive [J/m2]. ---------------------------------!
         csite%fmean_sfcw_depth           (recp) = ( csite%fmean_sfcw_depth        (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_sfcw_depth        (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_sfcw_energy          (recp) = ( csite%fmean_sfcw_energy       (recp)  &
                                                   * csite%fmean_sfcw_mass         (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_sfcw_energy       (donp)  &
                                                   * csite%fmean_sfcw_mass         (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         csite%fmean_sfcw_mass            (recp) = ( csite%fmean_sfcw_mass         (recp)  &
                                                   * csite%area                    (recp)  &
                                                   + csite%fmean_sfcw_mass         (donp)  &
                                                   * csite%area                    (donp)) &
                                                 *   newareai
         !----- Check whether there is enough surface water. ------------------------------!
         if (csite%fmean_sfcw_mass(recp) > tiny_sfcwater_mass) then
            csite%fmean_sfcw_energy       (recp) =   csite%fmean_sfcw_energy       (recp)  &
                                                 /   csite%fmean_sfcw_mass         (recp)
            call uint2tl( csite%fmean_sfcw_energy(recp)                                    &
                        , csite%fmean_sfcw_temp  (recp)                                    &
                        , csite%fmean_sfcw_fliq  (recp) )
         else
            csite%fmean_sfcw_mass  (recp)  = 0.
            csite%fmean_sfcw_depth (recp)  = 0.
            csite%fmean_sfcw_energy(recp)  = 0.
            csite%fmean_sfcw_temp  (recp)  = csite%fmean_soil_temp(mzg,recp)
            csite%fmean_sfcw_fliq  (recp)  = csite%fmean_soil_fliq(mzg,recp)
         end if
         !---------------------------------------------------------------------------------!

      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------! 
      !    Daily means.                                                                    !
      !------------------------------------------------------------------------------------! 
      if (writing_long .and.  (.not. fuse_initial) ) then
         csite%dmean_A_decomp           (recp) = ( csite%dmean_A_decomp           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_A_decomp           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_Af_decomp          (recp) = ( csite%dmean_Af_decomp          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_Af_decomp          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_co2_residual       (recp) = ( csite%dmean_co2_residual       (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_co2_residual       (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_energy_residual    (recp) = ( csite%dmean_energy_residual    (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_energy_residual    (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_water_residual     (recp) = ( csite%dmean_water_residual     (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_water_residual     (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_rh                 (recp) = ( csite%dmean_rh                 (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_rh                 (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_cwd_rh             (recp) = ( csite%dmean_cwd_rh             (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_cwd_rh             (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_nep                (recp) = ( csite%dmean_nep                (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_nep                (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_rk4step            (recp) = ( csite%dmean_rk4step            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_rk4step            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_available_water    (recp) = ( csite%dmean_available_water    (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_available_water    (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_can_theiv          (recp) = ( csite%dmean_can_theiv          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_can_theiv          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_can_theta          (recp) = ( csite%dmean_can_theta          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_can_theta          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_can_vpdef          (recp) = ( csite%dmean_can_vpdef          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_can_vpdef          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_can_shv            (recp) = ( csite%dmean_can_shv            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_can_shv            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_can_co2            (recp) = ( csite%dmean_can_co2            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_can_co2            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_can_prss           (recp) = ( csite%dmean_can_prss           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_can_prss           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_gnd_temp           (recp) = ( csite%dmean_gnd_temp           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_gnd_temp           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_gnd_shv            (recp) = ( csite%dmean_gnd_shv            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_gnd_shv            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_can_ggnd           (recp) = ( csite%dmean_can_ggnd           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_can_ggnd           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_rshort_gnd         (recp) = ( csite%dmean_rshort_gnd         (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_rshort_gnd         (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_par_gnd            (recp) = ( csite%dmean_par_gnd            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_par_gnd            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_rlong_gnd          (recp) = ( csite%dmean_rlong_gnd          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_rlong_gnd          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_rlongup            (recp) = ( csite%dmean_rlongup            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_rlongup            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_parup              (recp) = ( csite%dmean_parup              (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_parup              (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_nirup              (recp) = ( csite%dmean_nirup              (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_nirup              (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_rshortup           (recp) = ( csite%dmean_rshortup           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_rshortup           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_rnet               (recp) = ( csite%dmean_rnet               (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_rnet               (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_albedo             (recp) = ( csite%dmean_albedo             (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_albedo             (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_albedo_par         (recp) = ( csite%dmean_albedo_par         (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_albedo_par         (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_albedo_nir         (recp) = ( csite%dmean_albedo_nir         (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_albedo_nir         (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_rlong_albedo       (recp) = ( csite%dmean_rlong_albedo       (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_rlong_albedo       (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_ustar              (recp) = ( csite%dmean_ustar              (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_ustar              (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_tstar              (recp) = ( csite%dmean_tstar              (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_tstar              (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_qstar              (recp) = ( csite%dmean_qstar              (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_qstar              (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_cstar              (recp) = ( csite%dmean_cstar              (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_cstar              (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_carbon_ac          (recp) = ( csite%dmean_carbon_ac          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_carbon_ac          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_carbon_st          (recp) = ( csite%dmean_carbon_st          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_carbon_st          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_vapor_gc           (recp) = ( csite%dmean_vapor_gc           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_vapor_gc           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_vapor_ac           (recp) = ( csite%dmean_vapor_ac           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_vapor_ac           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_throughfall        (recp) = ( csite%dmean_throughfall        (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_throughfall        (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_runoff             (recp) = ( csite%dmean_runoff             (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_runoff             (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_drainage           (recp) = ( csite%dmean_drainage           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_drainage           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_sensible_gc        (recp) = ( csite%dmean_sensible_gc        (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_sensible_gc        (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_sensible_ac        (recp) = ( csite%dmean_sensible_ac        (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_sensible_ac        (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_qthroughfall       (recp) = ( csite%dmean_qthroughfall       (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_qthroughfall       (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_qrunoff            (recp) = ( csite%dmean_qrunoff            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_qrunoff            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_qdrainage          (recp) = ( csite%dmean_qdrainage          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_qdrainage          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_smoist_gg        (:,recp) = ( csite%dmean_smoist_gg        (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_smoist_gg        (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_transloss        (:,recp) = ( csite%dmean_transloss        (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_transloss        (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_sensible_gg      (:,recp) = ( csite%dmean_sensible_gg      (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_sensible_gg      (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Now we find the derived properties for the canopy air space.               !
         !---------------------------------------------------------------------------------!
         can_exner                   = press2exner (csite%dmean_can_prss(recp))
         csite%dmean_can_temp (recp) = extheta2temp(can_exner,csite%dmean_can_theta(recp))
         csite%dmean_can_rhos (recp) = idealdenssh ( csite%dmean_can_prss  (recp)          &
                                                   , csite%dmean_can_temp  (recp)          &
                                                   , csite%dmean_can_shv   (recp)          )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Find the soil mean temperature, liquid water fraction, and matric          !
         ! potential.                                                                      !
         !---------------------------------------------------------------------------------!
         do iii=lsl,mzg
            nsoil = ntext_soil(iii)
            call uextcm2tl( csite%dmean_soil_energy(iii,recp)                              &
                          , csite%dmean_soil_water (iii,recp) * wdns                       &
                          , soil(nsoil)%slcpd                                              &
                          , csite%dmean_soil_temp  (iii,recp)                              &
                          , csite%dmean_soil_fliq  (iii,recp))

            csite%dmean_soil_mstpot(iii,recp)  =                                           &
                                   matric_potential(nsoil,csite%dmean_soil_water(iii,recp))
         end do
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the temporary surface water properties.  They may not be available at  !
         ! all times, so we must check.                                                    !
         !---------------------------------------------------------------------------------!
         !----- Temporarily make energy extensive [J/m2]. ---------------------------------!
         csite%dmean_sfcw_depth         (recp) = ( csite%dmean_sfcw_depth         (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_sfcw_depth         (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_sfcw_energy        (recp) = ( csite%dmean_sfcw_energy        (recp)   &
                                                 * csite%dmean_sfcw_mass          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_sfcw_energy        (donp)   &
                                                 * csite%dmean_sfcw_mass          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%dmean_sfcw_mass          (recp) = ( csite%dmean_sfcw_mass          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%dmean_sfcw_mass          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         !----- Check whether there is enough surface water. ------------------------------!
         if (csite%dmean_sfcw_mass(recp) > tiny_sfcwater_mass) then
            csite%dmean_sfcw_energy     (recp) =   csite%dmean_sfcw_energy        (recp)   &
                                               /   csite%dmean_sfcw_mass          (recp)
            call uint2tl( csite%dmean_sfcw_energy(recp)                                    &
                        , csite%dmean_sfcw_temp  (recp)                                    &
                        , csite%dmean_sfcw_fliq  (recp) )
         else
            csite%dmean_sfcw_mass  (recp)  = 0.
            csite%dmean_sfcw_depth (recp)  = 0.
            csite%dmean_sfcw_energy(recp)  = 0.
            csite%dmean_sfcw_temp  (recp)  = csite%dmean_soil_temp(mzg,recp)
            csite%dmean_sfcw_fliq  (recp)  = csite%dmean_soil_fliq(mzg,recp)
         end if
         !------------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------! 
      !    Monthly means.                                                                  !
      !------------------------------------------------------------------------------------! 
      if (writing_eorq .and. (.not. fuse_initial)) then

         !---------------------------------------------------------------------------------!
         !    First we find the mean sum of squares, because they depend on the means too, !
         ! and the receptor values are lost after fusion.  All variables are intensive at  !
         ! the patch level.                                                                !
         !---------------------------------------------------------------------------------!

         csite%mmsqu_rh         (recp) = fuse_msqu( csite%mmean_rh         (recp)          &
                                                  , csite%mmsqu_rh         (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_rh         (donp)          &
                                                  , csite%mmsqu_rh         (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_cwd_rh     (recp) = fuse_msqu( csite%mmean_cwd_rh     (recp)          &
                                                  , csite%mmsqu_cwd_rh     (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_cwd_rh     (donp)          &
                                                  , csite%mmsqu_cwd_rh     (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_nep        (recp) = fuse_msqu( csite%mmean_nep        (recp)          &
                                                  , csite%mmsqu_nep        (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_nep        (donp)          &
                                                  , csite%mmsqu_nep        (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_rlongup    (recp) = fuse_msqu( csite%mmean_rlongup    (recp)          &
                                                  , csite%mmsqu_rlongup    (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_rlongup    (donp)          &
                                                  , csite%mmsqu_rlongup    (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_parup      (recp) = fuse_msqu( csite%mmean_parup      (recp)          &
                                                  , csite%mmsqu_parup      (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_parup      (donp)          &
                                                  , csite%mmsqu_parup      (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_nirup      (recp) = fuse_msqu( csite%mmean_nirup      (recp)          &
                                                  , csite%mmsqu_nirup      (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_nirup      (donp)          &
                                                  , csite%mmsqu_nirup      (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_rshortup   (recp) = fuse_msqu( csite%mmean_rshortup   (recp)          &
                                                  , csite%mmsqu_rshortup   (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_rshortup   (donp)          &
                                                  , csite%mmsqu_rshortup   (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_rnet       (recp) = fuse_msqu( csite%mmean_rnet       (recp)          &
                                                  , csite%mmsqu_rnet       (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_rnet       (donp)          &
                                                  , csite%mmsqu_rnet       (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_albedo     (recp) = fuse_msqu( csite%mmean_albedo     (recp)          &
                                                  , csite%mmsqu_albedo     (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_albedo     (donp)          &
                                                  , csite%mmsqu_albedo     (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_ustar      (recp) = fuse_msqu( csite%mmean_ustar      (recp)          &
                                                  , csite%mmsqu_ustar      (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_ustar      (donp)          &
                                                  , csite%mmsqu_ustar      (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_carbon_ac  (recp) = fuse_msqu( csite%mmean_carbon_ac  (recp)          &
                                                  , csite%mmsqu_carbon_ac  (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_carbon_ac  (donp)          &
                                                  , csite%mmsqu_carbon_ac  (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_carbon_st  (recp) = fuse_msqu( csite%mmean_carbon_st  (recp)          &
                                                  , csite%mmsqu_carbon_st  (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_carbon_st  (donp)          &
                                                  , csite%mmsqu_carbon_st  (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_vapor_gc   (recp) = fuse_msqu( csite%mmean_vapor_gc   (recp)          &
                                                  , csite%mmsqu_vapor_gc   (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_vapor_gc   (donp)          &
                                                  , csite%mmsqu_vapor_gc   (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_vapor_ac   (recp) = fuse_msqu( csite%mmean_vapor_ac   (recp)          &
                                                  , csite%mmsqu_vapor_ac   (recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_vapor_ac   (donp)          &
                                                  , csite%mmsqu_vapor_ac   (donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_sensible_gc(recp) = fuse_msqu( csite%mmean_sensible_gc(recp)          &
                                                  , csite%mmsqu_sensible_gc(recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_sensible_gc(donp)          &
                                                  , csite%mmsqu_sensible_gc(donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         csite%mmsqu_sensible_ac(recp) = fuse_msqu( csite%mmean_sensible_ac(recp)          &
                                                  , csite%mmsqu_sensible_ac(recp)          &
                                                  , csite%area             (recp)          &
                                                  , csite%mmean_sensible_ac(donp)          &
                                                  , csite%mmsqu_sensible_ac(donp)          &
                                                  , csite%area             (donp)          &
                                                  , corr_patch, .false.)
         !---------------------------------------------------------------------------------! 


         csite%mmean_rh                 (recp) = ( csite%mmean_rh                 (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_rh                 (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_cwd_rh             (recp) = ( csite%mmean_cwd_rh             (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_cwd_rh             (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_nep                (recp) = ( csite%mmean_nep                (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_nep                (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_rk4step            (recp) = ( csite%mmean_rk4step            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_rk4step            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_available_water    (recp) = ( csite%mmean_available_water    (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_available_water    (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_can_theiv          (recp) = ( csite%mmean_can_theiv          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_can_theiv          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_can_theta          (recp) = ( csite%mmean_can_theta          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_can_theta          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_can_vpdef          (recp) = ( csite%mmean_can_vpdef          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_can_vpdef          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_can_shv            (recp) = ( csite%mmean_can_shv            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_can_shv            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_can_co2            (recp) = ( csite%mmean_can_co2            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_can_co2            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_can_prss           (recp) = ( csite%mmean_can_prss           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_can_prss           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_gnd_temp           (recp) = ( csite%mmean_gnd_temp           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_gnd_temp           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_gnd_shv            (recp) = ( csite%mmean_gnd_shv            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_gnd_shv            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_can_ggnd           (recp) = ( csite%mmean_can_ggnd           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_can_ggnd           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_rshort_gnd         (recp) = ( csite%mmean_rshort_gnd         (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_rshort_gnd         (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_par_gnd            (recp) = ( csite%mmean_par_gnd            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_par_gnd            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_rlong_gnd          (recp) = ( csite%mmean_rlong_gnd          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_rlong_gnd          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_rlongup            (recp) = ( csite%mmean_rlongup            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_rlongup            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_parup              (recp) = ( csite%mmean_parup              (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_parup              (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_nirup              (recp) = ( csite%mmean_nirup              (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_nirup              (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_rshortup           (recp) = ( csite%mmean_rshortup           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_rshortup           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_rnet               (recp) = ( csite%mmean_rnet               (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_rnet               (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_albedo             (recp) = ( csite%mmean_albedo             (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_albedo             (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_albedo_par         (recp) = ( csite%mmean_albedo_par         (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_albedo_par         (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_albedo_nir         (recp) = ( csite%mmean_albedo_nir         (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_albedo_nir         (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_rlong_albedo       (recp) = ( csite%mmean_rlong_albedo       (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_rlong_albedo       (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_ustar              (recp) = ( csite%mmean_ustar              (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_ustar              (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_tstar              (recp) = ( csite%mmean_tstar              (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_tstar              (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_qstar              (recp) = ( csite%mmean_qstar              (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_qstar              (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_cstar              (recp) = ( csite%mmean_cstar              (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_cstar              (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_carbon_ac          (recp) = ( csite%mmean_carbon_ac          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_carbon_ac          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_carbon_st          (recp) = ( csite%mmean_carbon_st          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_carbon_st          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_vapor_gc           (recp) = ( csite%mmean_vapor_gc           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_vapor_gc           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_vapor_ac           (recp) = ( csite%mmean_vapor_ac           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_vapor_ac           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_throughfall        (recp) = ( csite%mmean_throughfall        (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_throughfall        (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_runoff             (recp) = ( csite%mmean_runoff             (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_runoff             (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_drainage           (recp) = ( csite%mmean_drainage           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_drainage           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_sensible_gc        (recp) = ( csite%mmean_sensible_gc        (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_sensible_gc        (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_sensible_ac        (recp) = ( csite%mmean_sensible_ac        (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_sensible_ac        (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_qthroughfall       (recp) = ( csite%mmean_qthroughfall       (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_qthroughfall       (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_qrunoff            (recp) = ( csite%mmean_qrunoff            (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_qrunoff            (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_qdrainage          (recp) = ( csite%mmean_qdrainage          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_qdrainage          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_fast_soil_c        (recp) = ( csite%mmean_fast_soil_c        (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_fast_soil_c        (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_slow_soil_c        (recp) = ( csite%mmean_slow_soil_c        (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_slow_soil_c        (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_struct_soil_c      (recp) = ( csite%mmean_struct_soil_c      (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_struct_soil_c      (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_struct_soil_l      (recp) = ( csite%mmean_struct_soil_l      (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_struct_soil_l      (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_fast_soil_n        (recp) = ( csite%mmean_fast_soil_n        (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_fast_soil_n        (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_mineral_soil_n     (recp) = ( csite%mmean_mineral_soil_n     (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_mineral_soil_n     (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_A_decomp           (recp) = ( csite%mmean_A_decomp           (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_A_decomp           (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_Af_decomp          (recp) = ( csite%mmean_Af_decomp          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_Af_decomp          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_co2_residual       (recp) = ( csite%mmean_co2_residual       (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_co2_residual       (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_energy_residual    (recp) = ( csite%mmean_energy_residual    (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_energy_residual    (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_water_residual     (recp) = ( csite%mmean_water_residual     (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_water_residual     (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_smoist_gg        (:,recp) = ( csite%mmean_smoist_gg        (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_smoist_gg        (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_transloss        (:,recp) = ( csite%mmean_transloss        (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_transloss        (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_sensible_gg      (:,recp) = ( csite%mmean_sensible_gg      (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_sensible_gg      (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Now we find the derived properties for the canopy air space.               !
         !---------------------------------------------------------------------------------!
         can_exner                   = press2exner (csite%mmean_can_prss(recp))
         csite%mmean_can_temp (recp) = extheta2temp(can_exner,csite%mmean_can_theta(recp))
         csite%mmean_can_rhos (recp) = idealdenssh ( csite%mmean_can_prss  (recp)          &
                                                   , csite%mmean_can_temp  (recp)          &
                                                   , csite%mmean_can_shv   (recp)          )
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Find the soil mean temperature, liquid water fraction, and matric          !
         ! potential.                                                                      !
         !---------------------------------------------------------------------------------!
         do iii=lsl,mzg
            nsoil = ntext_soil(iii)
            call uextcm2tl( csite%mmean_soil_energy(iii,recp)                              &
                          , csite%mmean_soil_water (iii,recp) * wdns                       &
                          , soil(nsoil)%slcpd                                              &
                          , csite%mmean_soil_temp  (iii,recp)                              &
                          , csite%mmean_soil_fliq  (iii,recp))

            csite%mmean_soil_mstpot(iii,recp) =                                            &
                                   matric_potential(nsoil,csite%mmean_soil_water(iii,recp))
         end do
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Find the temporary surface water properties.  They may not be available at  !
         ! all times, so we must check.                                                    !
         !---------------------------------------------------------------------------------!
         !----- Temporarily make energy extensive [J/m2]. ---------------------------------!
         csite%mmean_sfcw_depth         (recp) = ( csite%mmean_sfcw_depth         (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_sfcw_depth         (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_sfcw_energy        (recp) = ( csite%mmean_sfcw_energy        (recp)   &
                                                 * csite%mmean_sfcw_mass          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_sfcw_energy        (donp)   &
                                                 * csite%mmean_sfcw_mass          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%mmean_sfcw_mass          (recp) = ( csite%mmean_sfcw_mass          (recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%mmean_sfcw_mass          (donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         !----- Check whether there is enough surface water. ------------------------------!
         if (csite%mmean_sfcw_mass(recp) > tiny_sfcwater_mass) then
            csite%mmean_sfcw_energy     (recp) =   csite%mmean_sfcw_energy        (recp)   &
                                               /   csite%mmean_sfcw_mass          (recp)
            call uint2tl( csite%mmean_sfcw_energy(recp)                                    &
                        , csite%mmean_sfcw_temp  (recp)                                    &
                        , csite%mmean_sfcw_fliq  (recp) )
         else
            csite%mmean_sfcw_mass  (recp)  = 0.
            csite%mmean_sfcw_depth (recp)  = 0.
            csite%mmean_sfcw_energy(recp)  = 0.
            csite%mmean_sfcw_temp  (recp)  = csite%mmean_soil_temp(mzg,recp)
            csite%mmean_sfcw_fliq  (recp)  = csite%mmean_soil_fliq(mzg,recp)
         end if
         !------------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------! 
      !    Mean diel.                                                                      !
      !------------------------------------------------------------------------------------! 
      if (writing_dcyc .and. (.not. fuse_initial)) then
         !---------------------------------------------------------------------------------!
         !      First we solve the mean sum of squares as they depend on the mean and the  !
         ! original receptor data is lost after fusion takes place.                        !
         !---------------------------------------------------------------------------------!
         do t=1,ndcycle
            csite%qmsqu_rh         (t,recp) = fuse_msqu( csite%qmean_rh         (t,recp)   &
                                                       , csite%qmsqu_rh         (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_rh         (t,donp)   &
                                                       , csite%qmsqu_rh         (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_cwd_rh     (t,recp) = fuse_msqu( csite%qmean_cwd_rh     (t,recp)   &
                                                       , csite%qmsqu_cwd_rh     (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_cwd_rh     (t,donp)   &
                                                       , csite%qmsqu_cwd_rh     (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_nep        (t,recp) = fuse_msqu( csite%qmean_nep        (t,recp)   &
                                                       , csite%qmsqu_nep        (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_nep        (t,donp)   &
                                                       , csite%qmsqu_nep        (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_rlongup    (t,recp) = fuse_msqu( csite%qmean_rlongup    (t,recp)   &
                                                       , csite%qmsqu_rlongup    (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_rlongup    (t,donp)   &
                                                       , csite%qmsqu_rlongup    (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_parup      (t,recp) = fuse_msqu( csite%qmean_parup      (t,recp)   &
                                                       , csite%qmsqu_parup      (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_parup      (t,donp)   &
                                                       , csite%qmsqu_parup      (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_nirup      (t,recp) = fuse_msqu( csite%qmean_nirup      (t,recp)   &
                                                       , csite%qmsqu_nirup      (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_nirup      (t,donp)   &
                                                       , csite%qmsqu_nirup      (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_rshortup   (t,recp) = fuse_msqu( csite%qmean_rshortup   (t,recp)   &
                                                       , csite%qmsqu_rshortup   (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_rshortup   (t,donp)   &
                                                       , csite%qmsqu_rshortup   (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_rnet       (t,recp) = fuse_msqu( csite%qmean_rnet       (t,recp)   &
                                                       , csite%qmsqu_rnet       (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_rnet       (t,donp)   &
                                                       , csite%qmsqu_rnet       (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_albedo     (t,recp) = fuse_msqu( csite%qmean_albedo     (t,recp)   &
                                                       , csite%qmsqu_albedo     (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_albedo     (t,donp)   &
                                                       , csite%qmsqu_albedo     (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_ustar      (t,recp) = fuse_msqu( csite%qmean_ustar      (t,recp)   &
                                                       , csite%qmsqu_ustar      (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_ustar      (t,donp)   &
                                                       , csite%qmsqu_ustar      (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_carbon_ac  (t,recp) = fuse_msqu( csite%qmean_carbon_ac  (t,recp)   &
                                                       , csite%qmsqu_carbon_ac  (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_carbon_ac  (t,donp)   &
                                                       , csite%qmsqu_carbon_ac  (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_carbon_st  (t,recp) = fuse_msqu( csite%qmean_carbon_st  (t,recp)   &
                                                       , csite%qmsqu_carbon_st  (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_carbon_st  (t,donp)   &
                                                       , csite%qmsqu_carbon_st  (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_vapor_gc   (t,recp) = fuse_msqu( csite%qmean_vapor_gc   (t,recp)   &
                                                       , csite%qmsqu_vapor_gc   (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_vapor_gc   (t,donp)   &
                                                       , csite%qmsqu_vapor_gc   (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_vapor_ac   (t,recp) = fuse_msqu( csite%qmean_vapor_ac   (t,recp)   &
                                                       , csite%qmsqu_vapor_ac   (t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_vapor_ac   (t,donp)   &
                                                       , csite%qmsqu_vapor_ac   (t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_sensible_gc(t,recp) = fuse_msqu( csite%qmean_sensible_gc(t,recp)   &
                                                       , csite%qmsqu_sensible_gc(t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_sensible_gc(t,donp)   &
                                                       , csite%qmsqu_sensible_gc(t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
            csite%qmsqu_sensible_ac(t,recp) = fuse_msqu( csite%qmean_sensible_ac(t,recp)   &
                                                       , csite%qmsqu_sensible_ac(t,recp)   &
                                                       , csite%area               (recp)   &
                                                       , csite%qmean_sensible_ac(t,donp)   &
                                                       , csite%qmsqu_sensible_ac(t,donp)   &
                                                       , csite%area               (donp)   &
                                                       , corr_patch, .false.)
         end do
         !---------------------------------------------------------------------------------!


         csite%qmean_rh               (:,recp) = ( csite%qmean_rh               (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_rh               (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_cwd_rh           (:,recp) = ( csite%qmean_cwd_rh           (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_cwd_rh           (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_nep              (:,recp) = ( csite%qmean_nep              (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_nep              (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_rk4step          (:,recp) = ( csite%qmean_rk4step          (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_rk4step          (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_available_water  (:,recp) = ( csite%qmean_available_water  (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_available_water  (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_can_theiv        (:,recp) = ( csite%qmean_can_theiv        (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_can_theiv        (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_can_theta        (:,recp) = ( csite%qmean_can_theta        (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_can_theta        (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_can_vpdef        (:,recp) = ( csite%qmean_can_vpdef        (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_can_vpdef        (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_can_shv          (:,recp) = ( csite%qmean_can_shv          (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_can_shv          (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_can_co2          (:,recp) = ( csite%qmean_can_co2          (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_can_co2          (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_can_prss         (:,recp) = ( csite%qmean_can_prss         (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_can_prss         (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_gnd_temp         (:,recp) = ( csite%qmean_gnd_temp         (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_gnd_temp         (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_gnd_shv          (:,recp) = ( csite%qmean_gnd_shv          (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_gnd_shv          (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_can_ggnd         (:,recp) = ( csite%qmean_can_ggnd         (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_can_ggnd         (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_rshort_gnd       (:,recp) = ( csite%qmean_rshort_gnd       (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_rshort_gnd       (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_par_gnd          (:,recp) = ( csite%qmean_par_gnd          (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_par_gnd          (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_rlong_gnd        (:,recp) = ( csite%qmean_rlong_gnd        (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_rlong_gnd        (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_rlongup          (:,recp) = ( csite%qmean_rlongup          (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_rlongup          (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_parup            (:,recp) = ( csite%qmean_parup            (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_parup            (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_nirup            (:,recp) = ( csite%qmean_nirup            (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_nirup            (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_rshortup         (:,recp) = ( csite%qmean_rshortup         (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_rshortup         (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_rnet             (:,recp) = ( csite%qmean_rnet             (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_rnet             (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_albedo           (:,recp) = ( csite%qmean_albedo           (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_albedo           (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_albedo_par       (:,recp) = ( csite%qmean_albedo_par       (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_albedo_par       (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_albedo_nir       (:,recp) = ( csite%qmean_albedo_nir       (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_albedo_nir       (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_rlong_albedo     (:,recp) = ( csite%qmean_rlong_albedo     (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_rlong_albedo     (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_ustar            (:,recp) = ( csite%qmean_ustar            (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_ustar            (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_tstar            (:,recp) = ( csite%qmean_tstar            (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_tstar            (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_qstar            (:,recp) = ( csite%qmean_qstar            (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_qstar            (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_cstar            (:,recp) = ( csite%qmean_cstar            (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_cstar            (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_carbon_ac        (:,recp) = ( csite%qmean_carbon_ac        (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_carbon_ac        (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_carbon_st        (:,recp) = ( csite%qmean_carbon_st        (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_carbon_st        (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_vapor_gc         (:,recp) = ( csite%qmean_vapor_gc         (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_vapor_gc         (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_vapor_ac         (:,recp) = ( csite%qmean_vapor_ac         (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_vapor_ac         (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_throughfall      (:,recp) = ( csite%qmean_throughfall      (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_throughfall      (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_runoff           (:,recp) = ( csite%qmean_runoff           (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_runoff           (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_drainage         (:,recp) = ( csite%qmean_drainage         (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_drainage         (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_sensible_gc      (:,recp) = ( csite%qmean_sensible_gc      (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_sensible_gc      (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_sensible_ac      (:,recp) = ( csite%qmean_sensible_ac      (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_sensible_ac      (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_qthroughfall     (:,recp) = ( csite%qmean_qthroughfall     (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_qthroughfall     (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_qrunoff          (:,recp) = ( csite%qmean_qrunoff          (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_qrunoff          (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_qdrainage        (:,recp) = ( csite%qmean_qdrainage        (:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_qdrainage        (:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_smoist_gg      (:,:,recp) = ( csite%qmean_smoist_gg      (:,:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_smoist_gg      (:,:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_transloss      (:,:,recp) = ( csite%qmean_transloss      (:,:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_transloss      (:,:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai
         csite%qmean_sensible_gg    (:,:,recp) = ( csite%qmean_sensible_gg    (:,:,recp)   &
                                                 * csite%area                     (recp)   &
                                                 + csite%qmean_sensible_gg    (:,:,donp)   &
                                                 * csite%area                     (donp) ) &
                                               *   newareai


         do t=1,ndcycle
            !------------------------------------------------------------------------------!
            !      Now we find the derived properties for the canopy air space.            !
            !------------------------------------------------------------------------------!
            can_exner                     = press2exner (csite%qmean_can_prss(t,recp))
            csite%qmean_can_temp (t,recp) = extheta2temp( can_exner                        &
                                                        , csite%qmean_can_theta(t,recp))
            csite%qmean_can_rhos (t,recp) = idealdenssh ( csite%qmean_can_prss (t,recp)    &
                                                        , csite%qmean_can_temp (t,recp)    &
                                                        , csite%qmean_can_shv  (t,recp) )
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Find the soil mean temperature, liquid water fraction, and matric       !
            ! potential.                                                                   !
            !------------------------------------------------------------------------------!
            do iii=lsl,mzg
               nsoil = ntext_soil(iii)
               call uextcm2tl( csite%qmean_soil_energy(iii,t,recp)                         &
                             , csite%qmean_soil_water (iii,t,recp) * wdns                  &
                             , soil(nsoil)%slcpd                                           &
                             , csite%qmean_soil_temp  (iii,t,recp)                         &
                             , csite%qmean_soil_fliq  (iii,t,recp) )

               csite%qmean_soil_mstpot(iii,t,recp) =                                       &
                                 matric_potential(nsoil,csite%qmean_soil_water(iii,t,recp))
            end do
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find the temporary surface water properties.  They may not be available  !
            ! at all times, so we must check.                                              !
            !------------------------------------------------------------------------------!
            !----- Temporarily make energy extensive [J/m2]. ------------------------------!
            csite%qmean_sfcw_depth       (t,recp) = ( csite%qmean_sfcw_depth    (t,recp)   &
                                                    * csite%area                  (recp)   &
                                                    + csite%qmean_sfcw_depth    (t,donp)   &
                                                    * csite%area                  (donp) ) &
                                                  *   newareai
            csite%qmean_sfcw_energy      (t,recp) = ( csite%qmean_sfcw_energy   (t,recp)   &
                                                    * csite%qmean_sfcw_mass     (t,recp)   &
                                                    * csite%area                  (recp)   &
                                                    + csite%qmean_sfcw_energy   (t,donp)   &
                                                    * csite%qmean_sfcw_mass     (t,donp)   &
                                                    * csite%area                  (donp) ) &
                                                  *   newareai
            csite%qmean_sfcw_mass        (t,recp) = ( csite%qmean_sfcw_mass     (t,recp)   &
                                                    * csite%area                  (recp)   &
                                                    + csite%qmean_sfcw_mass     (t,donp)   &
                                                    * csite%area                  (donp) ) &
                                                  *   newareai
            !----- Check whether there is enough surface water. ---------------------------!
            if (csite%qmean_sfcw_mass(t,recp) > tiny_sfcwater_mass) then
               csite%qmean_sfcw_energy   (t,recp) =   csite%qmean_sfcw_energy   (t,recp)   &
                                                  /   csite%qmean_sfcw_mass     (t,recp)
               call uint2tl( csite%qmean_sfcw_energy(t,recp)                               &
                           , csite%qmean_sfcw_temp  (t,recp)                               &
                           , csite%qmean_sfcw_fliq  (t,recp) )
            else
               csite%qmean_sfcw_mass  (t,recp)  = 0.
               csite%qmean_sfcw_depth (t,recp)  = 0.
               csite%qmean_sfcw_energy(t,recp)  = 0.
               csite%qmean_sfcw_temp  (t,recp)  = csite%qmean_soil_temp(mzg,t,recp)
               csite%qmean_sfcw_fliq  (t,recp)  = csite%qmean_soil_fliq(mzg,t,recp)
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    We now update the canopy thermodynamic propeties:                               !
      ! + csite%can_temp(recp)                                                             !
      ! + csite%can_rhos(recp)                                                             !
      !------------------------------------------------------------------------------------!
      call update_patch_thermo_props(csite,recp,recp,mzg,mzs,ntext_soil)

      !------------------------------------------------------------------------------------!
      !     Now we need to adjust the densities of cohorts. Because the patch area         !
      ! increased we want to retain the same total amount of mass and energy.              !
      !------------------------------------------------------------------------------------!
      !----- 1. Adjust densities of cohorts in recipient patch ----------------------------!
      cpatch => csite%patch(recp)
      nrc = cpatch%ncohorts
      area_scale = csite%area(recp) * newareai
      call update_cohort_extensive_props(cpatch,1,nrc,area_scale)
      !----- 2. Adjust densities of cohorts in donor patch --------------------------------!
      cpatch => csite%patch(donp)
      ndc = cpatch%ncohorts
      area_scale = csite%area(donp) * newareai
      call update_cohort_extensive_props(cpatch,1,ndc,area_scale)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Fill a new patch with the donor and recipient cohort vectors.                   !
      !------------------------------------------------------------------------------------!
      !----- Allocate the temporary patch with room for all cohorts. ----------------------!
      if (ndc + nrc > 0 ) then
         nullify(temppatch)
         allocate(temppatch)
         call allocate_patchtype(temppatch,ndc + nrc )
         !----- Copy the recipient and donor cohorts to the temporary patch. --------------!
         call copy_patchtype(csite%patch(recp),temppatch,1,nrc,1,nrc)
         call copy_patchtype(csite%patch(donp),temppatch,1,ndc,nrc+1,nrc+ndc)
         !----- Reallocate the current recipient patch with all cohorts -------------------!
         call deallocate_patchtype(csite%patch(recp))
         call allocate_patchtype(csite%patch(recp),ndc+nrc)
         !----- Copy the temporary patch back to the recipient patch. ---------------------!
         call copy_patchtype(temppatch,csite%patch(recp),1,nrc+ndc,1,nrc+ndc)
         !----- Get rid of the temporary patch --------------------------------------------!
         call deallocate_patchtype(temppatch)
         deallocate(temppatch)
         !----- Sort cohorts in the new patch ---------------------------------------------!
         cpatch => csite%patch(recp)
         call sort_cohorts(cpatch)
         !---------------------------------------------------------------------------------!
         !    We just combined two patches, so we may be able to fuse some cohorts and/or  !
         ! eliminate others.                                                               !
         !---------------------------------------------------------------------------------!
         if (cpatch%ncohorts > 0 .and. maxcohort >= 0) then
            call fuse_cohorts(csite,recp,green_leaf_factor,lsl,fuse_initial)
            call terminate_cohorts(csite,recp,elim_nplant,elim_lai)
            call split_cohorts(cpatch,green_leaf_factor,lsl)
         end if
         !---------------------------------------------------------------------------------!
      end if

      !------------------------------------------------------------------------------------!
      !    Now we update some variables that depend on cohort statistics, namely:          !
      ! + csite%veg_height(recp)                                                           !
      ! + csite%veg_displace(recp)                                                         !
      ! + csite%disp_height(recp)                                                          !
      ! + csite%veg_rough(recp)                                                            !
      ! + csite%total_sfcw_depth(recp)                                                     !
      ! + csite%snowfac(recp)                                                              !
      ! + csite%opencan_frac(recp)                                                         !
      !------------------------------------------------------------------------------------!
      call update_patch_derived_props(csite,recp)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Now we update the budget variables:                                             !
      ! + csite%wbudget_initialstorage(recp)                                               !
      ! + csite%ebudget_initialstorage(recp)                                               !
      ! + csite%co2budget_initialstorage(recp)                                             !
      !------------------------------------------------------------------------------------!
      call update_budget(csite,lsl,recp,recp)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    This subroutine will update the size profile within patch.                      !
      ! + csite%cumlai_profile(:,:,recp)                                                   !
      !------------------------------------------------------------------------------------!
      call patch_pft_size_profile(csite,recp)
      !------------------------------------------------------------------------------------!

      !----- Last, but not the least, we update the patch area ----------------------------!
      csite%area(recp) = newarea

      return

   end subroutine fuse_2_patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine patch_pft_size_profile(csite,ipa)
      use ed_state_vars       , only : sitetype   & ! structure
                                     , patchtype  ! ! structure
      use fusion_fission_coms , only : ff_nhgt    & ! intent(in)
                                     , hgt_class  ! ! intent(in)
      use allometry           , only : size2bl    ! ! intent(in)
      use ed_max_dims         , only : n_pft      ! ! intent(in)
      use pft_coms            , only : hgt_min    & ! intent(in)
                                     , is_grass   ! ! intent(in)
      use ed_misc_coms        , only : igrass     ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target     :: csite     ! Current site
      integer                , intent(in) :: ipa       ! Current patch index
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer    :: cpatch    ! Current patch
      integer                             :: ipft      ! PFT index
      integer                             :: ihgt      ! Height class index
      integer                             :: ico       ! Counters
      real                                :: lai_pot   ! Potential LAI
      !------------------------------------------------------------------------------------!


      !----- Reset all bins to zero. ------------------------------------------------------!
      do ipft=1,n_pft
         do ihgt=1,ff_nhgt
            csite%cumlai_profile(ipft,ihgt,ipa)=0.0
         end do
      end do
      !------------------------------------------------------------------------------------!



      !----- Update bins ------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      cohortloop: do ico = 1,cpatch%ncohorts

         !----- Find the PFT class. -------------------------------------------------------!
         ipft    = cpatch%pft(ico)
         ihgt    = min(ff_nhgt,max(1,count(hgt_class < cpatch%hite(ico))))
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Check whether this cohort is almost at the minimum height given its PFT.    !
         ! If it is, then we will skip it.                                                 !
         !---------------------------------------------------------------------------------!
         if (cpatch%hite(ico) < hgt_min(ipft) + 0.2) cycle cohortloop
         !---------------------------------------------------------------------------------!


         !----- Find the height class. ----------------------------------------------------!
         ihgt    = min(ff_nhgt,max(1,count(hgt_class < cpatch%hite(ico))))
         !---------------------------------------------------------------------------------!


         !----- Find the potential (on-allometry) leaf area index. ------------------------!
         if (is_grass(ipft) .and. igrass==1) then
             !--use actual bleaf for grass
             lai_pot = cpatch%nplant(ico) * cpatch%sla(ico) * cpatch%bleaf(ico)
         else
             !--use dbh for trees
             lai_pot = cpatch%nplant(ico) * cpatch%sla(ico)                                &
                     * size2bl(cpatch%dbh(ico),cpatch%hite(ico),ipft)
         end if
         !---------------------------------------------------------------------------------!


         !----- Add the potential LAI to the bin. -----------------------------------------!
         csite%cumlai_profile(ipft,ihgt,ipa) = lai_pot                                     &
                                             + csite%cumlai_profile(ipft,ihgt,ipa)
         !---------------------------------------------------------------------------------!
      end do cohortloop
      !------------------------------------------------------------------------------------!



      !----- Integrate the leaf area index from top to bottom. ----------------------------!
      do ihgt=ff_nhgt-1,1,-1
         do ipft=1,n_pft
            csite%cumlai_profile(ipft,ihgt,ipa) = csite%cumlai_profile(ipft,ihgt  ,ipa)    &
                                                + csite%cumlai_profile(ipft,ihgt+1,ipa)
         end do
      end do
      !------------------------------------------------------------------------------------!

      return
   end subroutine patch_pft_size_profile
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine combines the mean sum of squares of two quantities (x and y).     !
   !                                                                                       !
   ! xmean, ymean -- the mean values of x and y                                            !
   ! xmsqu, ymsqu -- the mean sum of squares of x and y                                    !
   ! xwght, ywght -- the weights for x and y                                               !
   ! r_xy         -- correlation between x and y                                           !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function fuse_msqu(xmean,xmsqu,xwght,ymean,ymsqu,ywght,r_xy,extensive)
      use rk4_coms       , only : tiny_offset              ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: xmean
      real(kind=4), intent(in) :: xmsqu
      real(kind=4), intent(in) :: xwght
      real(kind=4), intent(in) :: ymean
      real(kind=4), intent(in) :: ymsqu
      real(kind=4), intent(in) :: ywght
      real(kind=4), intent(in) :: r_xy
      logical     , intent(in) :: extensive
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: xwmean
      real(kind=8)             :: xwmsqu
      real(kind=8)             :: ywmean
      real(kind=8)             :: ywmsqu
      real(kind=8)             :: w2sumi
      real(kind=8)             :: xwpywp
      real(kind=8)             :: r_xy8
      real(kind=8)             :: fuse_msqu8
      !----- External function. -----------------------------------------------------------!
      real(kind=4)    , external      :: sngloff     ! Safe double -> single precision
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Correct the scales.  If the property is extensive, we aggregate them, other-  !
      ! wise we find the weighted sum of squares.                                          !
      !------------------------------------------------------------------------------------!
      if (extensive) then
         xwmean = dble(xmean)
         xwmsqu = dble(xmsqu)
         ywmean = dble(ymean)
         ywmsqu = dble(ymsqu)
         w2sumi = 1.d0
      else
         xwmean = dble(xmean) * dble(xwght)
         xwmsqu = dble(xmsqu) * dble(xwght) * dble(xwght)
         ywmean = dble(ymean) * dble(ywght)
         ywmsqu = dble(ymsqu) * dble(ywght) * dble(ywght)
         w2sumi = 1.d0 / ( (dble(xwght) + dble(ywght)) * (dble(xwght) + dble(ywght)) )
      end if
      r_xy8     = dble(r_xy)
      !------------------------------------------------------------------------------------!


      !----- Find the terms. --------------------------------------------------------------!
      xwpywp = r_xy8 * sqrt( ( xwmsqu + xwmean*xwmean) * (ywmsqu + ywmean*ywmean) )
      !------------------------------------------------------------------------------------!


      !----- Add the terms to the answer. -------------------------------------------------!
      fuse_msqu8 = ( xwmsqu + ywmsqu + 2.d0 * ( xwmean * ywmean + xwpywp ) ) * w2sumi
      fuse_msqu  = sngloff(fuse_msqu8,tiny_offset)
      !------------------------------------------------------------------------------------!

      return
   end function fuse_msqu
   !=======================================================================================!
   !=======================================================================================!
end module fuse_fiss_utils
!==========================================================================================!
!==========================================================================================!
