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
      integer                            :: ico, inew    ! Counters
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
      use disturb_coms , only : min_new_patch_area ! ! intent(in)
      use ed_misc_coms , only : iqoutput           & ! intent(in)
                              , imoutput           & ! intent(in)
                              , idoutput           ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)       , target      :: csite        ! Current site
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)       , pointer     :: tempsite     ! Scratch site
      type(patchtype)      , pointer     :: cpatch       ! Pointer to current site
      integer                            :: ipa,ico      ! Counters
      logical, dimension(:), allocatable :: remain_table ! Flag: this patch will remain.
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
      elim_area = 0.0
      do ipa = 1,csite%npatches
         if (csite%area(ipa) < min_new_patch_area) then
            elim_area = elim_area + csite%area(ipa)
            remain_table(ipa) = .false.
         end if
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

      !------------------------------------------------------------------------------------!
      !    Renormalize the total area.  We must also rescale all extensive properties from !
      ! cohorts, since they are per unit area and we are effectively changing the area.    !
      ! IMPORTANT: Only cohort-level variables that have units per area (m2) should be     !
      !            rescaled.  Variables whose units are per plant or per leaf area         !
      !            (m2_leaf) should _NOT_ be included here.                                !
      !------------------------------------------------------------------------------------!
      new_area=0.
      area_scale = 1./(1. - elim_area)
      do ipa = 1,csite%npatches
         csite%area(ipa) = csite%area(ipa) * area_scale
         new_area = new_area + csite%area(ipa)

         cpatch => csite%patch(ipa)
         do ico = 1, cpatch%ncohorts
            cpatch%nplant               (ico) = cpatch%nplant            (ico) * area_scale
            cpatch%lai                  (ico) = cpatch%lai               (ico) * area_scale
            cpatch%wpa                  (ico) = cpatch%wpa               (ico) * area_scale
            cpatch%wai                  (ico) = cpatch%wai               (ico) * area_scale
            cpatch%mean_gpp             (ico) = cpatch%mean_gpp          (ico) * area_scale
            cpatch%mean_leaf_resp       (ico) = cpatch%mean_leaf_resp    (ico) * area_scale
            cpatch%mean_root_resp       (ico) = cpatch%mean_root_resp    (ico) * area_scale
            cpatch%mean_growth_resp     (ico) = cpatch%mean_growth_resp  (ico) * area_scale
            cpatch%mean_storage_resp    (ico) = cpatch%mean_storage_resp (ico) * area_scale
            cpatch%mean_vleaf_resp      (ico) = cpatch%mean_vleaf_resp   (ico) * area_scale
            cpatch%gpp                  (ico) = cpatch%gpp               (ico) * area_scale
            cpatch%leaf_respiration     (ico) = cpatch%leaf_respiration  (ico) * area_scale
            cpatch%root_respiration     (ico) = cpatch%root_respiration  (ico) * area_scale
            cpatch%leaf_water           (ico) = cpatch%leaf_water        (ico) * area_scale
            cpatch%leaf_hcap            (ico) = cpatch%leaf_hcap         (ico) * area_scale
            cpatch%leaf_energy          (ico) = cpatch%leaf_energy       (ico) * area_scale
            cpatch%wood_water           (ico) = cpatch%wood_water        (ico) * area_scale
            cpatch%wood_hcap            (ico) = cpatch%wood_hcap         (ico) * area_scale
            cpatch%wood_energy          (ico) = cpatch%wood_energy       (ico) * area_scale
            cpatch%monthly_dndt         (ico) = cpatch%monthly_dndt      (ico) * area_scale
            cpatch%today_gpp            (ico) = cpatch%today_gpp         (ico) * area_scale
            cpatch%today_nppleaf        (ico) = cpatch%today_nppleaf     (ico) * area_scale
            cpatch%today_nppfroot       (ico) = cpatch%today_nppfroot    (ico) * area_scale
            cpatch%today_nppsapwood     (ico) = cpatch%today_nppsapwood  (ico) * area_scale
            cpatch%today_nppcroot       (ico) = cpatch%today_nppcroot    (ico) * area_scale
            cpatch%today_nppseeds       (ico) = cpatch%today_nppseeds    (ico) * area_scale
            cpatch%today_nppwood        (ico) = cpatch%today_nppwood     (ico) * area_scale
            cpatch%today_nppdaily       (ico) = cpatch%today_nppdaily    (ico) * area_scale
            cpatch%today_gpp_pot        (ico) = cpatch%today_gpp_pot     (ico) * area_scale
            cpatch%today_gpp_max        (ico) = cpatch%today_gpp_max     (ico) * area_scale
            cpatch%today_leaf_resp      (ico) = cpatch%today_leaf_resp   (ico) * area_scale
            cpatch%today_root_resp      (ico) = cpatch%today_root_resp   (ico) * area_scale
                     
            !----- Crown area shall not exceed one. ---------------------------------------!
            cpatch%crown_area           (ico) = min(1.,cpatch%crown_area (ico) * area_scale)
            if (idoutput > 0 .or. imoutput > 0 .or. iqoutput > 0) then
               cpatch%dmean_par_l       (ico) = cpatch%dmean_par_l       (ico) * area_scale
               cpatch%dmean_par_l_beam  (ico) = cpatch%dmean_par_l_beam  (ico) * area_scale
               cpatch%dmean_par_l_diff  (ico) = cpatch%dmean_par_l_diff  (ico) * area_scale
            end if
            if (imoutput > 0 .or. iqoutput > 0) then
               cpatch%mmean_par_l       (ico) = cpatch%mmean_par_l       (ico) * area_scale
               cpatch%mmean_par_l_beam  (ico) = cpatch%mmean_par_l_beam  (ico) * area_scale
               cpatch%mmean_par_l_diff  (ico) = cpatch%mmean_par_l_diff  (ico) * area_scale
            end if
            if (iqoutput > 0) then
               cpatch%qmean_par_l     (:,ico) = cpatch%qmean_par_l     (:,ico) * area_scale
               cpatch%qmean_par_l_beam(:,ico) = cpatch%qmean_par_l_beam(:,ico) * area_scale
               cpatch%qmean_par_l_diff(:,ico) = cpatch%qmean_par_l_diff(:,ico) * area_scale
            end if
         end do
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
   !   This subroutine will perform cohort fusion based on various similarity criteria to  !
   ! determine whether they can be fused with no significant loss of information. The user !
   ! is welcome to set up a benchmark, but should be aware that no miracles will happen    !
   ! here. If there are more very distinct cohorts than maxcohort, then the user will need !
   ! to live with that and accept life is not always fair with those with limited          !
   ! computational resources.                                                              !
   !---------------------------------------------------------------------------------------!

   subroutine fuse_cohorts(csite,ipa, green_leaf_factor, lsl)

      use ed_state_vars       , only : sitetype            & ! Structure
                                     , patchtype           ! ! Structure
      use pft_coms            , only : rho                 & ! intent(in)
                                     , b1Ht                & ! intent(in)
                                     , hgt_max             & ! intent(in)
                                     , sla                 & ! intent(in)
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
                                     , dbh2bl              ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target      :: csite             ! Current site
      integer                , intent(in)  :: ipa               ! Current patch ID
      real, dimension(n_pft) , intent(in)  :: green_leaf_factor ! 
      integer                , intent(in)  :: lsl               ! Lowest soil level
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
                  lai_max = ( cpatch%nplant(recc)                                          &
                            * dbh2bl(cpatch%dbh(recc),cpatch%pft(recc))                    &
                            + cpatch%nplant(donc)                                          &
                            * dbh2bl(cpatch%dbh(donc),cpatch%pft(donc)))                   &
                          * cpatch%sla(recc)

                  !----- Checking the total size of this cohort before and after fusion. --!
                  total_size = cpatch%nplant(donc) * ( cpatch%balive(donc)                 &
                                                     + cpatch%bdead(donc)                  &
                                                     + cpatch%bstorage(donc) )             &
                             + cpatch%nplant(recc) * ( cpatch%balive(recc)                 &
                                                     + cpatch%bdead(recc)                  &
                                                     + cpatch%bstorage(recc) )


                  !------------------------------------------------------------------------!
                  !    Five conditions must be met to allow two cohorts to be fused:       !
                  ! 1. Both cohorts must have the same PFT;                                !
                  ! 2. Combined LAI won't be too large.                                    !
                  ! 3. Both cohorts must have the same status with respect to the first    !
                  !    census.                                                             !
                  ! 4. Both cohorts must have the same recruit status with respect to the  !
                  !    first census.                                                       !
                  ! 5. Both cohorts must have the same phenology status.                   !
                  !------------------------------------------------------------------------!
                  if (     cpatch%pft(donc)              == cpatch%pft(recc)               &
                     .and. lai_max                        < lai_fuse_tol*tolerance_mult    &
                     .and. cpatch%first_census(donc)     == cpatch%first_census(recc)      &
                     .and. cpatch%new_recruit_flag(donc) == cpatch%new_recruit_flag(recc)  &
                     .and. cpatch%phenology_status(donc) == cpatch%phenology_status(recc)  &
                     ) then

                     !----- Proceed with fusion -------------------------------------------!
                     call fuse_2_cohorts(cpatch,donc,recc,newn                             &
                                        ,green_leaf_factor(cpatch%pft(donc))               &
                                        ,csite%can_prss(ipa),lsl)

                     !----- Flag donating cohort as gone, so it won't be checked again. ---!
                     fuse_table(donc) = .false.
                     
                     !----- Checking whether total size and LAI are conserved. ------------!
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

      use ed_state_vars        , only : patchtype              ! ! structure
      use pft_coms             , only : q                      & ! intent(in), lookup table
                                      , qsw                    ! ! intent(in), lookup table
      use fusion_fission_coms  , only : lai_tol                ! ! intent(in)
      use ed_max_dims          , only : n_pft                  ! ! intent(in)
      use allometry            , only : dbh2h                  & ! function
                                      , bd2dbh                 & ! function
                                      , dbh2bd                 ! ! function
      use ed_misc_coms         , only : iqoutput               & ! intent(in)
                                      , imoutput               & ! intent(in)
                                      , idoutput               ! ! intent(in)
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
      integer                              :: ipa,ico,inew      ! Counters
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
               ! need to be rescaled.                                                      !
               ! IMPORTANT: Only cohort-level variables that have units per area           !
               !            (m2_ground) should be rescaled.  Variables whose units are per !
               !            plant or per leaf area (m2_leaf) should _NOT_ be included      !
               !            here.                                                          !
               !---------------------------------------------------------------------------!
               cpatch%lai                  (ico) = cpatch%lai               (ico) * 0.5
               cpatch%wpa                  (ico) = cpatch%wpa               (ico) * 0.5
               cpatch%wai                  (ico) = cpatch%wai               (ico) * 0.5
               cpatch%crown_area           (ico) = cpatch%crown_area        (ico) * 0.5
               cpatch%nplant               (ico) = cpatch%nplant            (ico) * 0.5
               cpatch%mean_gpp             (ico) = cpatch%mean_gpp          (ico) * 0.5
               cpatch%mean_leaf_resp       (ico) = cpatch%mean_leaf_resp    (ico) * 0.5
               cpatch%mean_root_resp       (ico) = cpatch%mean_root_resp    (ico) * 0.5
               cpatch%mean_growth_resp     (ico) = cpatch%mean_growth_resp  (ico) * 0.5
               cpatch%mean_storage_resp    (ico) = cpatch%mean_storage_resp (ico) * 0.5
               cpatch%mean_vleaf_resp      (ico) = cpatch%mean_vleaf_resp   (ico) * 0.5
               cpatch%today_gpp            (ico) = cpatch%today_gpp         (ico) * 0.5
               cpatch%today_nppleaf        (ico) = cpatch%today_nppleaf     (ico) * 0.5
               cpatch%today_nppfroot       (ico) = cpatch%today_nppfroot    (ico) * 0.5
               cpatch%today_nppsapwood     (ico) = cpatch%today_nppsapwood  (ico) * 0.5
               cpatch%today_nppcroot       (ico) = cpatch%today_nppcroot    (ico) * 0.5
               cpatch%today_nppseeds       (ico) = cpatch%today_nppseeds    (ico) * 0.5
               cpatch%today_nppdaily       (ico) = cpatch%today_nppdaily    (ico) * 0.5
               cpatch%today_gpp_pot        (ico) = cpatch%today_gpp_pot     (ico) * 0.5
               cpatch%today_gpp_max        (ico) = cpatch%today_gpp_max     (ico) * 0.5
               cpatch%today_leaf_resp      (ico) = cpatch%today_leaf_resp   (ico) * 0.5
               cpatch%today_root_resp      (ico) = cpatch%today_root_resp   (ico) * 0.5
               cpatch%gpp                  (ico) = cpatch%gpp               (ico) * 0.5
               cpatch%leaf_respiration     (ico) = cpatch%leaf_respiration  (ico) * 0.5
               cpatch%root_respiration     (ico) = cpatch%root_respiration  (ico) * 0.5
               cpatch%monthly_dndt         (ico) = cpatch%monthly_dndt      (ico) * 0.5
               cpatch%leaf_water           (ico) = cpatch%leaf_water        (ico) * 0.5
               cpatch%leaf_hcap            (ico) = cpatch%leaf_hcap         (ico) * 0.5
               cpatch%leaf_energy          (ico) = cpatch%leaf_energy       (ico) * 0.5
               cpatch%wood_water           (ico) = cpatch%wood_water        (ico) * 0.5
               cpatch%wood_hcap            (ico) = cpatch%wood_hcap         (ico) * 0.5
               cpatch%wood_energy          (ico) = cpatch%wood_energy       (ico) * 0.5
               if (idoutput > 0 .or. imoutput > 0 .or. iqoutput > 0 ) then
                  cpatch%dmean_par_l       (ico) = cpatch%dmean_par_l     (ico)   * 0.5
                  cpatch%dmean_par_l_beam  (ico) = cpatch%dmean_par_l_beam(ico)   * 0.5
                  cpatch%dmean_par_l_diff  (ico) = cpatch%dmean_par_l_diff(ico)   * 0.5
               end if
               if (imoutput > 0 .or. iqoutput > 0  ) then
                  cpatch%mmean_par_l       (ico) = cpatch%mmean_par_l     (ico)   * 0.5
                  cpatch%mmean_par_l_beam  (ico) = cpatch%mmean_par_l_beam(ico)   * 0.5
                  cpatch%mmean_par_l_diff  (ico) = cpatch%mmean_par_l_diff(ico)   * 0.5
               end if
               if (iqoutput > 0  ) then
                  cpatch%qmean_par_l     (:,ico) = cpatch%qmean_par_l     (:,ico) * 0.5
                  cpatch%qmean_par_l_beam(:,ico) = cpatch%qmean_par_l_beam(:,ico) * 0.5
                  cpatch%qmean_par_l_diff(:,ico) = cpatch%qmean_par_l_diff(:,ico) * 0.5
               end if

               !---------------------------------------------------------------------------!


               !----- Apply these values to the new cohort. -------------------------------!
               inew = inew+1
               call clone_cohort(cpatch,ico,inew)
               !---------------------------------------------------------------------------!

               !----- Tweaking bdead, to ensure carbon is conserved. ----------------------!
               cpatch%bdead(ico)  = cpatch%bdead(ico) * (1.-epsilon)
               cpatch%dbh  (ico)  = bd2dbh(cpatch%pft(ico), cpatch%bdead(ico))
               cpatch%hite (ico)  = dbh2h(cpatch%pft(ico), cpatch%dbh(ico))

               cpatch%bdead(inew) = cpatch%bdead(inew) * (1.+epsilon)
               cpatch%dbh  (inew) = bd2dbh(cpatch%pft(inew), cpatch%bdead(inew))
               cpatch%hite (inew) = dbh2h(cpatch%pft(inew), cpatch%dbh(inew))
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
   !   This subroutine will clone one cohort.                                              !
   !---------------------------------------------------------------------------------------!
   subroutine clone_cohort(cpatch,isc,idt)
   
      use ed_max_dims  , only : n_mort     ! ! intent(in)
      use ed_state_vars, only : patchtype  & ! Structure
                              , stoma_data ! ! Structure
      use ed_misc_coms , only : iqoutput   & ! intent(in)
                              , idoutput   & ! intent(in)
                              , imoutput   ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype) , target     :: cpatch ! Current patch
      integer         , intent(in) :: isc    ! Index of "Source" cohort
      integer         , intent(in) :: idt    ! Index of "Destination" cohort"
      !----- Local variables --------------------------------------------------------------!
      integer                      :: imonth
      type(stoma_data), pointer    :: osdt,ossc
      !------------------------------------------------------------------------------------!

      cpatch%pft(idt)                  = cpatch%pft(isc)
      cpatch%nplant(idt)               = cpatch%nplant(isc)
      cpatch%hite(idt)                 = cpatch%hite(isc)
      cpatch%dbh(idt)                  = cpatch%dbh(isc)
      cpatch%bdead(idt)                = cpatch%bdead(isc)
      cpatch%bleaf(idt)                = cpatch%bleaf(isc)
      cpatch%broot(idt)                = cpatch%broot(isc)
      cpatch%bsapwood(idt)             = cpatch%bsapwood(isc)
      cpatch%phenology_status(idt)     = cpatch%phenology_status(isc)
      cpatch%balive(idt)               = cpatch%balive(isc)
      cpatch%lai(idt)                  = cpatch%lai(isc)
      cpatch%wpa(idt)                  = cpatch%wpa(isc)
      cpatch%wai(idt)                  = cpatch%wai(isc)
      cpatch%crown_area(idt)           = cpatch%crown_area(isc)
      cpatch%bstorage(idt)             = cpatch%bstorage(isc)
      cpatch%leaf_resolvable(idt)      = cpatch%leaf_resolvable(isc)
      cpatch%wood_resolvable(idt)      = cpatch%wood_resolvable(isc)

      do imonth = 1,13
         cpatch%cb(imonth,idt)         = cpatch%cb(imonth,isc)
         cpatch%cb_max(imonth,idt)     = cpatch%cb_max(imonth,isc)
      enddo

      cpatch%cbr_bar(idt)              = cpatch%cbr_bar(isc)
      cpatch%leaf_energy(idt)          = cpatch%leaf_energy(isc)
      cpatch%leaf_hcap(idt)            = cpatch%leaf_hcap(isc)
      cpatch%leaf_temp(idt)            = cpatch%leaf_temp(isc)
      cpatch%leaf_temp_pv(idt)         = cpatch%leaf_temp_pv(isc)
      cpatch%leaf_fliq(idt)            = cpatch%leaf_fliq(isc)
      cpatch%leaf_water(idt)           = cpatch%leaf_water(isc)
      cpatch%wood_energy(idt)          = cpatch%wood_energy(isc)
      cpatch%wood_hcap(idt)            = cpatch%wood_hcap(isc)
      cpatch%wood_temp(idt)            = cpatch%wood_temp(isc)
      cpatch%wood_temp_pv(idt)         = cpatch%wood_temp_pv(isc)
      cpatch%wood_fliq(idt)            = cpatch%wood_fliq(isc)
      cpatch%wood_water(idt)           = cpatch%wood_water(isc)
      cpatch%veg_wind(idt)             = cpatch%veg_wind(isc)
      cpatch%lsfc_shv_open(idt)        = cpatch%lsfc_shv_open(isc)
      cpatch%lsfc_shv_closed(idt)      = cpatch%lsfc_shv_closed(isc)
      cpatch%lsfc_co2_open(idt)        = cpatch%lsfc_co2_open(isc)
      cpatch%lsfc_co2_closed(idt)      = cpatch%lsfc_co2_closed(isc)
      cpatch%lint_shv(idt)             = cpatch%lint_shv(isc)
      cpatch%lint_co2_open(idt)        = cpatch%lint_co2_open(isc)
      cpatch%lint_co2_closed(idt)      = cpatch%lint_co2_closed(isc)
      cpatch%mean_gpp(idt)             = cpatch%mean_gpp(isc)
      cpatch%mean_leaf_resp(idt)       = cpatch%mean_leaf_resp(isc)
      cpatch%mean_root_resp(idt)       = cpatch%mean_root_resp(isc)
      cpatch%mean_storage_resp(idt)    = cpatch%mean_storage_resp(isc)
      cpatch%mean_growth_resp(idt)     = cpatch%mean_growth_resp(isc)
      cpatch%mean_vleaf_resp(idt)      = cpatch%mean_vleaf_resp(isc)
      cpatch%today_leaf_resp(idt)      = cpatch%today_leaf_resp(isc)
      cpatch%today_root_resp(idt)      = cpatch%today_root_resp(isc)
      cpatch%today_gpp(idt)            = cpatch%today_gpp(isc)
      cpatch%today_nppleaf(idt)        = cpatch%today_nppleaf(isc)
      cpatch%today_nppfroot(idt)       = cpatch%today_nppfroot(isc)
      cpatch%today_nppsapwood(idt)     = cpatch%today_nppsapwood(isc)
      cpatch%today_nppcroot(idt)       = cpatch%today_nppcroot(isc)
      cpatch%today_nppseeds(idt)       = cpatch%today_nppseeds(isc)
      cpatch%today_nppwood(idt)        = cpatch%today_nppwood(isc)
      cpatch%today_nppdaily(idt)       = cpatch%today_nppdaily(isc)
      cpatch%today_gpp_pot(idt)        = cpatch%today_gpp_pot(isc)
      cpatch%today_gpp_max(idt)        = cpatch%today_gpp_max(isc)
      cpatch%growth_respiration(idt)   = cpatch%growth_respiration(isc)
      cpatch%storage_respiration(idt)  = cpatch%storage_respiration(isc)
      cpatch%vleaf_respiration(idt)    = cpatch%vleaf_respiration(isc)
      cpatch%fsn(idt)                  = cpatch%fsn(isc)
      cpatch%monthly_dndt(idt)         = cpatch%monthly_dndt(isc)
      cpatch%agb(idt)                  = cpatch%agb(isc)
      cpatch%basarea(idt)              = cpatch%basarea(isc)
      cpatch%dagb_dt(idt)              = cpatch%dagb_dt(isc)
      cpatch%dba_dt(idt)               = cpatch%dba_dt(isc)
      cpatch%ddbh_dt(idt)              = cpatch%ddbh_dt(isc)
      cpatch%Psi_open(idt)             = cpatch%Psi_open(isc)
      cpatch%krdepth(idt)              = cpatch%krdepth(isc)
      cpatch%first_census(idt)         = cpatch%first_census(isc)
      cpatch%new_recruit_flag(idt)     = cpatch%new_recruit_flag(isc)
      cpatch%par_l(idt)                = cpatch%par_l(isc)
      cpatch%par_l_beam(idt)           = cpatch%par_l_beam(isc)
      cpatch%par_l_diffuse(idt)        = cpatch%par_l_diffuse(isc)
      cpatch%rshort_l(idt)             = cpatch%rshort_l(isc)
      cpatch%rshort_l_beam(idt)        = cpatch%rshort_l_beam(isc)
      cpatch%rshort_l_diffuse(idt)     = cpatch%rshort_l_diffuse(isc)
      cpatch%rlong_l(idt)              = cpatch%rlong_l(isc)
      cpatch%rlong_l_surf(idt)         = cpatch%rlong_l_surf(isc)
      cpatch%rlong_l_incid(idt)        = cpatch%rlong_l_incid(isc)
      cpatch%rshort_w(idt)             = cpatch%rshort_w(isc)
      cpatch%rshort_w_beam(idt)        = cpatch%rshort_w_beam(isc)
      cpatch%rshort_w_diffuse(idt)     = cpatch%rshort_w_diffuse(isc)
      cpatch%rlong_w(idt)              = cpatch%rlong_w(isc)
      cpatch%rlong_w_surf(idt)         = cpatch%rlong_w_surf(isc)
      cpatch%rlong_w_incid(idt)        = cpatch%rlong_w_incid(isc)
      cpatch%light_level(idt)          = cpatch%light_level(isc)
      cpatch%light_level_beam(idt)     = cpatch%light_level_beam(isc)
      cpatch%light_level_diff(idt)     = cpatch%light_level_diff(isc)
      cpatch%lambda_light(idt)         = cpatch%lambda_light(isc)
      cpatch%beamext_level(idt)        = cpatch%beamext_level(isc)
      cpatch%diffext_level(idt)        = cpatch%diffext_level(isc)
      cpatch%leaf_gbh(idt)             = cpatch%leaf_gbh(isc)
      cpatch%leaf_gbw(idt)             = cpatch%leaf_gbw(isc)
      cpatch%wood_gbh(idt)             = cpatch%wood_gbh(isc)
      cpatch%wood_gbw(idt)             = cpatch%wood_gbw(isc)
      cpatch%A_open(idt)               = cpatch%A_open(isc)
      cpatch%A_closed(idt)             = cpatch%A_closed(isc)
      cpatch%Psi_closed(idt)           = cpatch%Psi_closed(isc)
      cpatch%gsw_open(idt)             = cpatch%gsw_open(isc)
      cpatch%gsw_closed(idt)           = cpatch%gsw_closed(isc)
      cpatch%fsw(idt)                  = cpatch%fsw(isc)
      cpatch%fs_open(idt)              = cpatch%fs_open(isc)
      cpatch%water_supply(idt)         = cpatch%water_supply(isc)
      cpatch%stomatal_conductance(idt) = cpatch%stomatal_conductance(isc)
      cpatch%leaf_maintenance(idt)     = cpatch%leaf_maintenance(isc)
      cpatch%root_maintenance(idt)     = cpatch%root_maintenance(isc)
      cpatch%leaf_drop(idt)            = cpatch%leaf_drop(isc)
      cpatch%bseeds(idt)               = cpatch%bseeds(isc)
      cpatch%leaf_respiration(idt)     = cpatch%leaf_respiration(isc)
      cpatch%root_respiration(idt)     = cpatch%root_respiration(isc)
      cpatch%mort_rate(:,idt)          = cpatch%mort_rate(:,isc)

      cpatch%gpp(idt)                  = cpatch%gpp(isc)
      cpatch%paw_avg(idt)              = cpatch%paw_avg(isc)
      cpatch%elongf(idt)               = cpatch%elongf(isc)

      cpatch%turnover_amp(idt)         = cpatch%turnover_amp(isc)     
      cpatch%llspan(idt)               = cpatch%llspan(isc)     
      cpatch%vm_bar(idt)               = cpatch%vm_bar(isc)  
      cpatch%sla(idt)                  = cpatch%sla(isc)  

      cpatch%old_stoma_vector(:,idt)   = cpatch%old_stoma_vector(:,isc)

      osdt => cpatch%old_stoma_data(idt)
      ossc => cpatch%old_stoma_data(isc)

      osdt%recalc           = ossc%recalc
      osdt%T_L              = ossc%T_L
      osdt%e_A              = ossc%e_A
      osdt%PAR              = ossc%PAR
      osdt%rb_factor        = ossc%rb_factor
      osdt%prss             = ossc%prss
      osdt%phenology_factor = ossc%phenology_factor
      osdt%gsw_open         = ossc%gsw_open
      osdt%ilimit           = ossc%ilimit
      osdt%T_L_residual     = ossc%T_L_residual
      osdt%e_a_residual     = ossc%e_a_residual
      osdt%par_residual     = ossc%par_residual
      osdt%rb_residual      = ossc%rb_residual
      osdt%leaf_residual    = ossc%leaf_residual
      osdt%gsw_residual     = ossc%gsw_residual
     
     
      if (idoutput > 0 .or. imoutput > 0 .or. iqoutput > 0) then
         cpatch%dmean_par_l           (idt) = cpatch%dmean_par_l           (isc) 
         cpatch%dmean_par_l_beam      (idt) = cpatch%dmean_par_l_beam      (isc) 
         cpatch%dmean_par_l_diff      (idt) = cpatch%dmean_par_l_diff      (isc) 
         cpatch%dmean_gpp             (idt) = cpatch%dmean_gpp             (isc)
         cpatch%dmean_nppleaf         (idt) = cpatch%dmean_nppleaf         (isc)
         cpatch%dmean_nppfroot        (idt) = cpatch%dmean_nppfroot        (isc)
         cpatch%dmean_nppsapwood      (idt) = cpatch%dmean_nppsapwood      (isc)
         cpatch%dmean_nppcroot        (idt) = cpatch%dmean_nppcroot        (isc)
         cpatch%dmean_nppseeds        (idt) = cpatch%dmean_nppseeds        (isc)
         cpatch%dmean_nppwood         (idt) = cpatch%dmean_nppwood         (isc)
         cpatch%dmean_nppdaily        (idt) = cpatch%dmean_nppdaily        (isc)
         cpatch%dmean_leaf_resp       (idt) = cpatch%dmean_leaf_resp       (isc)
         cpatch%dmean_root_resp       (idt) = cpatch%dmean_root_resp       (isc)
         cpatch%dmean_fs_open         (idt) = cpatch%dmean_fs_open         (isc)
         cpatch%dmean_fsw             (idt) = cpatch%dmean_fsw             (isc)
         cpatch%dmean_fsn             (idt) = cpatch%dmean_fsn             (isc)
         cpatch%dmean_psi_open        (idt) = cpatch%dmean_psi_open        (isc)
         cpatch%dmean_psi_closed      (idt) = cpatch%dmean_psi_closed      (isc)
         cpatch%dmean_water_supply    (idt) = cpatch%dmean_water_supply    (isc)
         cpatch%dmean_lambda_light    (idt) = cpatch%dmean_lambda_light    (isc)
         cpatch%dmean_light_level     (idt) = cpatch%dmean_light_level     (isc)
         cpatch%dmean_light_level_beam(idt) = cpatch%dmean_light_level_beam(isc)
         cpatch%dmean_light_level_diff(idt) = cpatch%dmean_light_level_diff(isc)
         cpatch%dmean_beamext_level   (idt) = cpatch%dmean_beamext_level   (isc)
         cpatch%dmean_diffext_level   (idt) = cpatch%dmean_diffext_level   (isc)
      end if

      if (imoutput > 0 .or. iqoutput > 0) then
         cpatch%mmean_par_l             (idt) = cpatch%mmean_par_l             (isc) 
         cpatch%mmean_par_l_beam        (idt) = cpatch%mmean_par_l_beam        (isc) 
         cpatch%mmean_par_l_diff        (idt) = cpatch%mmean_par_l_diff        (isc) 
         cpatch%mmean_fs_open           (idt) = cpatch%mmean_fs_open           (isc)
         cpatch%mmean_fsw               (idt) = cpatch%mmean_fsw               (isc)
         cpatch%mmean_fsn               (idt) = cpatch%mmean_fsn               (isc)
         cpatch%mmean_psi_open          (idt) = cpatch%mmean_psi_open          (isc)
         cpatch%mmean_psi_closed        (idt) = cpatch%mmean_psi_closed        (isc)
         cpatch%mmean_water_supply      (idt) = cpatch%mmean_water_supply      (isc)
         cpatch%mmean_leaf_maintenance  (idt) = cpatch%mmean_leaf_maintenance  (isc)
         cpatch%mmean_root_maintenance  (idt) = cpatch%mmean_root_maintenance  (isc)
         cpatch%mmean_leaf_drop         (idt) = cpatch%mmean_leaf_drop         (isc)
         cpatch%mmean_cb                (idt) = cpatch%mmean_cb                (isc)
         cpatch%mmean_lambda_light      (idt) = cpatch%mmean_lambda_light      (isc)
         cpatch%mmean_light_level       (idt) = cpatch%mmean_light_level       (isc)
         cpatch%mmean_light_level_beam  (idt) = cpatch%mmean_light_level_beam  (isc)
         cpatch%mmean_light_level_diff  (idt) = cpatch%mmean_light_level_diff  (isc)
         cpatch%mmean_beamext_level     (idt) = cpatch%mmean_beamext_level     (isc)
         cpatch%mmean_diffext_level     (idt) = cpatch%mmean_diffext_level     (isc)
         cpatch%mmean_gpp               (idt) = cpatch%mmean_gpp               (isc)
         cpatch%mmean_nppleaf           (idt) = cpatch%mmean_nppleaf           (isc)
         cpatch%mmean_nppfroot          (idt) = cpatch%mmean_nppfroot          (isc)
         cpatch%mmean_nppsapwood        (idt) = cpatch%mmean_nppsapwood        (isc)
         cpatch%mmean_nppcroot          (idt) = cpatch%mmean_nppcroot          (isc)
         cpatch%mmean_nppseeds          (idt) = cpatch%mmean_nppseeds          (isc)
         cpatch%mmean_nppwood           (idt) = cpatch%mmean_nppwood           (isc)
         cpatch%mmean_nppdaily          (idt) = cpatch%mmean_nppdaily          (isc)
         
         cpatch%mmean_leaf_resp         (idt) = cpatch%mmean_leaf_resp         (isc)
         cpatch%mmean_root_resp         (idt) = cpatch%mmean_root_resp         (isc)
         cpatch%mmean_growth_resp       (idt) = cpatch%mmean_growth_resp       (isc)
         cpatch%mmean_storage_resp      (idt) = cpatch%mmean_storage_resp      (isc)
         cpatch%mmean_vleaf_resp        (idt) = cpatch%mmean_vleaf_resp        (isc)
         cpatch%mmean_mort_rate       (:,idt) = cpatch%mmean_mort_rate       (:,isc)
      end if

      if (iqoutput > 0) then
         cpatch%qmean_par_l        (:,idt) = cpatch%qmean_par_l        (:,isc)
         cpatch%qmean_par_l_beam   (:,idt) = cpatch%qmean_par_l_beam   (:,isc)
         cpatch%qmean_par_l_diff   (:,idt) = cpatch%qmean_par_l_diff   (:,isc)
         cpatch%qmean_fs_open      (:,idt) = cpatch%qmean_fs_open      (:,isc)
         cpatch%qmean_fsw          (:,idt) = cpatch%qmean_fsw          (:,isc)
         cpatch%qmean_fsn          (:,idt) = cpatch%qmean_fsn          (:,isc)
         cpatch%qmean_psi_open     (:,idt) = cpatch%qmean_psi_open     (:,isc)
         cpatch%qmean_psi_closed   (:,idt) = cpatch%qmean_psi_closed   (:,isc)
         cpatch%qmean_water_supply (:,idt) = cpatch%qmean_water_supply (:,isc)
         cpatch%qmean_gpp          (:,idt) = cpatch%qmean_gpp          (:,isc)
         cpatch%qmean_leaf_resp    (:,idt) = cpatch%qmean_leaf_resp    (:,isc)
         cpatch%qmean_root_resp    (:,idt) = cpatch%qmean_root_resp    (:,isc)
      end if

      return
   end subroutine clone_cohort
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will merge two cohorts into 1. The donating cohort (donc) is the  !
   ! one that will be deallocated, while the receptor cohort (recc) will contain the       !
   !  information from both cohorts.                                                       !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_cohorts(cpatch,donc,recc, newn,green_leaf_factor, can_prss,lsl)
      use ed_state_vars , only : patchtype              ! ! Structure
      use pft_coms      , only : q                      & ! intent(in), lookup table
                               , qsw                    ! ! intent(in), lookup table
      use therm_lib     , only : uextcm2tl              & ! subroutine
                               , qslif                  ! ! function
      use allometry     , only : dbh2krdepth            & ! function
                               , bd2dbh                 & ! function
                               , dbh2h                  ! ! function
      use ed_max_dims   , only : n_mort                 ! ! intent(in)
      use ed_misc_coms  , only : imoutput               & ! intent(in)
                               , iqoutput               & ! intent(in)
                               , idoutput               & ! intent(in)
                               , ndcycle                ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype) , target     :: cpatch            ! Current patch
      integer                      :: donc              ! Donating cohort.
      integer                      :: recc              ! Receptor cohort.
      real            , intent(in) :: newn              ! New nplant
      real            , intent(in) :: green_leaf_factor ! Green leaf factor
      real            , intent(in) :: can_prss          ! Canopy air pressure
      integer         , intent(in) :: lsl               ! Lowest soil level
      !----- Local variables --------------------------------------------------------------!
      integer                      :: imon              ! Month for cb loop
      integer                      :: icyc              ! Time of day for dcycle loop
      integer                      :: imty              ! Mortality type
      real                         :: newni             ! Inverse of new nplants
      real                         :: newlaii           ! Inverse of new LAI
      real                         :: cb_act            !
      real                         :: cb_max            !
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Find the scaling factor for variables that are not "extensive".                 !
      !  - If the unit is X/plant, then we scale by nplant.                                !
      !  - If the unit is X/m2_leaf, then we scale by LAI.                                 !
      !  - If the unit is X/m2_gnd, then we add, since they are "extensive".               !
      !------------------------------------------------------------------------------------!
      newni   = 1.0 / newn
      if (cpatch%lai(recc) + cpatch%lai(donc) > 0.0) then
         newlaii = 1.0 / (cpatch%lai(recc) + cpatch%lai(donc))
      else
         newlaii = 0.0
      end if
      !------------------------------------------------------------------------------------!



      !----- Conserve carbon by calculating bdead first. ----------------------------------!
      cpatch%bdead(recc) = ( cpatch%nplant(recc) * cpatch%bdead(recc)                      &
                           + cpatch%nplant(donc) * cpatch%bdead(donc) ) * newni

      !----- Then get dbh and hite from bdead. --------------------------------------------!
      cpatch%dbh(recc)   = bd2dbh(cpatch%pft(recc), cpatch%bdead(recc))
      cpatch%hite(recc)  = dbh2h(cpatch%pft(recc),  cpatch%dbh(recc))


      !------------------------------------------------------------------------------------!
      !     Conserving carbon to get balive, bleaf, and bstorage.                          !
      !------------------------------------------------------------------------------------!
      cpatch%balive(recc)    = ( cpatch%nplant(recc) * cpatch%balive(recc)                 &
                               + cpatch%nplant(donc) * cpatch%balive(donc) ) *newni
      cpatch%broot(recc)     = ( cpatch%nplant(recc) * cpatch%broot(recc)                  &
                               + cpatch%nplant(donc) * cpatch%broot(donc) ) *newni
      cpatch%bsapwood(recc)  = ( cpatch%nplant(recc) * cpatch%bsapwood(recc)               &
                             + cpatch%nplant(donc) * cpatch%bsapwood(donc) ) *newni
      cpatch%bstorage(recc)  = ( cpatch%nplant(recc) * cpatch%bstorage(recc)               &
                               + cpatch%nplant(donc) * cpatch%bstorage(donc) ) * newni
      cpatch%bseeds(recc)    = ( cpatch%nplant(recc) * cpatch%bseeds(recc)                 &
                               + cpatch%nplant(donc) * cpatch%bseeds(donc) ) * newni
      cpatch%leaf_maintenance(recc) = newni                                                &
                            * ( cpatch%nplant(recc) * cpatch%leaf_maintenance(recc)        &
                              + cpatch%nplant(donc) * cpatch%leaf_maintenance(donc) )
      cpatch%root_maintenance(recc) = newni                                                &
                            * ( cpatch%nplant(recc) * cpatch%root_maintenance(recc)        &
                              + cpatch%nplant(donc) * cpatch%root_maintenance(donc) )
      cpatch%leaf_drop(recc) = newni                                                       &
                             * ( cpatch%nplant(recc) * cpatch%leaf_drop(recc)              &
                               + cpatch%nplant(donc) * cpatch%leaf_drop(donc) )  
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Bleaf must be zero if phenology status is 2.  This is probably done correctly   !
      ! throughout the code, but being safe here.                                          !
      !------------------------------------------------------------------------------------!
      if (cpatch%phenology_status(recc) < 2) then
         cpatch%bleaf(recc)  = ( cpatch%nplant(recc) * cpatch%bleaf(recc)                  &
                               + cpatch%nplant(donc) * cpatch%bleaf(donc) ) *newni
      else
         cpatch%bleaf(recc)      = 0.
      end if
      !------------------------------------------------------------------------------------!

      cpatch%wpa        (recc) = cpatch%wpa(recc)         + cpatch%wpa        (donc)
      cpatch%wai        (recc) = cpatch%wai(recc)         + cpatch%wai        (donc)
      cpatch%crown_area (recc) = min(1.,cpatch%crown_area(recc)  + cpatch%crown_area(donc))
      cpatch%leaf_energy(recc) = cpatch%leaf_energy(recc) + cpatch%leaf_energy(donc)
      cpatch%leaf_water (recc) = cpatch%leaf_water (recc) + cpatch%leaf_water (donc)
      cpatch%leaf_hcap  (recc) = cpatch%leaf_hcap  (recc) + cpatch%leaf_hcap  (donc)
      cpatch%wood_energy(recc) = cpatch%wood_energy(recc) + cpatch%wood_energy(donc)
      cpatch%wood_water (recc) = cpatch%wood_water (recc) + cpatch%wood_water (donc)
      cpatch%wood_hcap  (recc) = cpatch%wood_hcap  (recc) + cpatch%wood_hcap  (donc)

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
         cpatch%leaf_temp(recc)  = newni                                                   &
                                 * ( cpatch%leaf_temp(recc)  * cpatch%nplant(recc)         &
                                   + cpatch%leaf_temp(donc)  * cpatch%nplant(donc))
         cpatch%leaf_fliq(recc)  = 0.0
      end if
      
      !----- Simply set the previous time-steps temp as the current
      
      
      if ( cpatch%wood_hcap(recc) > 0. ) then
         !----- Update temperature using the standard thermodynamics. ---------------------!
         call uextcm2tl(cpatch%wood_energy(recc),cpatch%wood_water(recc)                   &
                       ,cpatch%wood_hcap(recc),cpatch%wood_temp(recc)                      &
                       ,cpatch%wood_fliq(recc))
      else 
         !----- Wood temperature cannot be found using uextcm2tl, this is a singularity. --!
         cpatch%wood_temp(recc)  = newni                                                   &
                                 * ( cpatch%wood_temp(recc)  * cpatch%nplant(recc)         &
                                   + cpatch%wood_temp(donc)  * cpatch%nplant(donc))
         cpatch%wood_fliq(recc)  = 0.0
      end if
      
      !----- Set time-steps temp as the current
      cpatch%leaf_temp_pv(recc) = cpatch%leaf_temp(recc)
      cpatch%wood_temp_pv(recc) = cpatch%wood_temp(recc)

      !------ Find the intercellular value assuming saturation. ---------------------------!
      cpatch%lint_shv(recc) = qslif(can_prss,cpatch%leaf_temp(recc))

      cb_act = 0.
      cb_max = 0.
      do imon = 1,12
         cpatch%cb(imon,recc)     = ( cpatch%cb(imon,recc) * cpatch%nplant(recc)           &
                                    + cpatch%cb(imon,donc) * cpatch%nplant(donc) ) * newni

         cpatch%cb_max(imon,recc) = ( cpatch%cb_max(imon,recc) * cpatch%nplant(recc)       &
                                    + cpatch%cb_max(imon,donc) * cpatch%nplant(donc))      &
                                    * newni
         cb_act = cb_act + cpatch%cb(imon,recc)
         cb_max = cb_max + cpatch%cb_max(imon,recc)
      end do
      cpatch%cb(13,recc)     = ( cpatch%cb(13,recc) * cpatch%nplant(recc)                  &
                               + cpatch%cb(13,donc) * cpatch%nplant(donc) ) * newni

      cpatch%cb_max(13,recc) = ( cpatch%cb_max(13,recc) * cpatch%nplant(recc)              &
                               + cpatch%cb_max(13,donc) * cpatch%nplant(donc))             &
                               * newni

      if(cb_max > 0.0)then
         cpatch%cbr_bar(recc) = cb_act / cb_max
      else
         cpatch%cbr_bar(recc) = 0.0
      end if



      !------------------------------------------------------------------------------------!
      !     Updating the mean carbon fluxes. They are fluxes per unit of area, so they     !
      ! should be added, not scaled.                                                       !
      !------------------------------------------------------------------------------------!
      cpatch%mean_gpp(recc) = cpatch%mean_gpp(recc) + cpatch%mean_gpp(donc)

      cpatch%mean_leaf_resp(recc)    = cpatch%mean_leaf_resp(recc)                         &
                                     + cpatch%mean_leaf_resp(donc)
      cpatch%mean_root_resp(recc)    = cpatch%mean_root_resp(recc)                         &
                                     + cpatch%mean_root_resp(donc)
      cpatch%mean_storage_resp(recc) = cpatch%mean_storage_resp(recc)                      &
                                     + cpatch%mean_storage_resp(donc)
      cpatch%mean_growth_resp(recc)  = cpatch%mean_growth_resp(recc)                       &
                                     + cpatch%mean_growth_resp(donc)
      cpatch%mean_vleaf_resp(recc)   = cpatch%mean_vleaf_resp(recc)                        &
                                     + cpatch%mean_vleaf_resp(donc)

       !------------------------------------------------------------------------------------!

      cpatch%today_gpp(recc)     = cpatch%today_gpp(recc)                                  &
                                 + cpatch%today_gpp(donc)
                                 
      cpatch%today_nppleaf(recc) = cpatch%today_nppleaf(recc)                              &
                                 + cpatch%today_nppleaf(donc)
                                 
      cpatch%today_nppfroot(recc)= cpatch%today_nppfroot(recc)                             &
                                 + cpatch%today_nppfroot(donc)
                                 
      cpatch%today_nppsapwood(recc) = cpatch%today_nppsapwood(recc)                        &
                                 + cpatch%today_nppsapwood(donc)
                                 
      cpatch%today_nppcroot(recc)= cpatch%today_nppcroot(recc)                             &
                                 + cpatch%today_nppcroot(donc)
                                 
      cpatch%today_nppseeds(recc)= cpatch%today_nppseeds(recc)                             &
                                 + cpatch%today_nppseeds(donc)
                                 
      cpatch%today_nppwood(recc) = cpatch%today_nppwood(recc)                              &
                                 + cpatch%today_nppwood(donc)
                                 
      cpatch%today_nppdaily(recc)= cpatch%today_nppdaily(recc)                             &
                                 + cpatch%today_nppdaily(donc)
                                 
      cpatch%today_nppleaf(recc) = cpatch%today_nppleaf(recc)                              &
                                 + cpatch%today_nppleaf(donc)
                                 
      cpatch%today_gpp_pot(recc) = cpatch%today_gpp_pot(recc)                              &
                                 + cpatch%today_gpp_pot(donc)

      cpatch%today_gpp_max(recc) = cpatch%today_gpp_max(recc)                              &
                                 + cpatch%today_gpp_max(donc)

      cpatch%today_leaf_resp(recc) = cpatch%today_leaf_resp(recc)                          &
                                  + cpatch%today_leaf_resp(donc)

      cpatch%today_root_resp(recc) = cpatch%today_root_resp(recc)                          &
                                   + cpatch%today_root_resp(donc)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Fuse the leaf surface and intenal properties.  Since they are intensive         !
      ! properties, they are scaled by the number of plants.  These numbers are diagnostic !
      ! and this should be used for the output only.                                       !
      !------------------------------------------------------------------------------------!
      cpatch%lsfc_shv_open(recc) = ( cpatch%lsfc_shv_open(recc) * cpatch%nplant(recc)      &
                                   + cpatch%lsfc_shv_open(donc) * cpatch%nplant(donc) )    &
                                   * newni
      cpatch%lsfc_shv_closed(recc) = ( cpatch%lsfc_shv_closed(recc) * cpatch%nplant(recc)  &
                                     + cpatch%lsfc_shv_closed(donc) * cpatch%nplant(donc)) &
                                   * newni
      cpatch%lsfc_co2_open(recc) = ( cpatch%lsfc_co2_open(recc) * cpatch%nplant(recc)      &
                                   + cpatch%lsfc_co2_open(donc) * cpatch%nplant(donc) )    &
                                   * newni
      cpatch%lsfc_co2_closed(recc) = ( cpatch%lsfc_co2_closed(recc) * cpatch%nplant(recc)  &
                                     + cpatch%lsfc_co2_closed(donc) * cpatch%nplant(donc)) &
                                   * newni
      cpatch%lint_co2_open(recc) = ( cpatch%lint_co2_open(recc) * cpatch%nplant(recc)      &
                                   + cpatch%lint_co2_open(donc) * cpatch%nplant(donc) )    &
                                   * newni
      cpatch%lint_co2_closed(recc) = ( cpatch%lint_co2_closed(recc) * cpatch%nplant(recc)  &
                                     + cpatch%lint_co2_closed(donc) * cpatch%nplant(donc)) &
                                   * newni
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Fusing the mortality rates.  The terms that are PFT-dependent but density-      !
      ! independent should be the same, so it doesn't matter which average we use.  The    !
      ! density-dependent should be averaged using nplant as the relative weight.          !
      !------------------------------------------------------------------------------------!
      do imty=1,n_mort
         cpatch%mort_rate(imty,recc) = ( cpatch%mort_rate(imty,recc) *cpatch%nplant(recc)  &
                                       + cpatch%mort_rate(imty,donc) *cpatch%nplant(donc)) &
                                     * newni
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Light level.  Using the intensive way of fusing.                                !
      !------------------------------------------------------------------------------------!
      cpatch%light_level(recc)      = ( cpatch%light_level(recc) *cpatch%nplant(recc)      &
                                      + cpatch%light_level(donc) *cpatch%nplant(donc) )    &
                                    * newni
      cpatch%light_level_beam(recc) = ( cpatch%light_level_beam(recc) *cpatch%nplant(recc) &
                                      + cpatch%light_level_beam(donc) *cpatch%nplant(donc))&
                                    * newni
      cpatch%light_level_diff(recc) = ( cpatch%light_level_diff(recc) *cpatch%nplant(recc) &
                                      + cpatch%light_level_diff(donc) *cpatch%nplant(donc))&
                                    * newni
      cpatch%beamext_level(recc)    = ( cpatch%beamext_level(recc) *cpatch%nplant(recc)    &
                                      + cpatch%beamext_level(donc) *cpatch%nplant(donc) )  &
                                    * newni
      cpatch%diffext_level(recc)    = ( cpatch%diffext_level(recc) *cpatch%nplant(recc)    &
                                      + cpatch%diffext_level(donc) *cpatch%nplant(donc) )  &
                                    * newni
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Not sure about the following variables.  From ed_state_vars, I would say that   !
      ! they should be averaged, not added because there it's written that these are per   !
      ! plant.  But from the fuse_2_patches subroutine here it seems they are per unit  !
      ! area.                                                                              !
      !------------------------------------------------------------------------------------!
      cpatch%growth_respiration(recc)  = newni *                                           &
                                ( cpatch%growth_respiration(recc)  * cpatch%nplant(recc)   &
                                + cpatch%growth_respiration(donc)  * cpatch%nplant(donc) )
     
      cpatch%storage_respiration(recc) = newni *                                           &
                                ( cpatch%storage_respiration(recc) * cpatch%nplant(recc)   &
                                + cpatch%storage_respiration(donc) * cpatch%nplant(donc) )
     
      cpatch%vleaf_respiration(recc)   = newni *                                           &
                                ( cpatch%vleaf_respiration(recc)   * cpatch%nplant(recc)   &
                                + cpatch%vleaf_respiration(donc)   * cpatch%nplant(donc) )
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Water demand is in kg/m2_leaf/s, so we scale them by LAI.  Water supply is in   !
      ! kg/m2_ground/s, so we just add them.                                               !
      !------------------------------------------------------------------------------------!
      cpatch%psi_open    (recc) = ( cpatch%psi_open  (recc) * cpatch%lai(recc)             &
                                  + cpatch%psi_open  (donc) * cpatch%lai(donc) ) * newlaii
      cpatch%psi_closed  (recc) = ( cpatch%psi_closed(recc) * cpatch%lai(recc)             & 
                                  + cpatch%psi_closed(donc) * cpatch%lai(donc) ) * newlaii
      cpatch%water_supply(recc) = cpatch%water_supply(recc) + cpatch%water_supply(donc)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Carbon demand is in kg_C/m2_leaf/s, so we scale them by LAI.  FSW and FSN are   !
      ! really related to leaves, so we scale them by LAI.                                 !
      !------------------------------------------------------------------------------------!
      cpatch%A_open  (recc)     = ( cpatch%A_open  (recc) * cpatch%lai(recc)               &
                                  + cpatch%A_open  (donc) * cpatch%lai(donc) ) * newlaii
      cpatch%A_closed(recc)     = ( cpatch%A_closed(recc) * cpatch%lai(recc)               &
                                  + cpatch%A_closed(donc) * cpatch%lai(donc) ) * newlaii
      cpatch%fsw     (recc)     = ( cpatch%fsw     (recc) * cpatch%lai(recc)               &
                                  + cpatch%fsw     (donc) * cpatch%lai(donc) ) * newlaii
      cpatch%fsn     (recc)     = ( cpatch%fsn     (recc) * cpatch%lai(recc)               &
                                  + cpatch%fsn     (donc) * cpatch%lai(donc) ) * newlaii
      cpatch%fs_open (recc)     = cpatch%fsw(recc) * cpatch%fsn(recc)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Merge biomass and basal area.  Contrary to the patch/site/polygon levels,       !
      ! these variables are "intensive" (or per plant) at the cohort level, so we must     !
      ! average them.                                                                      !
      !------------------------------------------------------------------------------------!
      cpatch%agb(recc)          = ( cpatch%agb(recc)         * cpatch%nplant(recc)         &
                                  + cpatch%agb(donc)         * cpatch%nplant(donc) )       &
                                * newni
      cpatch%basarea(recc)      = ( cpatch%basarea(recc)     * cpatch%nplant(recc)         &
                                  + cpatch%basarea(donc)     * cpatch%nplant(donc) )       &
                                * newni
      cpatch%dagb_dt(recc)      = ( cpatch%dagb_dt(recc)     * cpatch%nplant(recc)         &
                                  + cpatch%dagb_dt(donc)     * cpatch%nplant(donc) )       &
                                * newni
      cpatch%dba_dt(recc)       = ( cpatch%dba_dt(recc)      * cpatch%nplant(recc)         &
                                  + cpatch%dba_dt(donc)      * cpatch%nplant(donc) )       &
                                * newni
      cpatch%ddbh_dt(recc)      = ( cpatch%ddbh_dt(recc)     * cpatch%nplant(recc)         &
                                  + cpatch%ddbh_dt(donc)     * cpatch%nplant(donc) )       &
                                * newni
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Updating the tendency of plant density.  All variables are per unit of area,   !
      ! so they should be added, not scaled.                                               !
      !------------------------------------------------------------------------------------!
      cpatch%monthly_dndt(recc) = cpatch%monthly_dndt(recc) + cpatch%monthly_dndt(donc)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Update the carbon fluxes. They are fluxes per unit of area, so they should be  !
      ! added, not scaled.                                                                 !
      !------------------------------------------------------------------------------------!
      cpatch%gpp(recc) = cpatch%gpp(recc) + cpatch%gpp(donc)

      cpatch%leaf_respiration(recc) = cpatch%leaf_respiration(recc)                        &
                                    + cpatch%leaf_respiration(donc)
      cpatch%root_respiration(recc) = cpatch%root_respiration(recc)                        &
                                    + cpatch%root_respiration(donc)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Potential available water and elongation factor can be consider "intensive"    !
      ! variable water.                                                                    !
      !------------------------------------------------------------------------------------!
      cpatch%paw_avg(recc) = ( cpatch%paw_avg(recc)     * cpatch%nplant(recc)              &
                             + cpatch%paw_avg(donc)     * cpatch%nplant(donc) )            &
                           * newni
      cpatch%elongf(recc)  = ( cpatch%elongf(recc)     * cpatch%nplant(recc)               &
                             + cpatch%elongf(donc)     * cpatch%nplant(donc) )             &
                           * newni
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Light-phenology characteristics (MLO I am not sure if they should be scaled by  !
      ! nplant or LAI, it seems LAI would make more sense...)
      !------------------------------------------------------------------------------------!
      cpatch%turnover_amp(recc)  = ( cpatch%turnover_amp(recc) * cpatch%nplant(recc)       &
                                   + cpatch%turnover_amp(donc) * cpatch%nplant(donc) )     &
                                 * newni

      cpatch%llspan(recc)        = ( cpatch%llspan(recc)       * cpatch%nplant(recc)       &
                                   + cpatch%llspan(donc)       * cpatch%nplant(donc) )     &
                                 * newni

      cpatch%vm_bar(recc)        = ( cpatch%vm_bar(recc) * cpatch%nplant(recc)             &
                                   + cpatch%vm_bar(donc) * cpatch%nplant(donc) )           &
                                 * newni

      cpatch%sla(recc)           = ( cpatch%sla(recc) * cpatch%nplant(recc)                &
                                   + cpatch%sla(donc) * cpatch%nplant(donc) ) * newni
    
      cpatch%krdepth(recc)       = dbh2krdepth(cpatch%hite(recc),cpatch%dbh(recc)          &
                                              ,cpatch%pft(recc),lsl)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Now that we have daily and monthly means going to the cohort level, we must     !
      ! fuse them too.                                                                     !
      !------------------------------------------------------------------------------------!
      if (idoutput > 0 .or. imoutput > 0 .or. iqoutput > 0) then
         cpatch%dmean_light_level       (recc) = ( cpatch%dmean_light_level(recc)          &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_light_level(donc)          &
                                                 * cpatch%nplant(donc) ) * newni
         cpatch%dmean_light_level_beam  (recc) = ( cpatch%dmean_light_level_beam(recc)     &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_light_level_beam(donc)     &
                                                 * cpatch%nplant(donc) ) * newni
         cpatch%dmean_light_level_diff  (recc) = ( cpatch%dmean_light_level_diff(recc)     &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_light_level_diff(donc)     &
                                                 * cpatch%nplant(donc) ) * newni
         cpatch%dmean_beamext_level     (recc) = ( cpatch%dmean_beamext_level(recc)        &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_beamext_level(donc)        &
                                                 * cpatch%nplant(donc) ) * newni
         cpatch%dmean_diffext_level     (recc) = ( cpatch%dmean_diffext_level(recc)        &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_diffext_level(donc)        &
                                                 * cpatch%nplant(donc) ) * newni
         cpatch%dmean_lambda_light      (recc) = ( cpatch%dmean_lambda_light(recc)         &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_lambda_light(donc)         &
                                                 * cpatch%nplant(donc) ) * newni
         cpatch%dmean_gpp               (recc) = ( cpatch%dmean_gpp(recc)                  &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_gpp(donc)                  &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%dmean_nppleaf           (recc) = ( cpatch%dmean_nppleaf(recc)              &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_nppleaf(donc)              &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%dmean_nppfroot          (recc) = ( cpatch%dmean_nppfroot(recc)             &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_nppfroot(donc)             &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%dmean_nppsapwood        (recc) = ( cpatch%dmean_nppsapwood(recc)           &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_nppsapwood(donc)           &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%dmean_nppcroot          (recc) = ( cpatch%dmean_nppcroot(recc)             &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_nppcroot(donc)             &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%dmean_nppseeds          (recc) = ( cpatch%dmean_nppseeds(recc)             &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_nppseeds(donc)             &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%dmean_nppwood           (recc) = ( cpatch%dmean_nppwood(recc)              &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_nppwood(donc)              &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%dmean_nppdaily          (recc) = ( cpatch%dmean_nppdaily(recc)             &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_nppdaily(donc)             &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%dmean_leaf_resp         (recc) = ( cpatch%dmean_leaf_resp(recc)            &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_leaf_resp(donc)            &
                                                 * cpatch%nplant(donc) ) * newni
         cpatch%dmean_root_resp         (recc) = ( cpatch%dmean_root_resp(recc)            &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%dmean_root_resp(donc)            &
                                                 * cpatch%nplant(donc) ) * newni
         !----- The following variables depend on LAI more than nplant. -------------------!
         cpatch%dmean_fs_open           (recc) = ( cpatch%dmean_fs_open(recc)              &
                                                 * cpatch%lai(recc)                        &
                                                 + cpatch%dmean_fs_open(donc)              &
                                                 * cpatch%lai(donc) ) * newlaii
         cpatch%dmean_fsw               (recc) = ( cpatch%dmean_fsw(recc)                  &
                                                 * cpatch%lai(recc)                        &
                                                 + cpatch%dmean_fsw(donc)                  &
                                                 * cpatch%lai(donc) ) * newlaii
         cpatch%dmean_fsn               (recc) = ( cpatch%dmean_fsn(recc)                  &
                                                 * cpatch%lai(recc)                        &
                                                 + cpatch%dmean_fsn(donc)                  &
                                                 * cpatch%lai(donc) ) * newlaii

         !----- The following variables are "extensive", add them. ------------------------!
         cpatch%dmean_par_l             (recc) = cpatch%dmean_par_l       (recc)           &
                                               + cpatch%dmean_par_l       (donc)
         cpatch%dmean_par_l_beam        (recc) = cpatch%dmean_par_l_beam  (recc)           &
                                               + cpatch%dmean_par_l_beam  (donc)
         cpatch%dmean_par_l_diff        (recc) = cpatch%dmean_par_l_diff  (recc)           &
                                               + cpatch%dmean_par_l_diff  (donc)
         cpatch%dmean_psi_open          (recc) = cpatch%dmean_psi_open    (recc)           &
                                               + cpatch%dmean_psi_open    (donc)
         cpatch%dmean_psi_closed        (recc) = cpatch%dmean_psi_closed  (recc)           &
                                               + cpatch%dmean_psi_closed  (donc)
         cpatch%dmean_water_supply      (recc) = cpatch%dmean_water_supply(recc)           &
                                               + cpatch%dmean_water_supply(donc)
      end if
      if (imoutput > 0 .or. iqoutput > 0) then
         cpatch%mmean_light_level     (recc) = ( cpatch%mmean_light_level(recc)            &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_light_level(donc)            &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_light_level_beam(recc) = ( cpatch%mmean_light_level_beam(recc)       &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_light_level_beam(donc)       &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_light_level_diff(recc) = ( cpatch%mmean_light_level_diff(recc)       &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_light_level_diff(donc)       &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_beamext_level   (recc) = ( cpatch%mmean_beamext_level(recc)          &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_beamext_level(donc)          &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_diffext_level   (recc) = ( cpatch%mmean_diffext_level(recc)          &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_diffext_level(donc)          &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_lambda_light    (recc) = ( cpatch%mmean_lambda_light(recc)           &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_lambda_light(donc)           &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_leaf_maintenance(recc) = ( cpatch%mmean_leaf_maintenance(recc)       &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_leaf_maintenance(donc)       &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_root_maintenance(recc) = ( cpatch%mmean_root_maintenance(recc)       &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_root_maintenance(donc)       &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_leaf_drop       (recc) = ( cpatch%mmean_leaf_drop(recc)              &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_leaf_drop(donc)              &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_cb              (recc) = ( cpatch%mmean_cb(recc)                     &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_cb(donc)                     &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_gpp             (recc) = ( cpatch%mmean_gpp(recc)                    &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_gpp(donc)                    &
                                               * cpatch%nplant(donc) ) * newni
                                               
         cpatch%mmean_nppleaf           (recc) = ( cpatch%mmean_nppleaf(recc)              &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%mmean_nppleaf(donc)              &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%mmean_nppfroot          (recc) = ( cpatch%mmean_nppfroot(recc)             &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%mmean_nppfroot(donc)             &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%mmean_nppsapwood        (recc) = ( cpatch%mmean_nppsapwood(recc)           &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%mmean_nppsapwood(donc)           &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%mmean_nppcroot          (recc) = ( cpatch%mmean_nppcroot(recc)             &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%mmean_nppcroot(donc)             &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%mmean_nppseeds          (recc) = ( cpatch%mmean_nppseeds(recc)             &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%mmean_nppseeds(donc)             &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%mmean_nppwood           (recc) = ( cpatch%mmean_nppwood(recc)              &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%mmean_nppwood(donc)              &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%mmean_nppdaily          (recc) = ( cpatch%mmean_nppdaily(recc)             &
                                                 * cpatch%nplant(recc)                     &
                                                 + cpatch%mmean_nppdaily(donc)             &
                                                 * cpatch%nplant(donc) ) * newni
                                                 
         cpatch%mmean_leaf_resp       (recc) = ( cpatch%mmean_leaf_resp(recc)              &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_leaf_resp(donc)              &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_root_resp       (recc) = ( cpatch%mmean_root_resp(recc)              &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_root_resp(donc)              &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_growth_resp     (recc) = ( cpatch%mmean_growth_resp(recc)            &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_growth_resp(donc)            &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_storage_resp    (recc) = ( cpatch%mmean_storage_resp(recc)           &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_storage_resp(donc)           &
                                               * cpatch%nplant(donc) ) * newni
         cpatch%mmean_vleaf_resp      (recc) = ( cpatch%mmean_vleaf_resp(recc)             &
                                               * cpatch%nplant(recc)                       &
                                               + cpatch%mmean_vleaf_resp(donc)             &
                                               * cpatch%nplant(donc) ) * newni

         !---------------------------------------------------------------------------------!
         !    Fusing the mortality rates.  The terms that are PFT-dependent but density-   !
         ! independent should be the same, so it doesn't matter which average we use.  The !
         ! density-dependent should be averaged using nplant as the relative weight.       !
         !---------------------------------------------------------------------------------!
         do imty=1,n_mort
            cpatch%mmean_mort_rate(imty,recc) = ( cpatch%mmean_mort_rate(imty,recc)        &
                                                * cpatch%nplant(recc)                      &
                                                + cpatch%mort_rate(imty,donc)              &
                                                * cpatch%nplant(donc))                     &
                                              * newni
         end do

         !----- The following variables depend on LAI more than nplant. -------------------!
         cpatch%mmean_fs_open           (recc) = ( cpatch%mmean_fs_open(recc)              &
                                                 * cpatch%lai(recc)                        &
                                                 + cpatch%mmean_fs_open(donc)              &
                                                 * cpatch%lai(donc) ) * newlaii
         cpatch%mmean_fsw               (recc) = ( cpatch%mmean_fsw(recc)                  &
                                                 * cpatch%lai(recc)                        &
                                                 + cpatch%mmean_fsw(donc)                  &
                                                 * cpatch%lai(donc) ) * newlaii
         cpatch%mmean_fsn               (recc) = ( cpatch%mmean_fsn(recc)                  &
                                                 * cpatch%lai(recc)                        &
                                                 + cpatch%mmean_fsn(donc)                  &
                                                 * cpatch%lai(donc) ) * newlaii

         !----- The following variables are "extensive", add them. ------------------------!
         cpatch%mmean_par_l             (recc) = cpatch%mmean_par_l       (recc)           &
                                               + cpatch%mmean_par_l       (donc)
         cpatch%mmean_par_l_beam        (recc) = cpatch%mmean_par_l_beam  (recc)           &
                                               + cpatch%mmean_par_l_beam  (donc)
         cpatch%mmean_par_l_diff        (recc) = cpatch%mmean_par_l_diff  (recc)           &
                                               + cpatch%mmean_par_l_diff  (donc)
         cpatch%mmean_psi_open          (recc) = cpatch%mmean_psi_open    (recc)           &
                                               + cpatch%mmean_psi_open    (donc)
         cpatch%mmean_psi_closed        (recc) = cpatch%mmean_psi_closed  (recc)           &
                                               + cpatch%mmean_psi_closed  (donc)
         cpatch%mmean_water_supply      (recc) = cpatch%mmean_water_supply(recc)           &
                                               + cpatch%mmean_water_supply(donc)
      end if

      !------------------------------------------------------------------------------------!
      !    Fuse the mean diurnal cycle.                                                    !
      !------------------------------------------------------------------------------------!
      if (iqoutput > 0) then
         do icyc=1,ndcycle
            cpatch%qmean_gpp          (icyc,recc) = ( cpatch%qmean_gpp        (icyc,recc)  &
                                                    * cpatch%nplant                (recc)  &
                                                    + cpatch%qmean_gpp        (icyc,donc)  &
                                                    * cpatch%nplant                (donc)) &
                                                  * newni
            cpatch%qmean_leaf_resp    (icyc,recc) = ( cpatch%qmean_leaf_resp  (icyc,recc)  &
                                                    * cpatch%nplant                (recc)  &
                                                    + cpatch%qmean_leaf_resp  (icyc,donc)  &
                                                    * cpatch%nplant                (donc)) &
                                                  * newni
            cpatch%qmean_root_resp    (icyc,recc) = ( cpatch%qmean_root_resp  (icyc,recc)  &
                                                    * cpatch%nplant                (recc)  &
                                                    + cpatch%qmean_root_resp  (icyc,donc)  &
                                                    * cpatch%nplant                (donc)) &
                                                  * newni
            !----- The following variables depend on LAI more than nplant. ----------------!
            cpatch%qmean_fs_open      (icyc,recc) = ( cpatch%qmean_fs_open    (icyc,recc)  &
                                                    * cpatch%lai                   (recc)  &
                                                    + cpatch%qmean_fs_open    (icyc,donc)  &
                                                    * cpatch%lai                   (donc)) &
                                                  * newlaii
            cpatch%qmean_fsw          (icyc,recc) = ( cpatch%qmean_fsw        (icyc,recc)  &
                                                    * cpatch%lai                   (recc)  &
                                                    + cpatch%qmean_fsw        (icyc,donc)  &
                                                    * cpatch%lai                   (donc)) &
                                                  * newlaii
            cpatch%qmean_fsn          (icyc,recc) = ( cpatch%qmean_fsn        (icyc,recc)  &
                                                    * cpatch%lai                   (recc)  &
                                                    + cpatch%qmean_fsn        (icyc,donc)  &
                                                    * cpatch%lai                   (donc)) &
                                                  * newlaii

            !----- The following variables are "extensive", add them. ---------------------!
            cpatch%qmean_par_l        (icyc,recc) = cpatch%qmean_par_l        (icyc,recc)  &
                                                  + cpatch%qmean_par_l        (icyc,donc)
            cpatch%qmean_par_l_beam   (icyc,recc) = cpatch%qmean_par_l_beam   (icyc,recc)  &
                                                  + cpatch%qmean_par_l_beam   (icyc,donc)
            cpatch%qmean_par_l_diff   (icyc,recc) = cpatch%qmean_par_l_diff   (icyc,recc)  &
                                                  + cpatch%qmean_par_l_diff   (icyc,donc)
            cpatch%qmean_psi_open     (icyc,recc) = cpatch%qmean_psi_open     (icyc,recc)  &
                                                  + cpatch%qmean_psi_open     (icyc,donc)
            cpatch%qmean_psi_closed   (icyc,recc) = cpatch%qmean_psi_closed   (icyc,recc)  &
                                                  + cpatch%qmean_psi_closed   (icyc,donc)
            cpatch%qmean_water_supply (icyc,recc) = cpatch%qmean_water_supply (icyc,recc)  &
                                                  + cpatch%qmean_water_supply (icyc,donc)
         end do
      end if



      !------------------------------------------------------------------------------------!
      !     Lastly, we update nplant and LAI.                                              !
      !------------------------------------------------------------------------------------!
      cpatch%nplant(recc) = newn
      !------------------------------------------------------------------------------------!
      !    LAI must be zero if phenology status is 2.  This is probably done correctly     !
      ! throughout the code, but being safe here.                                          !
      !------------------------------------------------------------------------------------!
      cpatch%lai(recc) = cpatch%lai(recc) + cpatch%lai(donc)
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
   subroutine fuse_patches(cgrid,ifm)
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
      real   , dimension(ff_nhgt)         :: laimax_cum      ! Mean # of plants
      integer                             :: ipy             ! Counters
      integer                             :: isi             ! Counters
      integer                             :: jpy             ! Counters
      integer                             :: jsi             ! Counters
      integer                             :: ipa             ! Counters
      integer                             :: ico             ! Counters
      integer                             :: donp            ! Counters
      integer                             :: recp            ! Counters
      integer                             :: ipft            ! Counters
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
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     First time here.  Delete all files.                                            !
      !------------------------------------------------------------------------------------!
      if (first_time .and. print_fuse_details) then
         do jpy = 1, cgrid%npolygons
            jpoly => cgrid%polygon(jpy)
            do jsi = 1, jpoly%nsites
               write (fuse_fout,fmt='(a,2(a,i4.4),a)')                                     &
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
         end do
         first_time = .false.
      end if
      !---------------------------------------------------------------------------------------!



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
               
               !----- If patch is not empty, or has already been fused, move on. ----------!
               if ( (.not. fuse_table(donp)) .or. donpatch%ncohorts > 0) cycle donloope


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

                  !------------------------------------------------------------------------!
                  !     Skip the patch if it isn't empty, or it has already been fused, or !
                  ! if the donor and receptor have different disturbance types.            !
                  !------------------------------------------------------------------------!
                  if ( (.not. fuse_table(recp))                       .or.                 &
                       recpatch%ncohorts > 0                          .or.                 &
                       csite%dist_type(donp) /= csite%dist_type(recp)     ) then
                     cycle recloope
                  end if
                  !------------------------------------------------------------------------!

                  !----- Skip the patch if they don't have the same disturbance type. -----!
                  if ( csite%dist_type(donp) /= csite%dist_type(recp)) cycle recloope
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !     Take an average of the patch properties of donpatch and recpatch,  !
                  ! and assign the average recpatch.                                       !
                  !------------------------------------------------------------------------!
                  call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)%prss          &
                                     ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)               &
                                     ,cpoly%green_leaf_factor(:,isi),elim_nplant,elim_lai)


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
                     
                     !----- If patch is not empty, or has already been fused, move on. ----!
                     if ( (.not. fuse_table(donp)) .or. donpatch%ncohorts == 0) then
                        cycle donloopa
                     end if


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

                        !------------------------------------------------------------------!
                        !     Skip the patch if it isn't empty, or it has already been     !
                        ! fused, or if the donor and receptor have different disturbance   !
                        ! types.                                                           !
                        !------------------------------------------------------------------!
                        if ( (.not. fuse_table(recp))                       .or.           &
                             recpatch%ncohorts == 0                         .or.           &
                             csite%dist_type(donp) /= csite%dist_type(recp) .or.           &
                             csite%age(donp)      /=  csite%age(recp)    ) then
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
                           if (.not. fuse_flag) cycle recloopa
                           !---------------------------------------------------------------!
                        end do hgtloopa
                        !------------------------------------------------------------------!
                        !------------------------------------------------------------------!
                        !     Take an average of the patch properties of donpatch and      !
                        ! recpatch, and assign the average recpatch.                       !
                        !------------------------------------------------------------------!
                        call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)%prss    &
                                           ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)         &
                                           ,cpoly%green_leaf_factor(:,isi),elim_nplant     &
                                           ,elim_lai)


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

                  !----- Increment tolerance ----------------------------------------------!
                  sunny_toler =     sunny_toler * sunny_cumlai_mult
                  dark_toler  = max(dark_toler  * dark_cumlai_mult , dark_lai80 )
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
                     recp_found = csite%dist_type(donp) == csite%dist_type(recp) .and.     &
                                  fuse_table(recp) .and.                                   &
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
                  if (donpatch%ncohorts == 0 .and. recpatch%ncohorts == 0) cycle donloopp
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
                     if (.not. fuse_flag) cycle donloopp
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
                                     ,cpoly%green_leaf_factor(:,isi),elim_nplant,elim_lai)
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
               if (npatches_new <= abs(maxpatch)) exit mainfuseloop
               !---------------------------------------------------------------------------!

               !----- Increment tolerance -------------------------------------------------!
               sunny_toler =     sunny_toler * sunny_cumlai_mult
               dark_toler  = max(dark_toler  * dark_cumlai_mult , dark_lai80 )
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
                               ,elim_nplant,elim_lai)
      use ed_state_vars      , only : sitetype              & ! Structure 
                                    , patchtype             ! ! Structure
      use soil_coms          , only : soil                  ! ! intent(in), lookup table
      use ed_max_dims        , only : n_pft                 & ! intent(in)
                                    , n_dbh                 ! ! intent(in)
      use mem_polygons       , only : maxcohort             ! ! intent(in)
      use therm_lib          , only : uextcm2tl             ! ! function
      use ed_misc_coms       , only : iqoutput              & ! intent(in)
                                    , idoutput              & ! intent(in)
                                    , imoutput              & ! intent(in)
                                    , ndcycle               ! ! intent(in)
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
      real                   , intent(out) :: elim_nplant       ! Eliminated nplant 
      real                   , intent(out) :: elim_lai          ! Eliminated lai
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer     :: cpatch            ! Current patch
      type(patchtype)        , pointer     :: temppatch         ! Temporary patch
      integer                              :: ico,iii,icyc      ! Counters
      integer                              :: ndc               ! # of cohorts - donp patch
      integer                              :: nrc               ! # of cohorts - recp patch
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

      !----- Assign eliminated LAI and nplant to zero (everything stays) ------------------!
      elim_nplant = 0.
      elim_lai    = 0.

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
                                         * csite%sfcwater_mass(iii,recp)
      end do
      do iii=1,csite%nlev_sfcwater(donp)
         csite%sfcwater_energy(iii,donp) = csite%sfcwater_energy(iii,donp)                 &
                                         * csite%sfcwater_mass(iii,donp)
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

      csite%mean_rh(recp)                = newareai *                                      &
                                         ( csite%mean_rh(donp)        * csite%area(donp)   &
                                         + csite%mean_rh(recp)        * csite%area(recp) )

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
      !    Even though these variables are not prognostic, they need to be copied so the   !
      ! output will have the values.  Other variables will probably be scaled here as      !
      ! well.                                                                              !
      !------------------------------------------------------------------------------------!
      csite%avg_ustar     (recp)      = newareai *                                         &
                                      ( csite%avg_ustar     (donp)    * csite%area(donp)   &
                                      + csite%avg_ustar     (recp)    * csite%area(recp) )

      csite%avg_tstar     (recp)      = newareai *                                         &
                                      ( csite%avg_tstar     (donp)    * csite%area(donp)   &
                                      + csite%avg_tstar     (recp)    * csite%area(recp) )

      csite%avg_qstar     (recp)      = newareai *                                         &
                                      ( csite%avg_qstar     (donp)    * csite%area(donp)   &
                                      + csite%avg_qstar     (recp)    * csite%area(recp) )

      csite%avg_cstar     (recp)      = newareai *                                         &
                                      ( csite%avg_cstar     (donp)    * csite%area(donp)   &
                                      + csite%avg_cstar     (recp)    * csite%area(recp) )

      csite%avg_rshort_gnd(recp)      = newareai *                                         &
                                      ( csite%avg_rshort_gnd(donp)    * csite%area(donp)   &
                                      + csite%avg_rshort_gnd(recp)    * csite%area(recp) )

      csite%avg_rlong_gnd(recp)       = newareai *                                         &
                                      ( csite%avg_rlong_gnd(donp)     * csite%area(donp)   &
                                      + csite%avg_rlong_gnd(recp)     * csite%area(recp) )

      csite%avg_carbon_ac(recp)       = newareai *                                         &
                                      ( csite%avg_carbon_ac(donp)     * csite%area(donp)   &
                                      + csite%avg_carbon_ac(recp)     * csite%area(recp) )

      csite%avg_carbon_st(recp)       = newareai *                                         &
                                      ( csite%avg_carbon_st(donp)     * csite%area(donp)   &
                                      + csite%avg_carbon_st(recp)     * csite%area(recp) )

      csite%avg_vapor_lc(recp)        = newareai *                                         &
                                      ( csite%avg_vapor_lc(donp)      * csite%area(donp)   &
                                      + csite%avg_vapor_lc(recp)      * csite%area(recp) )  

      csite%avg_vapor_wc(recp)        = newareai *                                         &
                                      ( csite%avg_vapor_wc(donp)      * csite%area(donp)   &
                                      + csite%avg_vapor_wc(recp)      * csite%area(recp) )  

      csite%avg_vapor_gc(recp)        = newareai *                                         &
                                      ( csite%avg_vapor_gc(donp)      * csite%area(donp)   &
                                      + csite%avg_vapor_gc(recp)      * csite%area(recp) )  

      csite%avg_wshed_vg(recp)        = newareai *                                         &
                                      ( csite%avg_wshed_vg(donp)      * csite%area(donp)   &
                                      + csite%avg_wshed_vg(recp)      * csite%area(recp) )  

      csite%avg_intercepted(recp)     = newareai *                                         &
                                      ( csite%avg_intercepted(donp)   * csite%area(donp)   &
                                      + csite%avg_intercepted(recp)   * csite%area(donp) )

      csite%avg_throughfall(recp)     = newareai *                                         &
                                      ( csite%avg_throughfall(donp)   * csite%area(donp)   &
                                      + csite%avg_throughfall(recp)   * csite%area(donp) )

      csite%avg_vapor_ac(recp)        = newareai *                                         &
                                      ( csite%avg_vapor_ac(donp)      * csite%area(donp)   &
                                      + csite%avg_vapor_ac(recp)      * csite%area(recp) )  

      csite%avg_transp(recp)          = newareai *                                         &
                                      ( csite%avg_transp(donp)        * csite%area(donp)   &
                                      + csite%avg_transp(recp)        * csite%area(recp) )  

      csite%avg_evap(recp)            = newareai *                                         &
                                      ( csite%avg_evap(donp)          * csite%area(donp)   &
                                      + csite%avg_evap(recp)          * csite%area(recp) )  

      csite%avg_runoff(recp)          = newareai *                                         &
                                      ( csite%avg_runoff(donp)        * csite%area(donp)   &
                                      + csite%avg_runoff(recp)        * csite%area(recp) )  

      csite%avg_drainage(recp)        = newareai *                                         &
                                      ( csite%avg_drainage(donp)      * csite%area(donp)   &
                                      + csite%avg_drainage(recp)      * csite%area(recp) )  

      csite%aux(recp)                 = newareai *                                         &
                                      ( csite%aux(donp)               * csite%area(donp)   &
                                      + csite%aux(recp)               * csite%area(recp) )  

      csite%avg_sensible_lc(recp)     = newareai *                                         &
                                      ( csite%avg_sensible_lc(donp)   * csite%area(donp)   &
                                      + csite%avg_sensible_lc(recp)   * csite%area(recp) )  

      csite%avg_sensible_wc(recp)     = newareai *                                         &
                                      ( csite%avg_sensible_wc(donp)   * csite%area(donp)   &
                                      + csite%avg_sensible_wc(recp)   * csite%area(recp) )  

      csite%avg_qwshed_vg(recp)       = newareai *                                         &
                                      ( csite%avg_qwshed_vg(donp)     * csite%area(donp)   &
                                      + csite%avg_qwshed_vg(recp)     * csite%area(recp) )  

      csite%avg_qintercepted(recp)    = newareai *                                         &
                                      ( csite%avg_qintercepted(donp)  * csite%area(donp)   &
                                      + csite%avg_qintercepted(recp)  * csite%area(donp) )

      csite%avg_qthroughfall(recp)    = newareai *                                         &
                                      ( csite%avg_qthroughfall(donp)  * csite%area(donp)   &
                                      + csite%avg_qthroughfall(recp)  * csite%area(donp) )

      csite%avg_sensible_gc(recp)     = newareai *                                         &
                                      ( csite%avg_sensible_gc(donp)   * csite%area(donp)   &
                                      + csite%avg_sensible_gc(recp)   * csite%area(recp) )  

      csite%avg_sensible_ac(recp)     = newareai *                                         &
                                      ( csite%avg_sensible_ac(donp)   * csite%area(donp)   &
                                      + csite%avg_sensible_ac(recp)   * csite%area(recp) )  

      csite%avg_runoff_heat(recp)     = newareai *                                         &
                                      ( csite%avg_runoff_heat(donp)   * csite%area(donp)   &
                                      + csite%avg_runoff_heat(recp)   * csite%area(recp) )  

      csite%avg_drainage_heat(recp)   = newareai *                                         &
                                      ( csite%avg_drainage_heat(donp) * csite%area(donp)   &
                                      + csite%avg_drainage_heat(recp) * csite%area(recp) )  

      csite%avg_leaf_energy(recp)     = newareai *                                         &
                                      ( csite%avg_leaf_energy(donp)   * csite%area(donp)   &
                                      + csite%avg_leaf_energy(recp)   * csite%area(recp) )

      csite%avg_leaf_water(recp)      = newareai *                                         &
                                      ( csite%avg_leaf_water(donp)    * csite%area(donp)   &
                                      + csite%avg_leaf_water(recp)    * csite%area(recp) )

      csite%avg_leaf_hcap(recp)       = newareai *                                         &
                                      ( csite%avg_leaf_hcap(donp)     * csite%area(donp)   &
                                      + csite%avg_leaf_hcap(recp)     * csite%area(recp) )

      csite%avg_wood_energy(recp)     = newareai *                                         &
                                      ( csite%avg_wood_energy(donp)   * csite%area(donp)   &
                                      + csite%avg_wood_energy(recp)   * csite%area(recp) )  

      csite%avg_wood_water(recp)      = newareai *                                         &
                                      ( csite%avg_wood_water(donp)    * csite%area(donp)   &
                                      + csite%avg_wood_water(recp)    * csite%area(recp) )  

      csite%avg_wood_hcap(recp)       = newareai *                                         &
                                      ( csite%avg_wood_hcap(donp)     * csite%area(donp)   &
                                      + csite%avg_wood_hcap(recp)     * csite%area(recp) )

      csite%co2budget_residual(recp)  = newareai *                                         &
                                      ( csite%co2budget_residual(donp)* csite%area(donp)   &
                                      + csite%co2budget_residual(recp)* csite%area(recp) )  

      csite%co2budget_loss2atm(recp)  = newareai *                                         &
                                      ( csite%co2budget_loss2atm(donp)* csite%area(donp)   &
                                      + csite%co2budget_loss2atm(recp)* csite%area(recp) )  

      csite%co2budget_denseffect(recp)= newareai *                                         &
                      ( csite%co2budget_denseffect(donp) * csite%area(donp)                &
                      + csite%co2budget_denseffect(recp) * csite%area(recp) )

      csite%co2budget_gpp(recp)       = newareai *                                         &
                                      ( csite%co2budget_gpp(donp)     * csite%area(donp)   &
                                      + csite%co2budget_gpp(recp)     * csite%area(recp) )  

      csite%co2budget_plresp(recp)    = newareai *                                         &
                                      ( csite%co2budget_plresp(donp)  * csite%area(donp)   &
                                      + csite%co2budget_plresp(recp)  * csite%area(recp) )  

      csite%co2budget_rh(recp)        = newareai *                                         &
                                      ( csite%co2budget_rh(donp)      * csite%area(donp)   &
                                      + csite%co2budget_rh(recp)      * csite%area(recp) )  

      csite%ebudget_residual(recp)    = newareai *                                         &
                                      ( csite%ebudget_residual(donp)  * csite%area(donp)   &
                                      + csite%ebudget_residual(recp)  * csite%area(recp) )

      csite%ebudget_netrad(recp)      = newareai *                                         &
                                      ( csite%ebudget_netrad  (donp)  * csite%area(donp)   &
                                      + csite%ebudget_netrad  (recp)  * csite%area(recp) )

      csite%ebudget_loss2atm(recp)    = newareai *                                         &
                                      ( csite%ebudget_loss2atm(donp)  * csite%area(donp)   &
                                      + csite%ebudget_loss2atm(recp)  * csite%area(recp) )

      csite%ebudget_denseffect(recp)  = newareai *                                         &
                                      ( csite%ebudget_denseffect(donp) * csite%area(donp)  &
                                      + csite%ebudget_denseffect(recp) * csite%area(recp) )

      csite%ebudget_prsseffect(recp)  = newareai *                                         &
                                      ( csite%ebudget_prsseffect(donp) * csite%area(donp)  &
                                      + csite%ebudget_prsseffect(recp) * csite%area(recp) )

      csite%ebudget_loss2runoff(recp) = newareai *                                         &
                                     ( csite%ebudget_loss2runoff(donp) * csite%area(donp)  &
                                     + csite%ebudget_loss2runoff(recp) * csite%area(recp) )

      csite%ebudget_loss2drainage(recp) = newareai *                                       &
                                   ( csite%ebudget_loss2drainage(donp) * csite%area(donp)  &
                                   + csite%ebudget_loss2drainage(recp) * csite%area(recp) )


      csite%ebudget_precipgain(recp)  = newareai *                                         &
                                  ( csite%ebudget_precipgain(donp) * csite%area(donp)      &
                                  + csite%ebudget_precipgain(recp) * csite%area(recp) )

      csite%wbudget_residual(recp)    = newareai *                                         &
                                      ( csite%wbudget_residual(donp)  * csite%area(donp)   &
                                      + csite%wbudget_residual(recp)  * csite%area(recp) )

      csite%wbudget_loss2atm(recp)    = newareai *                                         &
                                      ( csite%wbudget_loss2atm(donp)  * csite%area(donp)   &
                                      + csite%wbudget_loss2atm(recp)  * csite%area(recp) )

      csite%wbudget_denseffect(recp)  = newareai *                                         &
                                      ( csite%wbudget_denseffect(donp) * csite%area(donp)  &
                                      + csite%wbudget_denseffect(recp) * csite%area(recp) )

      csite%wbudget_loss2runoff(recp) = newareai *                                         &
                                     ( csite%wbudget_loss2runoff(donp) * csite%area(donp)  &
                                     + csite%wbudget_loss2runoff(recp) * csite%area(recp) )

      csite%wbudget_loss2drainage(recp) = newareai *                                       &
                                   ( csite%wbudget_loss2drainage(donp) * csite%area(donp)  &
                                   + csite%wbudget_loss2drainage(recp) * csite%area(recp) )

      csite%wbudget_precipgain(recp)  = newareai *                                         &
                                  ( csite%wbudget_precipgain(donp) * csite%area(donp)      &
                                  + csite%wbudget_precipgain(recp) * csite%area(recp) )


      do iii=1,mzg
         csite%avg_smoist_gg(iii,recp)   = newareai *                                      &
              ( csite%avg_smoist_gg(iii,donp)       * csite%area(donp)                     &
              + csite%avg_smoist_gg(iii,recp)       * csite%area(recp) )

         csite%avg_transloss(iii,recp)   = newareai *                                      &
              ( csite%avg_transloss(iii,donp)       * csite%area(donp)                     &
              + csite%avg_transloss(iii,recp)       * csite%area(recp) )

         csite%aux_s(iii,recp)           = newareai *                                      &
              ( csite%aux_s(iii,donp)               * csite%area(donp)                     &
              + csite%aux_s(iii,recp)               * csite%area(recp) )

         csite%avg_sensible_gg(iii,recp) = newareai *                                      &
              ( csite%avg_sensible_gg(iii,donp)     * csite%area(donp)                     &
              + csite%avg_sensible_gg(iii,recp)     * csite%area(recp) )
      end do

      do iii=1,n_dbh
         csite%co2budget_gpp_dbh(iii,recp) = newareai *                                    &
              ( csite%co2budget_gpp_dbh(iii,donp)   * csite%area(donp)                     &
              + csite%co2budget_gpp_dbh(iii,recp)   * csite%area(recp) )
      end do
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------! 
      !    We must also check whether daily/monthly output variables exist.  If they do,   !
      ! then we must fuse them too.                                                        !
      !------------------------------------------------------------------------------------! 
      if (idoutput > 0 .or. imoutput > 0 .or. iqoutput > 0) then
         csite%dmean_rh(recp)           = newareai                                         &
                                        * ( csite%dmean_rh(donp) * csite%area(donp)        &
                                          + csite%dmean_rh(recp) * csite%area(recp) )
         csite%dmean_co2_residual(recp) = newareai                                         &
                                        * ( csite%dmean_co2_residual(donp)                 &
                                          * csite%area(donp)                               &
                                          + csite%dmean_co2_residual(recp)                 &
                                          * csite%area(recp) )
         csite%dmean_energy_residual(recp) = newareai                                      &
                                           * ( csite%dmean_energy_residual(donp)           &
                                             * csite%area(donp)                            &
                                             + csite%dmean_energy_residual(recp)           &
                                             * csite%area(recp) )
         csite%dmean_water_residual(recp)  = newareai                                      &
                                           * ( csite%dmean_water_residual(donp)            &
                                             * csite%area(donp)                            &
                                             + csite%dmean_water_residual(recp)            &
                                             * csite%area(recp) )
         csite%dmean_lambda_light(recp)    = newareai                                      &
                                           * ( csite%dmean_lambda_light(donp)              &
                                             * csite%area(donp)                            &
                                             + csite%dmean_lambda_light(recp)              &
                                             * csite%area(recp) )
         csite%dmean_A_decomp(recp)        = newareai                                      &
                                           * ( csite%dmean_A_decomp(donp)                  &
                                             * csite%area(donp)                            &
                                             + csite%dmean_A_decomp(recp)                  &
                                             * csite%area(recp) )

         csite%dmean_Af_decomp(recp)       = newareai                                      &
                                           * ( csite%dmean_Af_decomp(donp)                 &
                                             * csite%area(donp)                            &
                                             + csite%dmean_Af_decomp(recp)                 &
                                             * csite%area(recp) )
         csite%dmean_rk4step(recp)         = newareai                                      &
                                           * ( csite%dmean_rk4step(donp)                   &
                                             * csite%area(donp)                            &
                                             + csite%dmean_rk4step(recp)                   &
                                             * csite%area(recp) )
      end if
      if (imoutput > 0 .or. iqoutput > 0) then
         csite%mmean_rh(recp)           = newareai                                         &
                                        * ( csite%mmean_rh(donp) * csite%area(donp)        &
                                          + csite%mmean_rh(recp) * csite%area(recp) )
         csite%mmean_co2_residual(recp) = newareai                                         &
                                        * ( csite%mmean_co2_residual(donp)                 &
                                          * csite%area(donp)                               &
                                          + csite%mmean_co2_residual(recp)                 &
                                          * csite%area(recp) )
         csite%mmean_energy_residual(recp) = newareai                                      &
                                           * ( csite%mmean_energy_residual(donp)           &
                                             * csite%area(donp)                            &
                                             + csite%mmean_energy_residual(recp)           &
                                             * csite%area(recp) )
         csite%mmean_water_residual(recp)  = newareai                                      &
                                           * ( csite%mmean_water_residual(donp)            &
                                             * csite%area(donp)                            &
                                             + csite%mmean_water_residual(recp)            &
                                             * csite%area(recp) )
         csite%mmean_lambda_light(recp)    = newareai                                      &
                                           * ( csite%mmean_lambda_light(donp)              &
                                             * csite%area(donp)                            &
                                             + csite%mmean_lambda_light(recp)              &
                                             * csite%area(recp) )
         csite%mmean_A_decomp(recp)        = newareai                                      &
                                           * ( csite%mmean_A_decomp(donp)                  &
                                             * csite%area(donp)                            &
                                             + csite%mmean_A_decomp(recp)                  &
                                             * csite%area(recp) )
         csite%mmean_Af_decomp(recp)       = newareai                                      &
                                           * ( csite%mmean_Af_decomp(donp)                 &
                                             * csite%area(donp)                            &
                                             + csite%mmean_Af_decomp(recp)                 &
                                             * csite%area(recp) )
         csite%mmean_rk4step(recp)         = newareai                                      &
                                           * ( csite%mmean_rk4step(donp)                   &
                                             * csite%area(donp)                            &
                                             + csite%mmean_rk4step(recp)                   &
                                             * csite%area(recp) )
      end if

      if (iqoutput > 0) then
         do icyc=1,ndcycle
            csite%qmean_rh     (icyc,recp) = newareai                                      &
                                           * ( csite%qmean_rh               (icyc,donp)    &
                                             * csite%area                        (donp)    &
                                             + csite%qmean_rh               (icyc,recp)    &
                                             * csite%area                        (recp))
         end do
      end if

      !------------------------------------------------------------------------------------!
      !      Update the leaf and wood temperature and liquid fraction.  We must check      !
      ! whether we can solve the average or not, because it may be an empty patch or the   !
      ! user may have disabled branchwood thermodynamics.                                  !
      !------------------------------------------------------------------------------------!
      if (csite%avg_leaf_hcap(recp) > 0.) then
         call uextcm2tl(csite%avg_leaf_energy(recp),csite%avg_leaf_water(recp)             &
                       ,csite%avg_leaf_hcap(recp),csite%avg_leaf_temp(recp)                &
                       ,csite%avg_leaf_fliq(recp))
      else
         csite%avg_leaf_temp(recp) = newareai                                              &
                                   * ( csite%avg_leaf_temp(donp) * csite%area(donp)        &
                                     + csite%avg_leaf_temp(recp) * csite%area(recp) )
         csite%avg_leaf_fliq(recp) = newareai                                              &
                                   * ( csite%avg_leaf_fliq(donp) * csite%area(donp)        &
                                     + csite%avg_leaf_fliq(recp) * csite%area(recp) )
      end if
      if (csite%avg_wood_hcap(recp) > 0.) then
         call uextcm2tl(csite%avg_wood_energy(recp),csite%avg_wood_water(recp)             &
                       ,csite%avg_wood_hcap(recp),csite%avg_wood_temp(recp)                &
                       ,csite%avg_wood_fliq(recp))
      else
         csite%avg_wood_temp(recp) = newareai                                              &
                                   * ( csite%avg_wood_temp(donp) * csite%area(donp)        &
                                     + csite%avg_wood_temp(recp) * csite%area(recp) )
         csite%avg_wood_fliq(recp) = newareai                                              &
                                   * ( csite%avg_wood_fliq(donp) * csite%area(donp)        &
                                     + csite%avg_wood_fliq(recp) * csite%area(recp) )
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
      !------------------------------------------------------------------------------------!
      ! IMPORTANT: Only cohort-level variables that have units per area (m2) should be     !
      !            rescaled.  Variables whose units are per plant should _NOT_ be included !
      !            here.                                                                   !
      !------------------------------------------------------------------------------------!
      do ico = 1,nrc
         cpatch%lai                   (ico) = cpatch%lai                (ico)  * area_scale
         cpatch%wpa                   (ico) = cpatch%wpa                (ico)  * area_scale
         cpatch%wai                   (ico) = cpatch%wai                (ico)  * area_scale
         cpatch%nplant                (ico) = cpatch%nplant             (ico)  * area_scale
         cpatch%mean_gpp              (ico) = cpatch%mean_gpp           (ico)  * area_scale
         cpatch%mean_leaf_resp        (ico) = cpatch%mean_leaf_resp     (ico)  * area_scale
         cpatch%mean_root_resp        (ico) = cpatch%mean_root_resp     (ico)  * area_scale
         cpatch%mean_growth_resp      (ico) = cpatch%mean_growth_resp   (ico)  * area_scale
         cpatch%mean_storage_resp     (ico) = cpatch%mean_storage_resp  (ico)  * area_scale
         cpatch%mean_vleaf_resp       (ico) = cpatch%mean_vleaf_resp    (ico)  * area_scale
         cpatch%today_gpp             (ico) = cpatch%today_gpp          (ico)  * area_scale
         cpatch%today_nppleaf         (ico) = cpatch%today_nppleaf      (ico)  * area_scale
         cpatch%today_nppfroot        (ico) = cpatch%today_nppfroot     (ico)  * area_scale
         cpatch%today_nppsapwood      (ico) = cpatch%today_nppsapwood   (ico)  * area_scale
         cpatch%today_nppcroot        (ico) = cpatch%today_nppcroot     (ico)  * area_scale
         cpatch%today_nppseeds        (ico) = cpatch%today_nppseeds     (ico)  * area_scale
         cpatch%today_nppwood         (ico) = cpatch%today_nppwood      (ico)  * area_scale
         cpatch%today_nppdaily        (ico) = cpatch%today_nppdaily     (ico)  * area_scale
         cpatch%today_gpp_pot         (ico) = cpatch%today_gpp_pot      (ico)  * area_scale
         cpatch%today_gpp_max         (ico) = cpatch%today_gpp_max      (ico)  * area_scale
         cpatch%today_leaf_resp       (ico) = cpatch%today_leaf_resp    (ico)  * area_scale
         cpatch%today_root_resp       (ico) = cpatch%today_root_resp    (ico)  * area_scale
         cpatch%Psi_open              (ico) = cpatch%Psi_open           (ico)  * area_scale
         cpatch%gpp                   (ico) = cpatch%gpp                (ico)  * area_scale
         cpatch%leaf_respiration      (ico) = cpatch%leaf_respiration   (ico)  * area_scale
         cpatch%root_respiration      (ico) = cpatch%root_respiration   (ico)  * area_scale
         cpatch%monthly_dndt          (ico) = cpatch%monthly_dndt       (ico)  * area_scale
         cpatch%leaf_water            (ico) = cpatch%leaf_water         (ico)  * area_scale
         cpatch%leaf_hcap             (ico) = cpatch%leaf_hcap          (ico)  * area_scale
         cpatch%leaf_energy           (ico) = cpatch%leaf_energy        (ico)  * area_scale
         cpatch%wood_water            (ico) = cpatch%wood_water         (ico)  * area_scale
         cpatch%wood_hcap             (ico) = cpatch%wood_hcap          (ico)  * area_scale
         cpatch%wood_energy           (ico) = cpatch%wood_energy        (ico)  * area_scale
         !----- Crown area shall not exceed one. ---------------------------------------!
         cpatch%crown_area            (ico) = min(1.,cpatch%crown_area  (ico)  * area_scale)
         if (idoutput > 0 .or. imoutput > 0 .or. iqoutput > 0) then
            cpatch%dmean_par_l        (ico) = cpatch%dmean_par_l        (ico)  * area_scale
            cpatch%dmean_par_l_beam   (ico) = cpatch%dmean_par_l_beam   (ico)  * area_scale
            cpatch%dmean_par_l_diff   (ico) = cpatch%dmean_par_l_diff   (ico)  * area_scale
         end if
         if (imoutput > 0 .or. iqoutput > 0) then
            cpatch%mmean_par_l        (ico) = cpatch%mmean_par_l        (ico)  * area_scale
            cpatch%mmean_par_l_beam   (ico) = cpatch%mmean_par_l_beam   (ico)  * area_scale
            cpatch%mmean_par_l_diff   (ico) = cpatch%mmean_par_l_diff   (ico)  * area_scale
         end if
         if (iqoutput > 0) then
            cpatch%qmean_par_l      (:,ico) = cpatch%qmean_par_l      (:,ico)  * area_scale
            cpatch%qmean_par_l_beam (:,ico) = cpatch%qmean_par_l_beam (:,ico)  * area_scale
            cpatch%qmean_par_l_diff (:,ico) = cpatch%qmean_par_l_diff (:,ico)  * area_scale
         end if
      end do
      !----- 2. Adjust densities of cohorts in donor patch --------------------------------!
      cpatch => csite%patch(donp)
      ndc = cpatch%ncohorts
      area_scale = csite%area(donp) * newareai
      !------------------------------------------------------------------------------------!
      ! IMPORTANT: Only cohort-level variables that have units per area (m2) should be     !
      !            rescaled.  Variables whose units are per plant should _NOT_ be included !
      !            here.                                                                   !
      !------------------------------------------------------------------------------------!
      do ico = 1,ndc
         cpatch%lai                   (ico) = cpatch%lai                (ico)  * area_scale
         cpatch%wpa                   (ico) = cpatch%wpa                (ico)  * area_scale
         cpatch%wai                   (ico) = cpatch%wai                (ico)  * area_scale
         cpatch%nplant                (ico) = cpatch%nplant             (ico)  * area_scale
         cpatch%mean_gpp              (ico) = cpatch%mean_gpp           (ico)  * area_scale
         cpatch%mean_leaf_resp        (ico) = cpatch%mean_leaf_resp     (ico)  * area_scale
         cpatch%mean_root_resp        (ico) = cpatch%mean_root_resp     (ico)  * area_scale
         cpatch%mean_growth_resp      (ico) = cpatch%mean_growth_resp   (ico)  * area_scale
         cpatch%mean_storage_resp     (ico) = cpatch%mean_storage_resp  (ico)  * area_scale
         cpatch%mean_vleaf_resp       (ico) = cpatch%mean_vleaf_resp    (ico)  * area_scale
         cpatch%today_gpp             (ico) = cpatch%today_gpp          (ico)  * area_scale
         cpatch%today_nppleaf         (ico) = cpatch%today_nppleaf      (ico)  * area_scale
         cpatch%today_nppfroot        (ico) = cpatch%today_nppfroot     (ico)  * area_scale
         cpatch%today_nppsapwood      (ico) = cpatch%today_nppsapwood   (ico)  * area_scale
         cpatch%today_nppcroot        (ico) = cpatch%today_nppcroot     (ico)  * area_scale
         cpatch%today_nppseeds        (ico) = cpatch%today_nppseeds     (ico)  * area_scale
         cpatch%today_nppwood         (ico) = cpatch%today_nppwood      (ico)  * area_scale
         cpatch%today_nppdaily        (ico) = cpatch%today_nppdaily     (ico)  * area_scale
         cpatch%today_gpp_pot         (ico) = cpatch%today_gpp_pot      (ico)  * area_scale
         cpatch%today_gpp_max         (ico) = cpatch%today_gpp_max      (ico)  * area_scale
         cpatch%today_leaf_resp       (ico) = cpatch%today_leaf_resp    (ico)  * area_scale
         cpatch%today_root_resp       (ico) = cpatch%today_root_resp    (ico)  * area_scale
         cpatch%Psi_open              (ico) = cpatch%Psi_open           (ico)  * area_scale
         cpatch%gpp                   (ico) = cpatch%gpp                (ico)  * area_scale
         cpatch%leaf_respiration      (ico) = cpatch%leaf_respiration   (ico)  * area_scale
         cpatch%root_respiration      (ico) = cpatch%root_respiration   (ico)  * area_scale
         cpatch%monthly_dndt          (ico) = cpatch%monthly_dndt       (ico)  * area_scale
         cpatch%leaf_water            (ico) = cpatch%leaf_water         (ico)  * area_scale
         cpatch%leaf_hcap             (ico) = cpatch%leaf_hcap          (ico)  * area_scale
         cpatch%leaf_energy           (ico) = cpatch%leaf_energy        (ico)  * area_scale
         cpatch%wood_water            (ico) = cpatch%wood_water         (ico)  * area_scale
         cpatch%wood_hcap             (ico) = cpatch%wood_hcap          (ico)  * area_scale
         cpatch%wood_energy           (ico) = cpatch%wood_energy        (ico)  * area_scale
         !----- Crown area shall not exceed one. ---------------------------------------!
         cpatch%crown_area            (ico) = min(1.,cpatch%crown_area  (ico)  * area_scale)
         if (idoutput > 0 .or. imoutput > 0 .or. iqoutput > 0) then
            cpatch%dmean_par_l        (ico) = cpatch%dmean_par_l        (ico)  * area_scale
            cpatch%dmean_par_l_beam   (ico) = cpatch%dmean_par_l_beam   (ico)  * area_scale
            cpatch%dmean_par_l_diff   (ico) = cpatch%dmean_par_l_diff   (ico)  * area_scale
         end if
         if (imoutput > 0 .or. iqoutput > 0) then
            cpatch%mmean_par_l        (ico) = cpatch%mmean_par_l        (ico)  * area_scale
            cpatch%mmean_par_l_beam   (ico) = cpatch%mmean_par_l_beam   (ico)  * area_scale
            cpatch%mmean_par_l_diff   (ico) = cpatch%mmean_par_l_diff   (ico)  * area_scale
         end if
         if (iqoutput > 0) then
            cpatch%qmean_par_l      (:,ico) = cpatch%qmean_par_l      (:,ico)  * area_scale
            cpatch%qmean_par_l_beam (:,ico) = cpatch%qmean_par_l_beam (:,ico)  * area_scale
            cpatch%qmean_par_l_diff (:,ico) = cpatch%qmean_par_l_diff (:,ico)  * area_scale
         end if
      end do
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
            call fuse_cohorts(csite,recp,green_leaf_factor,lsl)
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
      ! + csite%lai(recp)                                                                  !
      ! + csite%veg_rough(recp)                                                            !
      ! + csite%total_sfcw_depth(recp)                                                     !
      ! + csite%snowfac(recp)                                                              !
      ! + csite%opencan_frac(recp)                                                         !
      !------------------------------------------------------------------------------------!
      call update_patch_derived_props(csite,lsl, prss,recp)
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
      use allometry           , only : dbh2bl     ! ! intent(in)
      use ed_max_dims         , only : n_pft      ! ! intent(in)
      use pft_coms            , only : hgt_min    ! ! intent(in)
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
         ihgt    = min(ff_nhgt,1 + count(hgt_class < cpatch%hite(ico)))
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Check whether this cohort is almost at the minimum height given its PFT.    !
         ! If it is, then we will skip it.                                                 !
         !---------------------------------------------------------------------------------!
         if (cpatch%hite(ico) < hgt_min(ipft) + 0.2) cycle cohortloop
         !---------------------------------------------------------------------------------!


         !----- Find the height class. ----------------------------------------------------!
         ihgt    = min(ff_nhgt,1 + count(hgt_class < cpatch%hite(ico)))
         !---------------------------------------------------------------------------------!


         !----- Find the potential (on-allometry) leaf area index. ------------------------!
         lai_pot = cpatch%nplant(ico) * cpatch%sla(ico)                                    &
                 * dbh2bl(cpatch%dbh(ico),ipft)
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
end module fuse_fiss_utils
!==========================================================================================!
!==========================================================================================!
