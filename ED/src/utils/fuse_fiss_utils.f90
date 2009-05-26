module fuse_fiss_utils_ar

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
   !   This subroutine will sort the cohorts by size (1st = tallest, last = shortest.)     !
   !---------------------------------------------------------------------------------------!
   subroutine sort_cohorts_ar(cpatch)

      use ed_state_vars,only :  patchtype   ! ! Structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target  :: cpatch     ! Current patch, to have cohorts sorted.
      !----- Local variables --------------------------------------------------------------!
      type(patchtype), pointer :: temppatch  ! Scratch patch structure
      integer                  :: ico, iico  ! Counters
      integer                  :: tallid     ! Identity of tallest cohort
      !------------------------------------------------------------------------------------!
      
      !----- No need to sort an empty patch or a patch with a single cohort. --------------!
      if (cpatch%ncohorts < 2) return

      !----- Assigning a scratch patch ----------------------------------------------------!
      nullify(temppatch)
      allocate(temppatch)
      call allocate_patchtype(temppatch,cpatch%ncohorts)
      
      iico = 1
      !---- Loop until all cohorts were sorted --------------------------------------------!
      do while(iico <= cpatch%ncohorts)
      
         !----- Finding the tallest cohort ------------------------------------------------!
         tallid = maxloc(cpatch%hite,dim=1)
         
         !----- Copying to the scratch structure ------------------------------------------!
         call copy_patchtype(cpatch,temppatch,tallid,tallid,iico,iico)
         
         !----- Putting a non-sense height so this will never "win" again. ----------------!
         cpatch%hite(tallid) = -huge(1.)

         iico = iico + 1
      end do

      !------ Copying the scratch patch to the regular one and deallocating it ------------!
      call copy_patchtype(temppatch,cpatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
      call deallocate_patchtype(temppatch)
      deallocate(temppatch)

      return

   end subroutine sort_cohorts_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will eliminate cohorts based on their sizes. This is intended to   !
   ! eliminate cohorts that have little contribution and thus we can speed up the run.     !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_cohorts_ar(csite,ipa,elim_nplant,elim_lai)
      use pft_coms           , only : min_recruit_size & ! intent(in)
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

         !----- Checking whether the cohort size is too small -----------------------------!
         csize = cpatch%nplant(ico)                                                        &
               * (cpatch%balive(ico) + cpatch%bdead(ico) + cpatch%bstorage(ico))

         if( csize < (0.1 * min_recruit_size(cpatch%pft(ico)))) then
            !----- Cohort is indeed too small, it won't remain ----------------------------!
            remain_table(ico) = .false.
            elim_nplant = elim_nplant + cpatch%nplant(ico)
            elim_lai    = elim_lai    + cpatch%lai(ico)

            !----- Update litter pools ----------------------------------------------------!
            csite%fsc_in(ipa) = csite%fsc_in(ipa) + cpatch%nplant(ico)                     &
                              * (f_labile(cpatch%pft(ico)) * cpatch%balive(ico)            &
                                + cpatch%bstorage(ico))

            csite%fsn_in(ipa) = csite%fsn_in(ipa) + cpatch%nplant(ico)                     &
                              * (f_labile(cpatch%pft(ico)) * cpatch%balive(ico)            &
                                 / c2n_leaf(cpatch%pft(ico))                               &
                                 + cpatch%bstorage(ico) / c2n_storage)
            
            csite%ssc_in(ipa) = csite%ssc_in(ipa) + cpatch%nplant(ico)                     &
                              * ((1.0 - f_labile(cpatch%pft(ico))) * cpatch%balive(ico)    &
                                 + cpatch%bdead(ico))
            
            csite%ssl_in(ipa) = csite%ssl_in(ipa) + cpatch%nplant(ico)                     &
                              * ((1.0 - f_labile(cpatch%pft(ico)))                         &
                                 *cpatch%balive(ico) + cpatch%bdead(ico))                  &
                              * l2n_stem/c2n_stem

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
   end subroutine terminate_cohorts_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine will eliminate tiny or empty patches. This is intended to          !
   ! eliminate patches that have little contribution and thus we can speed up the run.     !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_patches_ar(csite)

      use ed_state_vars, only :  polygontype       & ! Structure
                               , sitetype          & ! Structure
                               , patchtype         ! ! Structure
      use disturb_coms,  only : min_new_patch_area ! ! intent(in)

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
      !------------------------------------------------------------------------------------!
      new_area=0.
      area_scale = 1./(1. - elim_area)
      do ipa = 1,csite%npatches
         csite%area(ipa) = csite%area(ipa) * area_scale
         new_area = new_area + csite%area(ipa)

         cpatch => csite%patch(ipa)
         do ico = 1, cpatch%ncohorts
            cpatch%nplant(ico)              = cpatch%nplant(ico)              * area_scale
            cpatch%lai(ico)                 = cpatch%lai(ico)                 * area_scale
            cpatch%bai(ico)                 = cpatch%bai(ico)                 * area_scale
            cpatch%sai(ico)                 = cpatch%sai(ico)                 * area_scale
            cpatch%mean_gpp(ico)            = cpatch%mean_gpp(ico)            * area_scale
            cpatch%mean_leaf_resp(ico)      = cpatch%mean_leaf_resp(ico)      * area_scale
            cpatch%mean_root_resp(ico)      = cpatch%mean_root_resp(ico)      * area_scale
            cpatch%Psi_open(ico)            = cpatch%Psi_open(ico)            * area_scale
            cpatch%gpp(ico)                 = cpatch%gpp(ico)                 * area_scale
            cpatch%leaf_respiration(ico)    = cpatch%leaf_respiration(ico)    * area_scale
            cpatch%root_respiration(ico)    = cpatch%root_respiration(ico)    * area_scale
            cpatch%veg_water(ico)           = cpatch%veg_water(ico)           * area_scale
            cpatch%hcapveg(ico)             = cpatch%hcapveg(ico)             * area_scale
            cpatch%veg_energy(ico)          = cpatch%veg_energy(ico)          * area_scale
            cpatch%monthly_dndt(ico)        = cpatch%monthly_dndt(ico)        * area_scale
         end do
      end do

      if (abs(new_area-1.0) > 1.e-5) then
         write (unit=*,fmt='(a,1x,es12.5)') ' + ELIM_AREA:',elim_area
         write (unit=*,fmt='(a,1x,es12.5)') ' + NEW_AREA: ',new_area
         call fatal_error('New_area should be 1 but it isn''t!!!','terminate_patches_ar'   &
                       &,'fuse_fiss_utils.f90')
      end if 
      
      return
   end subroutine terminate_patches_ar
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
   subroutine fuse_cohorts_ar(csite,ipa, green_leaf_factor, lsl)

      use ed_state_vars       , only : sitetype            & ! Structure
                                     , patchtype           ! ! Structure
      use pft_coms            , only : rho                 & ! intent(in)
                                     , b1Ht                & ! intent(in)
                                     , max_dbh               ! intent(in)
      use fusion_fission_coms , only : fusetol_h           & ! intent(in)
                                     , fusetol             & ! intent(in)
                                     , lai_fuse_tol        & ! intent(in)
                                     , fuse_relax          & ! intent(in)
                                     , coh_tolerance_max   ! ! intent(in)
      use max_dims            , only : n_pft               ! ! intent(in)
      use mem_sites           , only : maxcohort           ! ! intent(in)
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
      real                                 :: hite_threshold ! Height threshold
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

      !----- Return if maxcohort is 0 (flag for no cohort fusion). ------------------------!
      if (maxcohort == 0) return

      !----- Start with no factor ---------------------------------------------------------!
      tolerance_mult = 1.0

      cpatch => csite%patch(ipa)
     
      !------------------------------------------------------------------------------------!
      !     Return if we don't have many cohorts anyway. Cohort fusion is just for         !
      ! computational efficiency. The abs is necessary because maxcohort can be negative.  !
      ! Negative is a flag to tell ED that cohort fusion is allowed to take place only     !
      ! during initialisation.                                                             !
      !------------------------------------------------------------------------------------!
      if(cpatch%ncohorts <= abs(maxcohort))return


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
         !----- Get fusion height threshold -----------------------------------------------!
         if (rho(cpatch%pft(ico3)) == 0.0)then
            hite_threshold = b1Ht(cpatch%pft(ico3))
         else
            hite_threshold = dbh2h(cpatch%pft(ico3),max_dbh(cpatch%pft(ico3)))
         end if

         if (cpatch%hite(ico3) < (0.95 * hite_threshold) ) then
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

               !----- Get fusion height threshold. ----------------------------------------!
               if(rho(cpatch%pft(donc)) == 0.0) then
                  hite_threshold = b1Ht(cpatch%pft(donc))
               else
                  hite_threshold = dbh2h(cpatch%pft(donc), max_dbh(cpatch%pft(donc)))
               end if

               !----- Test for similarity -------------------------------------------------!
               if (cpatch%hite(donc) >= (0.95 * hite_threshold )) then
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
                  lai_max = (cpatch%nplant(recc)*dbh2bl(cpatch%dbh(recc),cpatch%pft(recc)) &
                          + cpatch%nplant(donc)*dbh2bl(cpatch%dbh(donc),cpatch%pft(donc))) &
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
                     call fuse_2_cohorts_ar(cpatch,donc,recc,newn                          &
                                           ,green_leaf_factor(cpatch%pft(donc)),lsl)

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
                                        &,'fuse_2_cohorts_ar','fuse_fiss_utils.f90')
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
                        if (rho(cpatch%pft(ico3)) == 0.0) then
                           hite_threshold = b1Ht(cpatch%pft(ico3))
                        else
                           hite_threshold = dbh2h(cpatch%pft(ico3)                         &
                                                 , max_dbh(cpatch%pft(ico3)))
                        end if
                        if (cpatch%hite(ico3) < (0.95 * hite_threshold )) then
                           mean_hite = mean_hite + cpatch%hite(ico3)
                           nshort = nshort+1
                        else
                           if(cpatch%dbh(ico3).eq.0. ) then
                              print*,"dbh(ico3) is zero",cpatch%dbh(ico3)
                              call fatal_error('Zero DBH!','fuse_cohorts_ar'&
                                              &,'fuse_fiss_utils.f90')
                           end if
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
         call sort_cohorts_ar(cpatch)
         csite%cohort_count(ipa) = count(fuse_table)
      end if

      !----- Deallocate the aux. table ----------------------------------------------------!
      deallocate(fuse_table)
     
      return
   end subroutine fuse_cohorts_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will split two cohorts if its LAI has become too large.             !
   !---------------------------------------------------------------------------------------!
   subroutine split_cohorts_ar(cpatch, green_leaf_factor, lsl)

      use ed_state_vars        , only :  patchtype              ! ! structure
      use pft_coms             , only :  q                      & ! intent(in), lookup table
                                       , qsw                      ! intent(in), lookup table
      use fusion_fission_coms  , only :  lai_tol                ! ! intent(in)
      use max_dims             , only :  n_pft                  ! ! intent(in)
      use allometry            , only :  dbh2h                  & ! function
                                       , bd2dbh                 & ! function
                                       , dbh2bd                 ! ! function
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
      real                                 :: slai              ! LAI
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

         slai = cpatch%nplant(ico)                                                         &
              * cpatch%balive(ico)                                                         &
              * green_leaf_factor(cpatch%pft(ico))                                         &
              / ( 1.0 + q(cpatch%pft(ico)) + qsw(cpatch%pft(ico)) * cpatch%hite(ico) )     &
              * cpatch%sla(ico)

         !----- If the resulting LAI is too large, split this cohort. ---------------------!
         split_mask(ico) = slai > lai_tol
         
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
               !---------------------------------------------------------------------------!
               cpatch%lai(ico)                 = cpatch%lai(ico)                  * 0.5
               cpatch%bai(ico)                 = cpatch%bai(ico)                  * 0.5
               cpatch%sai(ico)                 = cpatch%sai(ico)                  * 0.5
               cpatch%nplant(ico)              = cpatch%nplant(ico)               * 0.5
               cpatch%mean_gpp(ico)            = cpatch%mean_gpp(ico)             * 0.5
               cpatch%mean_leaf_resp(ico)      = cpatch%mean_leaf_resp(ico)       * 0.5
               cpatch%mean_root_resp(ico)      = cpatch%mean_root_resp(ico)       * 0.5
               cpatch%dmean_gpp(ico)           = cpatch%dmean_gpp(ico)            * 0.5
               cpatch%dmean_gpp_pot(ico)       = cpatch%dmean_gpp_pot(ico)        * 0.5
               cpatch%dmean_gpp_max(ico)       = cpatch%dmean_gpp_max(ico)        * 0.5
               cpatch%dmean_leaf_resp(ico)     = cpatch%dmean_leaf_resp(ico)      * 0.5
               cpatch%dmean_root_resp(ico)     = cpatch%dmean_root_resp(ico)      * 0.5
               cpatch%Psi_open(ico)            = cpatch%Psi_open(ico)             * 0.5
               cpatch%gpp(ico)                 = cpatch%gpp(ico)                  * 0.5
               cpatch%leaf_respiration(ico)    = cpatch%leaf_respiration(ico)     * 0.5
               cpatch%root_respiration(ico)    = cpatch%root_respiration(ico)     * 0.5
               cpatch%monthly_dndt(ico)        = cpatch%monthly_dndt(ico)         * 0.5
               cpatch%veg_water(ico)           = cpatch%veg_water(ico)            * 0.5
               cpatch%hcapveg(ico)             = cpatch%hcapveg(ico)              * 0.5
               cpatch%veg_energy(ico)          = cpatch%veg_energy(ico)           * 0.5
               !---------------------------------------------------------------------------!


               !----- Apply these values to the new cohort. -------------------------------!
               inew = inew+1
               call clone_cohort_ar(cpatch,ico,inew)
               !---------------------------------------------------------------------------!

               !----- Tweaking bdead, to ensure carbon is conserved. ----------------------!
               cpatch%bdead(ico)  = cpatch%bdead(ico)*(1.-epsilon)
               cpatch%dbh(ico)    = bd2dbh(cpatch%pft(ico), cpatch%bdead(ico))
               cpatch%hite(ico)   = dbh2h(cpatch%pft(ico), cpatch%dbh(ico))

               cpatch%bdead(inew)  = cpatch%bdead(inew)*(1.+epsilon)
               cpatch%dbh(inew)    = bd2dbh(cpatch%pft(inew), cpatch%bdead(inew))
               cpatch%hite(inew)   = dbh2h(cpatch%pft(inew), cpatch%dbh(inew))
               !---------------------------------------------------------------------------!

            end if
         end do

         !----- After splitting, cohorts may need to be sorted again... -------------------!
         call sort_cohorts_ar(cpatch)

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
                                        &,'split_cohorts_ar','fuse_fiss_utils.f90')
         end if
         
      end if
      deallocate(split_mask)
      return
   end subroutine split_cohorts_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will clone one cohort.                                              !
   !---------------------------------------------------------------------------------------!
   subroutine clone_cohort_ar(cpatch,isc,idt)
   
      use ed_state_vars, only : patchtype  & ! Strucuture
                              , stoma_data ! ! Structure
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype) , target     :: cpatch ! Current patch
      integer         , intent(in) :: isc    ! Index of "Source" cohort
      integer         , intent(in) :: idt    ! Index of "Destination" cohort"
      !----- Local variables --------------------------------------------------------------!
      integer                      :: imonth
      type(stoma_data), pointer    :: osdt,ossc
      !------------------------------------------------------------------------------------!

      cpatch%pft(idt)                 = cpatch%pft(isc)
      cpatch%nplant(idt)              = cpatch%nplant(isc)
      cpatch%hite(idt)                = cpatch%hite(isc)
      cpatch%dbh(idt)                 = cpatch%dbh(isc)
      cpatch%bdead(idt)               = cpatch%bdead(isc)
      cpatch%bleaf(idt)               = cpatch%bleaf(isc)
      cpatch%phenology_status(idt)    = cpatch%phenology_status(isc)
      cpatch%balive(idt)              = cpatch%balive(isc)
      cpatch%lai(idt)                 = cpatch%lai(isc)
      cpatch%bai(idt)                 = cpatch%bai(isc)
      cpatch%sai(idt)                 = cpatch%sai(isc)
      cpatch%bstorage(idt)            = cpatch%bstorage(isc)
      
      !----- We copy this one, but it will be updated before it is really used. -----------!
      cpatch%solvable(idt)            = cpatch%solvable(isc)

      do imonth = 1,13
         cpatch%cb(imonth,idt)        = cpatch%cb(imonth,isc)
         cpatch%cb_max(imonth,idt)    = cpatch%cb_max(imonth,isc)
      enddo

      cpatch%cbr_bar(idt)             = cpatch%cbr_bar(isc)
      cpatch%veg_energy(idt)          = cpatch%veg_energy(isc)
      cpatch%veg_temp(idt)            = cpatch%veg_temp(isc)
      cpatch%veg_fliq(idt)            = cpatch%veg_fliq(isc)
      cpatch%veg_water(idt)           = cpatch%veg_water(isc)
      cpatch%mean_gpp(idt)            = cpatch%mean_gpp(isc)
      cpatch%mean_leaf_resp(idt)      = cpatch%mean_leaf_resp(isc)
      cpatch%mean_root_resp(idt)      = cpatch%mean_root_resp(isc)
      cpatch%dmean_leaf_resp(idt)     = cpatch%dmean_leaf_resp(isc)
      cpatch%dmean_root_resp(idt)     = cpatch%dmean_root_resp(isc)
      cpatch%dmean_gpp(idt)           = cpatch%dmean_gpp(isc)
      cpatch%dmean_gpp_pot(idt)       = cpatch%dmean_gpp_pot(isc)
      cpatch%dmean_gpp_max(idt)       = cpatch%dmean_gpp_max(isc)
      cpatch%growth_respiration(idt)  = cpatch%growth_respiration(isc)
      cpatch%storage_respiration(idt) = cpatch%storage_respiration(isc)
      cpatch%vleaf_respiration(idt)   = cpatch%vleaf_respiration(isc)
      cpatch%fsn(idt)                 = cpatch%fsn(isc)
      cpatch%monthly_dndt(idt)        = cpatch%monthly_dndt(isc)
      cpatch%Psi_open(idt)            = cpatch%Psi_open(isc)
      cpatch%krdepth(idt)             = cpatch%krdepth(isc)
      cpatch%first_census(idt)        = cpatch%first_census(isc)
      cpatch%new_recruit_flag(idt)    = cpatch%new_recruit_flag(isc)
      cpatch%par_v(idt)               = cpatch%par_v(isc)
      cpatch%par_v_beam(idt)          = cpatch%par_v_beam(isc)
      cpatch%par_v_diffuse(idt)       = cpatch%par_v_diffuse(isc)
      cpatch%rshort_v(idt)            = cpatch%rshort_v(isc)
      cpatch%rshort_v_beam(idt)       = cpatch%rshort_v_beam(isc)
      cpatch%rshort_v_diffuse(idt)    = cpatch%rshort_v_diffuse(isc)
      cpatch%rlong_v(idt)             = cpatch%rlong_v(isc)
      cpatch%rlong_v_surf(idt)        = cpatch%rlong_v_surf(isc)
      cpatch%rlong_v_incid(idt)       = cpatch%rlong_v_incid(isc)
      cpatch%rb(idt)                  = cpatch%rb(isc)
      cpatch%A_open(idt)              = cpatch%A_open(isc)
      cpatch%A_closed(idt)            = cpatch%A_closed(isc)
      cpatch%Psi_closed(idt)          = cpatch%Psi_closed(isc)
      cpatch%rsw_open(idt)            = cpatch%rsw_open(isc)
      cpatch%rsw_closed(idt)          = cpatch%rsw_closed(isc)
      cpatch%fsw(idt)                 = cpatch%fsw(isc)
      cpatch%fs_open(idt)             = cpatch%fs_open(isc)
      cpatch%stomatal_resistance(idt) = cpatch%stomatal_resistance(isc)
      cpatch%maintenance_costs(idt)   = cpatch%maintenance_costs(isc)
      cpatch%bseeds(idt)              = cpatch%bseeds(isc)
      cpatch%leaf_respiration(idt)    = cpatch%leaf_respiration(isc)
      cpatch%root_respiration(idt)    = cpatch%root_respiration(isc)
      cpatch%hcapveg(idt)             = cpatch%hcapveg(isc)

      cpatch%gpp(idt)                 = cpatch%gpp(isc)
      cpatch%paw_avg(idt)          = cpatch%paw_avg(isc)

      cpatch%turnover_amp(idt)  = cpatch%turnover_amp(isc)     
      cpatch%llspan(idt)  = cpatch%llspan(isc)     
      cpatch%vm_bar(idt)  = cpatch%vm_bar(isc)  
      cpatch%sla(idt)  = cpatch%sla(isc)  

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
     
      return
   end subroutine clone_cohort_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will merge two cohorts into 1. The donating cohort (donc) is the  !
   ! one that will be deallocated, while the receptor cohort (recc) will contain the       !
   !  information from both cohorts.                                                       !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_cohorts_ar(cpatch,donc,recc, newn,green_leaf_factor, lsl)
      use ed_state_vars , only :  patchtype             ! ! Structure
      use pft_coms      , only :  q                     & ! intent(in), lookup table
                                , qsw                     ! intent(in), lookup table
      use therm_lib     , only :  qwtk                  ! ! subroutine
      use allometry     , only :  calc_root_depth       & ! function
                                , assign_root_depth     & ! function
                                , bd2dbh                & ! function
                                , dbh2h                 ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype) , target     :: cpatch            ! Current patch
      integer                      :: donc              ! Donating cohort.
      integer                      :: recc              ! Receptor cohort.
      real            , intent(in) :: newn              ! New nplant
      real            , intent(in) :: green_leaf_factor ! Green leaf factor
      integer         , intent(in) :: lsl               ! Lowest soil level
      !----- Local variables --------------------------------------------------------------!
      integer                      :: imon              ! Month for cb loop
      real                         :: newni             ! Inverse of new nplants
      real                         :: cb_act            !
      real                         :: cb_max            !
      real                         :: root_depth        !
      !------------------------------------------------------------------------------------!

      newni = 1.0 / newn
     

      !----- Conserve carbon by calculating bdead first. ----------------------------------!
      cpatch%bdead(recc) = ( cpatch%nplant(recc) * cpatch%bdead(recc)                      &
                           + cpatch%nplant(donc) * cpatch%bdead(donc) ) * newni

      !----- Then get dbh and hite from bdead. --------------------------------------------!
      cpatch%dbh(recc)   = bd2dbh(cpatch%pft(recc), cpatch%bdead(recc))
      cpatch%hite(recc)  = dbh2h(cpatch%pft(recc),  cpatch%dbh(recc))


      !------------------------------------------------------------------------------------!
      !     Conserving carbon to get balive, bleaf, and bstorage.                          !
      !------------------------------------------------------------------------------------!
      cpatch%balive(recc) = ( cpatch%nplant(recc) * cpatch%balive(recc)                    &
                            + cpatch%nplant(donc) * cpatch%balive(donc) ) *newni
      cpatch%bstorage(recc) = ( cpatch%nplant(recc) * cpatch%bstorage(recc)                &
                              + cpatch%nplant(donc) * cpatch%bstorage(donc) ) * newni
      cpatch%bseeds(recc)   = ( cpatch%nplant(recc) * cpatch%bseeds(recc)                  &
                              + cpatch%nplant(donc) * cpatch%bseeds(donc) ) * newni
      cpatch%maintenance_costs(recc) = newni                                               &
                            * ( cpatch%nplant(recc) * cpatch%maintenance_costs(recc)       &
                              + cpatch%nplant(donc) * cpatch%maintenance_costs(donc) )
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Bleaf and LAI must be zero if phenology status is 2.  This is probably done     !
      ! correctly throughout the code, but being safe here.                                !
      !------------------------------------------------------------------------------------!
      if (cpatch%phenology_status(recc) < 2) then
         cpatch%bleaf(recc)  = ( cpatch%nplant(recc) * cpatch%bleaf(recc)                  &
                               + cpatch%nplant(donc) * cpatch%bleaf(donc) ) *newni
         cpatch%lai(recc)    = cpatch%lai(recc) + cpatch%lai(donc)
      else
         cpatch%bleaf(recc)      = 0.
         cpatch%lai(recc)        = 0.
      end if

      cpatch%bai(recc)    = cpatch%bai(recc)  + cpatch%bai(donc)
      cpatch%sai(recc)    = cpatch%sai(recc)  + cpatch%sai(donc)
      cpatch%veg_energy(recc) = cpatch%veg_energy(recc) + cpatch%veg_energy(donc)
      cpatch%veg_water(recc)  = cpatch%veg_water(recc)  + cpatch%veg_water(donc)
      cpatch%hcapveg(recc)    = cpatch%hcapveg(recc)    + cpatch%hcapveg(donc)
      if ( cpatch%hcapveg(recc) > 0. ) then !----- almost always the case. ----------------!
         !----- Updating temperature ------------------------------------------------------!
         call qwtk(cpatch%veg_energy(recc),cpatch%veg_water(recc),cpatch%hcapveg(recc)     &
                  ,cpatch%veg_temp(recc),cpatch%veg_fliq(recc))
      else !---- Veg_temp cannot be found using qwtk, it is a singularity. ----------------!
         cpatch%veg_temp(recc)  = newni *                                                  &
                                ( cpatch%veg_temp(recc)  * cpatch%nplant(recc)             &
                                + cpatch%veg_temp(donc)  * cpatch%nplant(donc))
         cpatch%veg_fliq(recc)  = 0.0
      end if

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

      cpatch%mean_leaf_resp(recc) = cpatch%mean_leaf_resp(recc)                            &
                                  + cpatch%mean_leaf_resp(donc)

      cpatch%mean_root_resp(recc) = cpatch%mean_root_resp(recc)                            &
                                  + cpatch%mean_root_resp(donc)

      cpatch%dmean_gpp(recc)     = cpatch%dmean_gpp(recc)                                  &
                                 + cpatch%dmean_gpp(donc)

      cpatch%dmean_gpp_pot(recc) = cpatch%dmean_gpp_pot(recc)                              &
                                 + cpatch%dmean_gpp_pot(donc)

      cpatch%dmean_gpp_max(recc) = cpatch%dmean_gpp_max(recc)                              &
                                 + cpatch%dmean_gpp_max(donc)

      cpatch%dmean_leaf_resp(recc) = cpatch%dmean_leaf_resp(recc)                          &
                                  + cpatch%dmean_leaf_resp(donc)

      cpatch%dmean_root_resp(recc) = cpatch%dmean_root_resp(recc)                          &
                                   + cpatch%dmean_root_resp(donc)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Not sure about the following variables.  From ed_state_vars, I would say that   !
      ! they should be averaged, not added because there it's written that these are per   !
      ! plant.  But from the fuse_2_patches_ar subroutine here it seems they are per unit  !
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

      !------ Psi_open is in kg/m2/s, so we add them. -------------------------------------!
      cpatch%Psi_open(recc)   = cpatch%Psi_open(recc)   + cpatch%Psi_open(donc)
      cpatch%Psi_closed(recc) = cpatch%Psi_closed(recc) + cpatch%Psi_closed(donc)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !   This one wasn't here originally, shouldn't it be? It is in split_cohorts, so my  !
      ! impression is that this should be always combined and rescaled everywhere in the   !
      ! patch and cohort fusion.  Not sure of how relevant this is...               MLO.   !
      !------------------------------------------------------------------------------------!
      cpatch%monthly_dndt(recc)        = cpatch%monthly_dndt(recc)                         &
                                       + cpatch%monthly_dndt(donc)
      !------------------------------------------------------------------------------------!
     
      cpatch%fsn(recc) = ( cpatch%fsn(recc) * cpatch%nplant(recc)                          &
                         + cpatch%fsn(donc) * cpatch%nplant(donc) ) * newni
     
      cpatch%fsw(recc) = ( cpatch%fsw(recc) * cpatch%nplant(recc)                          &
                         + cpatch%fsw(donc) * cpatch%nplant(donc) ) * newni


      !------------------------------------------------------------------------------------!
      !     Updating the carbon fluxes. They are fluxes per unit of area, so they should   !
      ! be added, not scaled.                                                              !
      !------------------------------------------------------------------------------------!
      cpatch%gpp(recc) = cpatch%gpp(recc) + cpatch%gpp(donc)

      cpatch%leaf_respiration(recc) = cpatch%leaf_respiration(recc)                        &
                                    + cpatch%leaf_respiration(donc)
      cpatch%root_respiration(recc) = cpatch%root_respiration(recc)                        &
                                    + cpatch%root_respiration(donc)
      !------------------------------------------------------------------------------------!
     
      cpatch%paw_avg(recc) = cpatch%paw_avg(recc) + cpatch%paw_avg(donc)

      cpatch%turnover_amp(recc)  = (cpatch%turnover_amp(recc) * cpatch%nplant(recc)        &
           + cpatch%turnover_amp(donc) * cpatch%nplant(donc) ) *newni   
      cpatch%llspan(recc)  = (cpatch%llspan(recc) * cpatch%nplant(recc)                    &
           + cpatch%llspan(donc) * cpatch%nplant(donc) ) *newni   
      cpatch%vm_bar(recc)  = (cpatch%vm_bar(recc) * cpatch%nplant(recc)                    &
           + cpatch%vm_bar(donc) * cpatch%nplant(donc) ) *newni   
      cpatch%sla(recc)  = (cpatch%sla(recc) * cpatch%nplant(recc)                          &
           + cpatch%sla(donc) * cpatch%nplant(donc) ) *newni   
    
      root_depth = calc_root_depth(cpatch%hite(recc), cpatch%dbh(recc), cpatch%pft(recc))
      cpatch%krdepth(recc) = assign_root_depth(root_depth, lsl)

      !----- Last, but not the least, we update nplant ------------------------------------!
      cpatch%nplant(recc) = newn

      return
   end subroutine fuse_2_cohorts_ar
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
   subroutine fuse_patches_ar(cgrid,ifm)
     
      use ed_state_vars       , only :  edtype            & ! structure
                                      , polygontype       & ! structure
                                      , sitetype          & ! structure
                                      , patchtype         ! ! structure
      use fusion_fission_coms , only :  ff_ndbh           & ! intent(in)
                                      , ntol              & ! intent(in)
                                      , profile_tol       & ! intent(in)
                                      , pat_tolerance_max ! ! intent(in)
      use max_dims            , only :  n_pft             ! ! intent(in)
      use mem_sites           , only :  maxpatch          & ! intent(in)
                                      , maxcohort         ! ! intent(in)
      use ed_node_coms        , only :  mynum

      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(edtype)         , target      :: cgrid           ! Current grid
      integer              , intent(in)  :: ifm             ! Current grid index
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)    , pointer     :: cpoly           ! Current polygon
      type(sitetype)       , pointer     :: csite           ! Current site
      type(patchtype)      , pointer     :: cpatch          ! Current patch
      type(sitetype)       , pointer     :: tempsite        ! Temporary site
      logical, dimension(:), allocatable :: fuse_table      ! Flag: this will remain.
      real   , dimension(n_pft,ff_ndbh)  :: mean_nplant     ! Mean # of plants
      integer                            :: ipy,isi         ! Counters
      integer                            :: ipa,ico         ! Counters
      integer                            :: donp,recp       ! Counters
      integer                            :: ipft,idbh       ! Counters
      integer                            :: npatches_new    ! New # of patches
      integer                            :: npatches_old    ! Old # of patches
      logical                            :: fuse_flag       ! Flag: I will perform fusion.
      real                               :: diff            !
      real                               :: refv            !
      real                               :: norm            !
      real                               :: tolerance_mult  ! Multiplying factor for tol.
      real                               :: old_area        ! For area conservation check
      real                               :: new_area        ! For area conservation check
      real                               :: old_lai_tot     ! Old total LAI
      real                               :: old_nplant_tot  ! Old total nplant
      real                               :: new_lai_tot     ! New total LAI
      real                               :: new_nplant_tot  ! New total nplant
      real                               :: elim_nplant     ! Elim. nplant during 1 fusion
      real                               :: elim_lai        ! Elim. LAI during 1 fusion
      real                               :: elim_nplant_tot ! Total eliminated nplant
      real                               :: elim_lai_tot    ! Elim. eliminated LAI
      integer                            :: tot_npolygons   ! Total # of polygons
      integer                            :: tot_nsites      ! Total # of sites
      integer                            :: tot_npatches    ! Total # of patches
      integer                            :: tot_ncohorts    ! Total # of cohorts
      !------------------------------------------------------------------------------------!

      !----- Return if maxpatch is 0, this is a flag for no patch fusion. -----------------!
      if (maxpatch == 0) return
      
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         
         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !----- Skip this site if it contains only one patch... ------------------------!
            if (csite%npatches < 2) cycle siteloop

            !----- Allocate the swapper patches in the site type and the fusion table. ----!
            nullify(tempsite)
            allocate(tempsite)
            call allocate_sitetype(tempsite, csite%npatches )
            allocate(fuse_table(csite%npatches))
            fuse_table(:) = .true.

            !----- This is for mass conservation check ------------------------------------!
            old_nplant_tot = 0.
            old_lai_tot    = 0.
            old_area       = 0.
            do ipa = 1,csite%npatches
               old_area  = old_area + csite%area(ipa)
               cpatch => csite%patch(ipa)
               do ico = 1, cpatch%ncohorts
                  old_nplant_tot = old_nplant_tot + cpatch%nplant(ico)*csite%area(ipa)
                  old_lai_tot    = old_lai_tot    + cpatch%lai(ico)*csite%area(ipa)
               end do
            end do

            !----- Initially we didn't eliminate any cohort... ----------------------------!
            elim_nplant_tot = 0.
            elim_lai_tot    = 0.

            !------------------------------------------------------------------------------!
            ! ALGORITHM                                                                    !
            !                                                                              !
            ! 1. Set all fusion flags to true                                              !
            ! 2. Create size profiles                                                      !
            ! 3. Go to every patch                                                         !
            ! 4. Find next older patch with same dist_type                                 !
            ! 5. Check fusion criterion. If within criterion, fuse, otherwise, skip        !
            ! 6. Loop from the youngest to oldest patch                                    !
            !------------------------------------------------------------------------------!
            if (csite%npatches > abs(maxpatch)) then


               mean_nplant = 0.0
               do ipa = csite%npatches,1,-1
                  call patch_pft_size_profile_ar(csite,ipa,ff_ndbh)
               
                  !------------------------------------------------------------------------!
                  !    Get a mean density profile for all of the patches. This will be     !
                  ! used for normalization.                                                !
                  !------------------------------------------------------------------------!
                  do ipft=1,n_pft
                     do idbh=1,ff_ndbh 
                        mean_nplant(ipft,idbh) = mean_nplant(ipft,idbh)                    &
                                               + csite%pft_density_profile(ipft,idbh,ipa)  &
                                               / real(csite%npatches)
                     end do
                  end do
               end do

               !----- Start with no multiplication factor. --------------------------------!
               tolerance_mult = 1.0
               max_patch: do
                  npatches_old = count(fuse_table)
               
                  !----- Loop from youngest to the second oldest patch --------------------!
                  do donp = csite%npatches,2,-1
                     cpatch => csite%patch(donp)

                     !----- If this patch was already merged, skip it. --------------------!
                     if (fuse_table(donp)) then
                        !------------------------------------------------------------------!
                        !    Cycle through the next patches and compare densities, but     !
                        ! only compare densities if the patches have the same disturbance  !
                        ! types. Of course, only existing patches (i.e. that weren't       !
                        ! merged yet) are compared.                                        !
                        !------------------------------------------------------------------!
                        next_patch: do recp = donp-1,1,-1

                           if ( csite%dist_type(donp) == csite%dist_type(recp)             &
                              .and. fuse_table(recp) ) then
                           
                              !------------------------------------------------------------!
                              !    Once we have identified the patch with the same         !
                              ! disturbance type and closest age (recp), determine if      !
                              ! it is similar enough to average (fuse) the two together.   !
                              !------------------------------------------------------------!
                              fuse_flag = .true.

                              !------------------------------------------------------------!
                              !     Testing.  If two patches are empty, I guess it's fine  !
                              ! to just fuse them.                                         !
                              !------------------------------------------------------------!
                              if (csite%patch(donp)%ncohorts > 0 .or.                      &
                                  csite%patch(recp)%ncohorts > 0) then
                                 !-----  Fusion criterion. --------------------------------!
                                 fuseloop:do ipft=1,n_pft
                                    do idbh=1,ff_ndbh

                                       if (csite%pft_density_profile(ipft,idbh,donp) >     &
                                           tolerance_mult*ntol                        .or. &
                                           csite%pft_density_profile(ipft,idbh,recp) >     &
                                           tolerance_mult*ntol                       ) then
                                          !------------------------------------------------!
                                          !     This is the normalized difference in their !
                                          ! biodensity profiles. If the normalized         !
                                          ! difference is greater than the tolerance for   !
                                          ! any of the pfts and dbh classes, then reject   !
                                          ! them as similar.                               !
                                          !                                                !
                                          ! Note: If one of the patches is missing any     !
                                          !       member of the profile it will force the  !
                                          !       norm to be 2.0. That is the highest the  !
                                          !       norm should be able to go.               !
                                          !------------------------------------------------!
                                          diff = abs(                                      &
                                                 csite%pft_density_profile(ipft,idbh,donp) &
                                               - csite%pft_density_profile(ipft,idbh,recp) )
                                          refv = 0.5                                       &
                                               *(csite%pft_density_profile(ipft,idbh,donp) &
                                               + csite%pft_density_profile(ipft,idbh,recp) )
                                          !norm = diff / mean_nplant(ipft,idbh)
                                          norm = diff / refv

                                          if (norm > profile_tol) then
                                             fuse_flag = .false.   ! reject
                                             exit fuseloop
                                          end if
                                       end if
                                    end do
                                 end do fuseloop
                              end if

                              !----- Create a mapping of the patches that fuse together. --!
                              if (fuse_flag) then

                                 !---------------------------------------------------------!
                                 !     Take an average of the patch properties at index    !
                                 ! donp and ipa_tp assign the average to index ipa_tp.     !
                                 !---------------------------------------------------------!
                                 call fuse_2_patches_ar(csite,donp,recp                    &
                                                       ,cpoly%met(isi)%rhos,cpoly%lsl(isi) &
                                                       ,cpoly%green_leaf_factor(:,isi)     &
                                                       ,elim_nplant,elim_lai)

                                 !----- Updating total eliminated nplant and LAI  ---------!
                                 elim_nplant_tot = elim_nplant_tot                         &
                                                 + elim_nplant * csite%area(recp)
                                 elim_lai_tot    = elim_lai_tot                            &
                                                 + elim_lai    * csite%area(recp)

                                 !---------------------------------------------------------!
                                 !     Recalculate the pft size profile for the averaged   !
                                 ! patch at donp_tp.                                       !
                                 !---------------------------------------------------------!
                                 call patch_pft_size_profile_ar(csite,recp,ff_ndbh)

                                 !---------------------------------------------------------!
                                 !     The patch at index donp is no longer valid, it      !
                                 ! should be flagged as such.                              !
                                 !---------------------------------------------------------!
                                 fuse_table(donp) = .false.

                                 !---------------------------------------------------------!
                                 !     If we have gotten to this point, we have found our  !
                                 ! donor patch and have performed the fusion. Exit the     !
                                 ! patch loop.                                             !
                                 !---------------------------------------------------------!
                                 exit next_patch
                              end if ! if( fuse_flag)
                           end if ! if(csite%dist_type(donp) == csite%dist_type(recp)...
                        end do next_patch       ! do recp
                     end if          ! if (.not. fuse_table(donp)) then

                     npatches_new = count(fuse_table)
                     if (npatches_new <= abs(maxpatch)) exit max_patch
                  end do          ! do donp = csite%npatches,2,-1

                  !------------------------------------------------------------------------!
                  !    If no fusion happened and it exceed the maximum tolerance, give up. !
                  !------------------------------------------------------------------------!
                  npatches_new = count(fuse_table)
                  if(npatches_new == npatches_old   .and.                                  &
                     tolerance_mult > pat_tolerance_max ) exit max_patch
                  
                  !----- Increment tolerance ----------------------------------------------!
                  tolerance_mult = tolerance_mult * 1.1
               end do max_patch
     
               !----- Set the number of patches in the site to "npatches_new" -------------!
               tempsite%npatches = npatches_new

               !----- If there was any patch fusion, need to shrink csite -----------------!
               if (npatches_new < csite%npatches) then
                  !------------------------------------------------------------------------!
                  !    Copy the selected data into the temporary space, args 1 and 3 must  !
                  ! be dimension of arg 4. Argument 2 must be the dimension of the sum of  !
                  ! the 3rd argument.                                                      !
                  !------------------------------------------------------------------------!
                  call copy_sitetype_mask(csite,tempsite,fuse_table,size(fuse_table)       &
                                         ,npatches_new)
                  call deallocate_sitetype(csite)

                  !----- Reallocate the current site. -------------------------------------!
                  call allocate_sitetype(csite,npatches_new)

                  !----- Copy the selected temporary data into the orignal site vectors. --!
                  fuse_table(:)              = .false.
                  fuse_table(1:npatches_new) = .true.
                  call copy_sitetype_mask(tempsite,csite,fuse_table,size(fuse_table)       &
                                         ,npatches_new)
                  !------------------------------------------------------------------------!
                  !     The new and fused csite is now complete, clean up the temporary    !
                  ! data. Deallocate it afterwards.                                        !
                  !------------------------------------------------------------------------!
                  call deallocate_sitetype(tempsite)

               end if
            end if
            !----- Deallocation should happen outside the "if" statement ------------------!
            deallocate(tempsite)
            deallocate(fuse_table)
            
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
            if (new_area       < 0.99 * old_area       .or.                                &
                new_area       > 1.01 * old_area       .or.                                &
                new_nplant_tot < 0.99 * (old_nplant_tot - elim_nplant_tot) .or.            &
                new_nplant_tot > 1.01 * (old_nplant_tot - elim_nplant_tot) .or.            &
                new_lai_tot    < 0.99 * (old_lai_tot    - elim_lai_tot   ) .or.            &
                new_lai_tot    > 1.01 * (old_lai_tot    - elim_lai_tot   )     ) then
               write (unit=*,fmt='(a,1x,es12.5)') 'OLD_AREA:       ',old_area
               write (unit=*,fmt='(a,1x,es12.5)') 'NEW_AREA:       ',new_area
               write (unit=*,fmt='(a,1x,es12.5)') 'NEW_LAI_TOT:    ',new_lai_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'OLD_LAI_TOT:    ',old_lai_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'ELIM_LAI_TOT:   ',elim_lai_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'NEW_NPLANT_TOT: ',new_nplant_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'OLD_NPLANT_TOT: ',old_nplant_tot
               write (unit=*,fmt='(a,1x,es12.5)') 'ELIM_NPLANT_TOT:',elim_nplant_tot
               call fatal_error('Conservation failed while fusing patches'                 &
                              &,'fuse_patches_ar','fuse_fiss_utils.f90')
            end if
            
         end do siteloop
      end do polyloop
      
      !------------------------------------------------------------------------------------!
      !     Printing a banner to inform the user how many patches and cohorts exist.  To   !
      ! avoid dumping too many information, display the message only when the user is not  !
      ! running regional runs.                                                             !
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
   end subroutine fuse_patches_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine will merge two patches into 1.                                      !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_patches_ar(csite,donp,recp,rhos,lsl,green_leaf_factor                 &
                               ,elim_nplant,elim_lai)
      use ed_state_vars      , only :  sitetype              & ! Structure 
                                     , patchtype             ! ! Structure
      use soil_coms          , only :  soil                  & ! intent(in), lookup table
                                     , min_sfcwater_mass     ! ! intent(in)
      use grid_coms          , only :  nzg                   & ! intent(in)
                                     , nzs                   ! ! intent(in)
      use fusion_fission_coms, only :  ff_ndbh               ! ! intent(in)
      use max_dims           , only :  n_pft                 & ! intent(in)
                                     , n_dbh                 ! ! intent(in)
      use mem_sites          , only :  maxcohort             ! ! intent(in)
      use consts_coms        , only :  cpi                   & ! intent(in)
                                     , cpor                  & ! intent(in)
                                     , p00                   ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target      :: csite             ! Current site
      integer                , intent(in)  :: donp              ! Donating patch
      integer                , intent(in)  :: recp              ! Receptor patch
      integer                , intent(in)  :: lsl               ! Lowest soil level
      real, dimension(n_pft) , intent(in)  :: green_leaf_factor ! Green leaf factor...
      real                   , intent(in)  :: rhos              ! Sfc. air density
      real                   , intent(out) :: elim_nplant       ! Eliminated nplant 
      real                   , intent(out) :: elim_lai          ! Eliminated lai
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer     :: cpatch            ! Current patch
      type(patchtype)        , pointer     :: temppatch         ! Temporary patch
      integer                              :: ico,iii           ! Counters
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

      csite%can_temp(recp)           = newareai *                                          &
                                     ( csite%can_temp(donp)           * csite%area(donp)   &
                                     + csite%can_temp(recp)           * csite%area(recp) )

      csite%can_shv(recp)            = newareai *                                          &
                                     ( csite%can_shv(donp)            * csite%area(donp)   &
                                     + csite%can_shv(recp)            * csite%area(recp) )

      csite%hcapveg(recp)            = newareai *                                          &
                                     ( csite%hcapveg(donp)            * csite%area(donp)   &
                                     + csite%hcapveg(recp)            * csite%area(recp) )


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
      !    (Both are done in new_patch_sfc_props_ar).                                      !
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!


      !----- Merging soil energy and water. -----------------------------------------------!
      do iii=1,nzg
         csite%soil_energy(iii,recp)     = newareai *                                      &
                                         ( csite%soil_energy(iii,donp) * csite%area(donp)  &
                                         + csite%soil_energy(iii,recp) * csite%area(recp))

         csite%soil_water(iii,recp)      = newareai *                                      &
                                         ( csite%soil_water(iii,recp)  * csite%area(recp)  &
                                         + csite%soil_water(iii,donp)  * csite%area(donp))
      end do

      !------------------------------------------------------------------------------------!
      !    This subroutine takes care of filling:                                          !
      !                                                                                    !
      ! + csite%ground_shv(recp)                                                           !
      ! + csite%surface_ssh(recp)                                                          !
      ! + csite%soil_tempk(k,recp)                                                         !
      ! + csite%soil_fracliq(k,recp)                                                       !
      ! + csite%nlev_sfcwater(recp)                                                        !
      ! + csite%sfcwater_energy(k,recp) (Just converting back to J/kg)                     !
      ! + csite%csite%sfcwater_tempk(k,recp)                                               !
      ! + csite%sfcwater_fracliq(k,recp)                                                   !
      !------------------------------------------------------------------------------------!
      call new_patch_sfc_props_ar(csite,recp,rhos)
      !------------------------------------------------------------------------------------!

      csite%mean_rh(recp)                = newareai *                                      &
                                         ( csite%mean_rh(donp)        * csite%area(donp)   &
                                         + csite%mean_rh(recp)        * csite%area(recp) )

      csite%dmean_A_decomp(recp)         = newareai *                                      &
                                         ( csite%dmean_A_decomp(donp) * csite%area(donp)   &
                                         + csite%dmean_A_decomp(recp) * csite%area(recp) )

      csite%dmean_Af_decomp(recp)        = newareai *                                      &
                                         ( csite%dmean_Af_decomp(donp)* csite%area(donp)   &
                                         + csite%dmean_Af_decomp(recp)* csite%area(recp) )

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
      csite%avg_carbon_ac(recp)       = newareai *                                         &
                                      ( csite%avg_carbon_ac(donp)     * csite%area(donp)   &
                                      + csite%avg_carbon_ac(recp)     * csite%area(recp) )

      csite%avg_vapor_vc(recp)        = newareai *                                         &
                                      ( csite%avg_vapor_vc(donp)      * csite%area(donp)   &
                                      + csite%avg_vapor_vc(recp)      * csite%area(recp) )  

      csite%avg_dew_cg(recp)          = newareai *                                         &
                                      ( csite%avg_dew_cg(donp)        * csite%area(donp)   &
                                      + csite%avg_dew_cg(recp)        * csite%area(recp) )  

      csite%avg_vapor_gc(recp)        = newareai *                                         &
                                      ( csite%avg_vapor_gc(donp)      * csite%area(donp)   &
                                      + csite%avg_vapor_gc(recp)      * csite%area(recp) )  

      csite%avg_wshed_vg(recp)        = newareai *                                         &
                                      ( csite%avg_wshed_vg(donp)      * csite%area(donp)   &
                                      + csite%avg_wshed_vg(recp)      * csite%area(recp) )  

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

      csite%aux(recp)                 = newareai *                                         &
                                      ( csite%aux(donp)               * csite%area(donp)   &
                                      + csite%aux(recp)               * csite%area(recp) )  

      csite%avg_sensible_vc(recp)     = newareai *                                         &
                                      ( csite%avg_sensible_vc(donp)   * csite%area(donp)   &
                                      + csite%avg_sensible_vc(recp)   * csite%area(recp) )  

      csite%avg_sensible_2cas(recp)   = newareai *                                         &
                                      ( csite%avg_sensible_2cas(donp) * csite%area(donp)   &
                                      + csite%avg_sensible_2cas(recp) * csite%area(recp) )  

      csite%avg_qwshed_vg(recp)       = newareai *                                         &
                                      ( csite%avg_qwshed_vg(donp)     * csite%area(donp)   &
                                      + csite%avg_qwshed_vg(recp)     * csite%area(recp) )  

      csite%avg_sensible_gc(recp)     = newareai *                                         &
                                      ( csite%avg_sensible_gc(donp)   * csite%area(donp)   &
                                      + csite%avg_sensible_gc(recp)   * csite%area(recp) )  

      csite%avg_sensible_ac(recp)     = newareai *                                         &
                                      ( csite%avg_sensible_ac(donp)   * csite%area(donp)   &
                                      + csite%avg_sensible_ac(recp)   * csite%area(recp) )  

      csite%avg_sensible_tot(recp)    = newareai *                                         &
                                      ( csite%avg_sensible_tot(donp)  * csite%area(donp)   &
                                      + csite%avg_sensible_tot(recp)  * csite%area(recp) )  

      csite%avg_runoff_heat(recp)     = newareai *                                         &
                                      ( csite%avg_runoff_heat(donp)   * csite%area(donp)   &
                                      + csite%avg_runoff_heat(recp)   * csite%area(recp) )  

      csite%avg_heatstor_veg(recp)    = newareai *                                         &
                                      ( csite%avg_heatstor_veg(donp)  * csite%area(donp)   &
                                      + csite%avg_heatstor_veg(recp)  * csite%area(recp) )  

      csite%avg_veg_energy(recp)      = newareai *                                         &
                                      ( csite%avg_veg_energy(donp)    * csite%area(donp)   &
                                      + csite%avg_veg_energy(recp)    * csite%area(recp) )  

      csite%avg_veg_temp(recp)        = newareai *                                         &
                                      ( csite%avg_veg_temp(donp)      * csite%area(donp)   &
                                      + csite%avg_veg_temp(recp)      * csite%area(recp) )  

      csite%avg_veg_fliq(recp)        = newareai *                                         &
                                      ( csite%avg_veg_fliq(donp)      * csite%area(donp)   &
                                      + csite%avg_veg_fliq(recp)      * csite%area(recp) )  

      csite%avg_veg_water(recp)       = newareai *                                         &
                                      ( csite%avg_veg_water(donp)     * csite%area(donp)   &
                                      + csite%avg_veg_water(recp)     * csite%area(recp) )  

      csite%co2budget_gpp(recp)       = newareai *                                         &
                                      ( csite%co2budget_gpp(donp)     * csite%area(donp)   &
                                      + csite%co2budget_gpp(recp)     * csite%area(recp) )  

      csite%co2budget_plresp(recp)    = newareai *                                         &
                                      ( csite%co2budget_plresp(donp)  * csite%area(donp)   &
                                      + csite%co2budget_plresp(recp)  * csite%area(recp) )  

      csite%co2budget_rh(recp)        = newareai *                                         &
                                      ( csite%co2budget_rh(donp)      * csite%area(donp)   &
                                      + csite%co2budget_rh(recp)      * csite%area(recp) )  

      do iii=1,nzg
         csite%avg_smoist_gg(iii,recp)   = newareai *                                      &
              ( csite%avg_smoist_gg(iii,donp)       * csite%area(donp)                     &
              + csite%avg_smoist_gg(iii,recp)       * csite%area(recp) )

         csite%avg_smoist_gc(iii,recp)   = newareai *                                      &
              ( csite%avg_smoist_gc(iii,donp)       * csite%area(donp)                     &
              + csite%avg_smoist_gc(iii,recp)       * csite%area(recp) )

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
      !     Now we need to adjust the densities of cohorts. Because the patch area         !
      ! increased we want to retain the same total amount of mass and energy.              !
      !------------------------------------------------------------------------------------!
      !----- 1. Adjust densities of cohorts in recipient patch ----------------------------!
      cpatch => csite%patch(recp)
      nrc = cpatch%ncohorts
      area_scale = csite%area(recp) * newareai
      do ico = 1,nrc
         cpatch%lai(ico)                 = cpatch%lai(ico)                  * area_scale
         cpatch%bai(ico)                 = cpatch%bai(ico)                  * area_scale
         cpatch%sai(ico)                 = cpatch%sai(ico)                  * area_scale
         cpatch%nplant(ico)              = cpatch%nplant(ico)               * area_scale
         cpatch%mean_gpp(ico)            = cpatch%mean_gpp(ico)             * area_scale
         cpatch%mean_leaf_resp(ico)      = cpatch%mean_leaf_resp(ico)       * area_scale
         cpatch%mean_root_resp(ico)      = cpatch%mean_root_resp(ico)       * area_scale
         cpatch%dmean_gpp(ico)           = cpatch%dmean_gpp(ico)            * area_scale
         cpatch%dmean_gpp_pot(ico)       = cpatch%dmean_gpp_pot(ico)        * area_scale
         cpatch%dmean_gpp_max(ico)       = cpatch%dmean_gpp_max(ico)        * area_scale
         cpatch%dmean_leaf_resp(ico)     = cpatch%dmean_leaf_resp(ico)      * area_scale
         cpatch%dmean_root_resp(ico)     = cpatch%dmean_root_resp(ico)      * area_scale
         cpatch%Psi_open(ico)            = cpatch%Psi_open(ico)             * area_scale
         cpatch%gpp(ico)                 = cpatch%gpp(ico)                  * area_scale
         cpatch%leaf_respiration(ico)    = cpatch%leaf_respiration(ico)     * area_scale
         cpatch%root_respiration(ico)    = cpatch%root_respiration(ico)     * area_scale
         cpatch%monthly_dndt(ico)        = cpatch%monthly_dndt(ico)         * area_scale
         cpatch%veg_water(ico)           = cpatch%veg_water(ico)            * area_scale
         cpatch%hcapveg(ico)             = cpatch%hcapveg(ico)              * area_scale
         cpatch%veg_energy(ico)          = cpatch%veg_energy(ico)           * area_scale
      end do
      !----- 2. Adjust densities of cohorts in donor patch --------------------------------!
      cpatch => csite%patch(donp)
      ndc = cpatch%ncohorts
      area_scale = csite%area(donp) * newareai
      do ico = 1,ndc
         cpatch%lai(ico)                 = cpatch%lai(ico)                  * area_scale
         cpatch%bai(ico)                 = cpatch%bai(ico)                  * area_scale
         cpatch%sai(ico)                 = cpatch%sai(ico)                  * area_scale
         cpatch%nplant(ico)              = cpatch%nplant(ico)               * area_scale
         cpatch%mean_gpp(ico)            = cpatch%mean_gpp(ico)             * area_scale
         cpatch%mean_leaf_resp(ico)      = cpatch%mean_leaf_resp(ico)       * area_scale
         cpatch%mean_root_resp(ico)      = cpatch%mean_root_resp(ico)       * area_scale
         cpatch%dmean_gpp(ico)           = cpatch%dmean_gpp(ico)            * area_scale
         cpatch%dmean_gpp_pot(ico)       = cpatch%dmean_gpp_pot(ico)        * area_scale
         cpatch%dmean_gpp_max(ico)       = cpatch%dmean_gpp_max(ico)        * area_scale
         cpatch%dmean_leaf_resp(ico)     = cpatch%dmean_leaf_resp(ico)      * area_scale
         cpatch%dmean_root_resp(ico)     = cpatch%dmean_root_resp(ico)      * area_scale
         cpatch%Psi_open(ico)            = cpatch%Psi_open(ico)             * area_scale
         cpatch%gpp(ico)                 = cpatch%gpp(ico)                  * area_scale
         cpatch%leaf_respiration(ico)    = cpatch%leaf_respiration(ico)     * area_scale
         cpatch%root_respiration(ico)    = cpatch%root_respiration(ico)     * area_scale
         cpatch%monthly_dndt(ico)        = cpatch%monthly_dndt(ico)         * area_scale
         cpatch%veg_water(ico)           = cpatch%veg_water(ico)            * area_scale
         cpatch%hcapveg(ico)             = cpatch%hcapveg(ico)              * area_scale
         cpatch%veg_energy(ico)          = cpatch%veg_energy(ico)           * area_scale
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
         call sort_cohorts_ar(cpatch)
         !---------------------------------------------------------------------------------!
         !    We just combined two patches, so we may be able to fuse some cohorts and/or  !
         ! eliminate others.                                                               !
         !---------------------------------------------------------------------------------!
         if (cpatch%ncohorts > 0 .and. maxcohort >= 0) then
            call fuse_cohorts_ar(csite,recp,green_leaf_factor,lsl)
            call terminate_cohorts_ar(csite,recp,elim_nplant,elim_lai)
            call split_cohorts_ar(cpatch,green_leaf_factor,lsl)
         end if
         !---------------------------------------------------------------------------------!
      end if

      !------------------------------------------------------------------------------------!
      !    Now we update some variables that depend on cohort statistics, namely:          !
      ! + csite%veg_height(recp)                                                           !
      ! + csite%lai(recp)                                                                  !
      ! + csite%veg_rough(recp)                                                            !
      ! + csite%wbudget_initialstorage(recp)                                               !
      ! + csite%ebudget_initialstorage(recp)                                               !
      ! + csite%co2budget_initialstorage(recp)                                             !
      !------------------------------------------------------------------------------------!
      call update_patch_derived_props_ar(csite,lsl, rhos,recp)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    This subroutine will update the size profile within patch.                      !
      ! + csite%pft_density_profile(:,:,recp)                                              !
      !------------------------------------------------------------------------------------!
      call patch_pft_size_profile_ar(csite,recp,ff_ndbh)
      !------------------------------------------------------------------------------------!

      !----- Last, but not the least, we update the patch area ----------------------------!
      csite%area(recp) = newarea

      return

   end subroutine fuse_2_patches_ar
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine patch_pft_size_profile_ar(csite,ipa,nbins)
      use ed_state_vars      , only :  sitetype   & ! structure
                                     , patchtype  ! ! structure
      use fusion_fission_coms, only :  maxdbh     ! ! intent(in)
      use max_dims           , only :  n_pft      ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target     :: csite             ! Current site
      integer                , intent(in) :: ipa               ! Current patch ID
      integer                , intent(in) :: nbins             ! # of DBH classes
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer    :: cpatch        ! Current patch
      integer                             :: ipft,idbh,ico ! Counters
      real                                :: ddbh          ! Class interval size
      !------------------------------------------------------------------------------------!

      !----- Finding the size of each DBH class interval ----------------------------------!
      ddbh = maxdbh/real(nbins)

      !----- Initialize bins --------------------------------------------------------------!
      do ipft=1,n_pft
         do idbh=1,nbins
            csite%pft_density_profile(ipft,idbh,ipa)=0.0
         end do
      end do

      !----- Update bins ------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      do ico = 1,cpatch%ncohorts

         ipft = cpatch%pft(ico)
         idbh = min(nbins,max(1,ceiling(cpatch%dbh(ico)/ddbh)))
         csite%pft_density_profile(ipft,idbh,ipa) = cpatch%nplant(ico)                     &
                                                  + csite%pft_density_profile(ipft,idbh,ipa)
      end do

      return
   end subroutine patch_pft_size_profile_ar
   !=======================================================================================!
   !=======================================================================================!
end module fuse_fiss_utils_ar
!==========================================================================================!
!==========================================================================================!
