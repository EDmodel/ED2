!==========================================================================================!
!==========================================================================================!
! MODULE: FUSE_FISS_UTILS
!> \brief Routines to terminate and fuse patches, and to terminate, fuse and split cohorts.
!> \author  Translated from ED1 by David Medivgy, Ryan Knox and Marcos Longo
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
   !  SUBROUTINE: SORT_COHORTS      
   !> \brief This subroutine will sort the cohorts by size (1st = tallest, last = shortest.)
   !> \details In case there is a tie (for example, when 2 cohorts have reached the
   !> maximum possible height, then we use DBH for tie breaking, and if they have the
   !> exact same DBH, then we simply pick the lowest index (as they are exactly the same).
   !> This could cause some problems when the new grass allometry is implemented, though.       
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


      !------------------------------------------------------------------------------------!
      !     Check whether this patch is already sorted.   We don't want to do the entire   !
      ! deallocating/copying/allocating thing if it's not needed as this takes up too much !
      ! time.                                                                              !
      !------------------------------------------------------------------------------------!
      sorted = .true.
      sortcheck: do ico=1,cpatch%ncohorts-1
         sorted = cpatch%height(ico) >= cpatch%dbh(ico+1) .and.                            &
                  cpatch%dbh(ico)  >= cpatch%dbh(ico+1)
         if (.not. sorted) exit sortcheck
      end do sortcheck
      if (sorted) return
      !------------------------------------------------------------------------------------!



      !----- Assign a scratch patch. ------------------------------------------------------!
      nullify(temppatch)
      allocate(temppatch)
      call allocate_patchtype(temppatch,cpatch%ncohorts)
      !------------------------------------------------------------------------------------!

      !----- Allocate the logical flag for tie-breaking. ----------------------------------!
      allocate(attop(cpatch%ncohorts))
      !------------------------------------------------------------------------------------!

      ico = 0
      !---- Loop until all cohorts were sorted. -------------------------------------------!
      do while(ico < cpatch%ncohorts)
         ico = ico + 1

         !----- Find the maximum height. --------------------------------------------------!
         tophgt = maxval(cpatch%height)

         !----- Find all cohorts that are at this height. ---------------------------------!
         attop  = cpatch%height == tophgt

         !----- Find the fattest cohort at a given height. --------------------------------!
         tallco = maxloc(cpatch%dbh,dim=1,mask=attop)

         !----- Copy to the scratch structure. --------------------------------------------!
         call copy_patchtype(cpatch,temppatch,tallco,tallco,ico,ico)

         !----- Put a non-sense DBH so this will never "win" again. -----------------------!
         cpatch%height(tallco) = -huge(1.)
         cpatch%dbh   (tallco) = -huge(1.)

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
   !  SUBROUTINE: TERMINATE_COHORTS 
   !> \brief This subroutine will eliminate cohorts based on their sizes.
   !> This is intended to eliminate cohorts that have little contribution
   !> and thus we can speed up the run.
   !---------------------------------------------------------------------------------------!
   subroutine terminate_cohorts(csite,ipa,cmet,is_initial,elim_nplant,elim_lai)
      use pft_coms           , only : min_cohort_size     & ! intent(in)
                                    , l2n_stem            & ! intent(in)
                                    , c2n_stem            & ! intent(in)
                                    , c2n_storage         & ! intent(in)
                                    , c2n_leaf            & ! intent(in)
                                    , f_labile_leaf       & ! intent(in)
                                    , f_labile_stem       & ! intent(in)
                                    , agf_bs              ! ! intent(in)
      use fusion_fission_coms, only : print_fuse_details  ! ! intent(in)
      use ed_misc_coms       , only : current_time        & ! intent(in)
                                    , frqsumi             ! ! intent(in)
      use rk4_coms           , only : checkbudget         ! ! intent(in)
      use ed_state_vars      , only : patchtype           & ! structure
                                    , sitetype            ! ! structure
      use met_driver_coms    , only : met_driv_state      ! ! structure
      use therm_lib          , only : tq2enthalpy         & ! function
                                    , idealdenssh         & ! function
                                    , reducedpress        ! ! function
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)       , target      :: csite          ! Current site
      integer              , intent(in)  :: ipa            ! Current patch ID
      type(met_driv_state) , target      :: cmet           ! Current met forcing
      logical              , intent(in)  :: is_initial     ! Initial call?
      real                 , intent(out) :: elim_nplant    ! Nplants eliminated here
      real                 , intent(out) :: elim_lai       ! LAI eliminated here
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)      , pointer     :: cpatch         ! Current patch
      type(patchtype)      , pointer     :: temppatch      ! Scratch patch structure
      logical, dimension(:), allocatable :: remain_table   ! Flag: this cohort will remain
      logical                            :: is_tiny        ! Cohort is too tiny.
      integer                            :: ico            ! Counter
      integer                            :: ipft           ! PFT size
      real                               :: csize          ! Size of current cohort
      real                               :: elim_e_hcap    ! Heat capacity loss  (energy)
      real                               :: elim_e_wcap    ! Water capacity loss (energy)
      real                               :: elim_w_wcap    ! Water capacity loss (water )
      real                               :: veg_energy_tot ! Total internal energy (cohort)
      real                               :: veg_energy_im2 ! Internal water energy (cohort)
      real                               :: veg_water_im2  ! Internal water mass (cohort)
      real                               :: veg_qboil      ! Leaf energy transferred to CAS
      real                               :: veg_boil_tot   ! Water transferred to CAS
      real                               :: can_prss       ! CAS pressure
      real                               :: can_rhos       ! CAS density
      !------ Debugging variables. --------------------------------------------------------!
      logical          , save      :: first_time    = .true.
      character(len=18), parameter :: terminate_log = 'end_cohort_log.txt'
      character(len=10), parameter :: fmti          = '(a,1x,i14)'
      character(len=13), parameter :: fmtf          = '(a,1x,es12.5)'
      character(len=12), parameter :: fmth          = '(a,11(1x,a))'
      character(len=27), parameter :: fmtt          = '(a,i4.4,2(1x,i2.2),1x,f6.0)'
      character(len=35), parameter :: fmtc          = '(i12,1x,i12,2(11x,l1),8(1x,es12.5))'
      !------------------------------------------------------------------------------------!

      cpatch        => csite%patch(ipa)
      elim_nplant   = 0.
      elim_lai      = 0.
      elim_e_hcap   = 0.
      elim_e_wcap   = 0.
      elim_w_wcap   = 0.
      veg_boil_tot  = 0.

      !----- Initialize the temporary patch structures and the remain/terminate table -----!
      nullify(temppatch)
      allocate(temppatch)
      allocate(remain_table(cpatch%ncohorts))
      remain_table(:) = .true.
      !------------------------------------------------------------------------------------!



      !----- Debugging message. -----------------------------------------------------------!
      if (first_time) then
         if (print_fuse_details) then
            open (unit=23,file=terminate_log,status='replace',action='write')
            write(unit=23,fmt='(a)') ' Cohort termination log'
            write(unit=23,fmt='(a)') ' '
            close(unit=23,status='keep')
         end if
         first_time = .false.
      end if
      !------------------------------------------------------------------------------------!



      !----- Debugging message. -----------------------------------------------------------!
      if (print_fuse_details) then
         open (unit=23,file=terminate_log,status='replace',action='write')
         write(unit=23,fmt='(a)') '======================================================='
         write(unit=23,fmt='(a)') '======================================================='
         write(unit=23,fmt=fmtt ) ' TIME               : ',current_time%year                &
                                                         ,current_time%month               &
                                                         ,current_time%date                &
                                                         ,current_time%time
         write(unit=23,fmt=fmti ) ' PATCH              : ',ipa
         write(unit=23,fmt=fmti ) ' DIST_TYPE          : ',csite%dist_type(ipa)
         write(unit=23,fmt=fmtf ) ' AGE                : ',csite%age      (ipa)
         write(unit=23,fmt='(a)') ' '
         write(unit=23,fmt='(a)') ' '
         write(unit=23,fmt='(a)') ' ----------------------------'
         write(unit=23,fmt='(a)') '  List of eliminated cohorts'
         write(unit=23,fmt='(a)') ' ----------------------------'

         write(unit=23,fmt=fmth ) '         ICO','        IPFT','     IS_TINY'             &
                                 ,' IS_INVIABLE','         LAI','      NPLANT'             &
                                 ,'      BALIVE','      BDEADA','      BDEADB'             &
                                 ,'    BSTORAGE','       CSIZE','MIN_COH_SIZE'

      end if
      !------------------------------------------------------------------------------------!



      !----- Main loop --------------------------------------------------------------------!
      do ico = 1,cpatch%ncohorts

         !----- Save the PFT type in a convenient alias. ----------------------------------!
         ipft = cpatch%pft(ico)

         !----- Check whether the cohort size is too small --------------------------------!
         csize = cpatch%nplant(ico)                                                        &
               * ( cpatch%balive(ico) + cpatch%bstorage(ico)                               &
                 + cpatch%bdeada(ico) + cpatch%bdeadb  (ico) )
         is_tiny = csize < min_cohort_size(ipft)
         !---------------------------------------------------------------------------------!

         if ( is_tiny .or. (.not. cpatch%is_viable(ico)) ) then
            !----- Cohort is indeed too small or it is not viable, terminate it. ----------!
            remain_table(ico) = .false.
            elim_nplant       = elim_nplant + cpatch%nplant(ico) * csite%area(ipa)
            elim_lai          = elim_lai    + cpatch%lai(ico)    * csite%area(ipa)
            !------------------------------------------------------------------------------!



            !----- Debugging message. -----------------------------------------------------!
            if (print_fuse_details) then
               write(unit=23,fmt=fmtc )  ico,ipft,is_tiny, .not. cpatch%is_viable(ico)     &
                                        ,cpatch%lai(ico),cpatch%nplant(ico)                &
                                        ,cpatch%balive(ico),cpatch%bdeada(ico)             &
                                        ,cpatch%bdeadb(ico),cpatch%bstorage(ico)           &
                                        ,csize,min_cohort_size(ipft)
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Check whether to eliminate internal and surface water (this is necessary !
            ! only if the cohort was flagged as resolvable).                               !
            !------------------------------------------------------------------------------!
            veg_water_im2  = 0.0
            veg_energy_im2 = 0.0
            !----- Leaves. ----------------------------------------------------------------!
            if (cpatch%leaf_resolvable(ico)) then
               veg_water_im2  = veg_water_im2  + cpatch%leaf_water_im2(ico)
               veg_energy_im2 = veg_energy_im2 + cpatch%leaf_water_im2(ico)                &
                              * tq2enthalpy(cpatch%leaf_temp(ico),1.0,.true.)
            end if
            !----- Wood. ------------------------------------------------------------------!
            if (cpatch%wood_resolvable(ico)) then
               veg_water_im2  = veg_water_im2  + cpatch%wood_water_im2(ico)
               veg_energy_im2 = veg_energy_im2 + cpatch%wood_water_im2(ico)                &
                              * tq2enthalpy(cpatch%wood_temp(ico),1.0,.true.)
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     In case the cohort has any intercepted water, exchange moisture with the !
            ! canopy air space by donating the total amount as "boiling" (fast evaporation !
            ! or sublimation).                                                             !
            !------------------------------------------------------------------------------!
            veg_qboil = 0.0
            !----- Leaves. ----------------------------------------------------------------!
            if (cpatch%leaf_resolvable(ico)) then
               veg_boil_tot = veg_boil_tot + cpatch%leaf_water(ico)
               veg_qboil    = veg_qboil    + cpatch%leaf_water(ico)                        &
                                           * tq2enthalpy(cpatch%leaf_temp(ico),1.0,.true.)
            end if
            !----- Wood. ------------------------------------------------------------------!
            if (cpatch%wood_resolvable(ico)) then
               veg_boil_tot = veg_boil_tot + cpatch%wood_water(ico)
               veg_qboil    = veg_qboil    + cpatch%wood_water(ico)                        &
                                           * tq2enthalpy(cpatch%wood_temp(ico),1.0,.true.)
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Account for loss of internal energy stored in oven-dry biomass.          !
            !------------------------------------------------------------------------------!
            veg_energy_tot = 0.0
            !----- Leaves. ----------------------------------------------------------------!
            if (cpatch%leaf_resolvable(ico)) then
               veg_energy_tot = veg_energy_tot + cpatch%leaf_energy(ico)
            end if
            !----- Wood. ------------------------------------------------------------------!
            if (cpatch%wood_resolvable(ico)) then
               veg_energy_tot = veg_energy_tot + cpatch%wood_energy(ico)
            end if
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Update the total water and energy losses due to termination.  Do this    !
            ! only if cohort was somehow still resolvable, otherwise the cohort was        !
            ! already out of the budget.                                                   !
            !------------------------------------------------------------------------------!
            elim_e_hcap = elim_e_hcap + veg_energy_tot - veg_energy_im2 - veg_qboil
            elim_e_wcap = elim_e_wcap + veg_energy_im2
            elim_w_wcap = elim_w_wcap + veg_water_im2
            !------------------------------------------------------------------------------!



            !----- Update litter pools ----------------------------------------------------!
            csite%fgc_in(ipa) = csite%fgc_in(ipa) + cpatch%nplant(ico)                     &
                              * ( agf_bs(ipft) * cpatch%bstorage(ico)                      &
                                + f_labile_leaf(ipft) * cpatch%bleaf(ico)                  &
                                + f_labile_stem(ipft)                                      &
                                * ( cpatch%bsapwooda(ico) + cpatch%bbarka   (ico)          &
                                  + cpatch%bdeada   (ico) ) )
            csite%fsc_in(ipa) = csite%fsc_in(ipa) + cpatch%nplant(ico)                     &
                              * ( (1.0 - agf_bs(ipft)) * cpatch%bstorage(ico)              &
                                + f_labile_leaf(ipft) * cpatch%broot(ico)                  &
                                + f_labile_stem(ipft)                                      &
                                * ( cpatch%bsapwoodb(ico) + cpatch%bbarkb   (ico)          &
                                  + cpatch%bdeadb   (ico) ) )

            csite%stgc_in(ipa) = csite%stgc_in(ipa) + cpatch%nplant(ico)                   &
                               * ( ( 1.0 - f_labile_leaf(ipft) ) * cpatch%bleaf(ico)       &
                                 + ( 1.0 - f_labile_stem(ipft) )                           &
                                 * ( cpatch%bsapwooda(ico) + cpatch%bbarka(ico)            &
                                   + cpatch%bdeada   (ico)                      ) )

            csite%stsc_in(ipa) = csite%stsc_in(ipa) + cpatch%nplant(ico)                   &
                               * ( ( 1.0 - f_labile_leaf(ipft) ) * cpatch%broot(ico)       &
                                 + ( 1.0 - f_labile_stem(ipft) )                           &
                                 * ( cpatch%bsapwoodb(ico) + cpatch%bbarkb(ico)            &
                                   + cpatch%bdeadb   (ico)                      ) )

            csite%stgl_in(ipa) = csite%stgl_in(ipa) + cpatch%nplant(ico)                   &
                               * ( ( 1.0 - f_labile_leaf(ipft) ) * cpatch%bleaf(ico)       &
                                 + ( 1.0 - f_labile_stem(ipft) )                           &
                                 * ( cpatch%bsapwooda(ico) + cpatch%bbarka   (ico)         &
                                   + cpatch%bdeada   (ico)                         ) )     &
                               * l2n_stem / c2n_stem(ipft)

            csite%stsl_in(ipa) = csite%stsl_in(ipa) + cpatch%nplant(ico)                   &
                               * ( ( 1.0 - f_labile_leaf(ipft) ) * cpatch%broot(ico)       &
                                 + ( 1.0 - f_labile_stem(ipft) )                           &
                                 * ( cpatch%bsapwoodb(ico) + cpatch%bbarkb   (ico)         &
                                   + cpatch%bdeadb   (ico)                         ) )     &
                               * l2n_stem / c2n_stem(ipft)

            csite%fgn_in(ipa) = csite%fgn_in(ipa) + cpatch%nplant(ico)                     &
                              * ( agf_bs(ipft) / c2n_storage * cpatch%bstorage(ico)        &
                                + f_labile_leaf(ipft) / c2n_leaf(ipft) * cpatch%bleaf(ico) &
                                + f_labile_stem(ipft) / c2n_leaf(ipft)                     &
                                * ( cpatch%bsapwooda(ico) + cpatch%bbarka(ico)             &
                                  + cpatch%bdeada   (ico)                      ) )

            csite%fsn_in(ipa) = csite%fsn_in(ipa) + cpatch%nplant(ico)                     &
                              * ( (1.0-agf_bs(ipft))  / c2n_storage * cpatch%bstorage(ico) &
                                + f_labile_leaf(ipft) / c2n_leaf(ipft) * cpatch%broot(ico) &
                                + f_labile_stem(ipft) / c2n_leaf(ipft)                     &
                                * ( cpatch%bsapwoodb(ico) + cpatch%bbarkb(ico)             &
                                  + cpatch%bdeadb   (ico) ) )

            csite%stgn_in(ipa) = csite%stgn_in(ipa) + cpatch%nplant(ico)                   &
                               * ( ( 1.0 - f_labile_leaf(ipft) ) * cpatch%bleaf(ico)       &
                                 + ( 1.0 - f_labile_stem(ipft) )                           &
                                 * ( cpatch%bsapwooda(ico) + cpatch%bbarka(ico)            &
                                   + cpatch%bdeada   (ico)                      ) )        &
                                 / c2n_stem(ipft)

            csite%stsn_in(ipa) = csite%stsn_in(ipa) + cpatch%nplant(ico)                   &
                               * ( ( 1.0 - f_labile_leaf(ipft) ) * cpatch%broot(ico)       &
                                 + ( 1.0 - f_labile_stem(ipft) )                           &
                                 * ( cpatch%bsapwoodb(ico) + cpatch%bbarkb(ico)            &
                                   + cpatch%bdeadb   (ico)                      ) )        &
                                 / c2n_stem(ipft)
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!




      !----- Debugging message. -----------------------------------------------------------!
      if (print_fuse_details) then
         write(unit=23,fmt='(a)') '======================================================='
         write(unit=23,fmt='(a)') '======================================================='
         write(unit=23,fmt='(a)') ' '
         write(unit=23,fmt='(a)') ' '
         write(unit=23,fmt='(a)') ' '
         write(unit=23,fmt='(a)') ' '
         write(unit=23,fmt='(a)') ' '
         close(unit=23,status='keep')
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Update total changes in energy due to change in vegetation mass and internal    !
      ! water, and send the intercepted water to the canopy air space.  Total enthalpy due !
      ! to the forced boiling will be consistently updated because the higher specific     !
      ! humidity will translate into more enthalpy.                                        !
      !------------------------------------------------------------------------------------!
      if (checkbudget .and. (.not. is_initial)) then
         csite%ebudget_hcapeffect(ipa) = csite%ebudget_hcapeffect(ipa)                     &
                                       - elim_e_hcap * frqsumi
         csite%ebudget_wcapeffect(ipa) = csite%ebudget_wcapeffect(ipa)                     &
                                       - elim_e_wcap * frqsumi
         csite%wbudget_wcapeffect(ipa) = csite%wbudget_wcapeffect(ipa)                     &
                                       - elim_w_wcap * frqsumi
      end if
      if (veg_boil_tot > 0.0) then
         can_prss           = reducedpress(cmet%prss,cmet%atm_theta,cmet%atm_shv           &
                                          ,cmet%geoht,csite%can_theta(ipa)                 &
                                          ,csite%can_shv(ipa),csite%can_depth(ipa) )
         can_rhos           = idealdenssh (can_prss,csite%can_temp(ipa),csite%can_shv(ipa))
         csite%can_shv(ipa) = csite%can_shv(ipa)                                           &
                            + veg_boil_tot / (csite%can_depth(ipa) * can_rhos)
      end if
      !------------------------------------------------------------------------------------!




      !----- Copy the remaining cohorts to a temporary patch ------------------------------!
      call allocate_patchtype(temppatch,count(remain_table))
      call copy_patchtype_mask(cpatch,temppatch,remain_table,size(remain_table)            &
                              ,count(remain_table))
      !------------------------------------------------------------------------------------!



      !----- Reallocate the new patch and populate with the saved cohorts -----------------!
      call deallocate_patchtype(cpatch)
      call allocate_patchtype(cpatch,count(remain_table))
      call copy_patchtype(temppatch,cpatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
      call sort_cohorts(cpatch)
      !------------------------------------------------------------------------------------!



      !----- Deallocate the temporary patch -----------------------------------------------!
      call deallocate_patchtype(temppatch)
      deallocate(temppatch)
      deallocate(remain_table)
      !------------------------------------------------------------------------------------!



      !----- Update the cohort census at the site level -----------------------------------!
      csite%cohort_count(ipa) = cpatch%ncohorts
      if (elim_lai < 0. .or. elim_nplant < 0.) then
         write (unit=*,fmt='(a,1x,es12.5)') 'TERMINATE: ELIM_LAI=',elim_lai
         write (unit=*,fmt='(a,1x,es12.5)') 'TERMINATE: ELIM_NPLANT=',elim_nplant
      end if
      !------------------------------------------------------------------------------------!


      return
   end subroutine terminate_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: TERMINATE_PATCHES 
   !> \brief This subroutine will eliminate tiny or empty patches. This is intended to
   !> eliminate patches that have little contribution and thus we can speed up the run.
   !---------------------------------------------------------------------------------------!
   subroutine terminate_patches(csite,lai_criterion)

      use ed_state_vars      , only : polygontype        & ! Structure
                                    , sitetype           & ! Structure
                                    , patchtype          ! ! Structure
      use allometry          , only : size2bl            ! ! function
      use fusion_fission_coms, only : pat_laimax_fine    ! ! intent(in)
      use disturb_coms       , only : min_patch_area     ! ! intent(in)
      use pft_coms           , only : SLA                ! ! intent(in)
      use rk4_coms           , only : tiny_offset        ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)       , target      :: csite         ! Current site
      logical              , optional    :: lai_criterion ! Use LAI to decide patch term.
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)       , pointer     :: tempsite      ! Scratch site
      type(patchtype)      , pointer     :: cpatch        ! Current patch
      integer                            :: ipa           ! Patch counter
      integer                            :: ico           ! Cohort counter
      integer                            :: ipft          ! PFT index
      logical, dimension(:), allocatable :: remain_table  ! Flag: this patch will remain.
      real(kind=8)                       :: total_area    ! Area of removed patches
      real(kind=8)                       :: elim_area     ! Area of removed patches
      real(kind=8)                       :: new_area      ! Just to make sure area is 1.
      real                               :: area_scale    ! Scaling area factor.
      real                               :: bleaf_max     ! Maximum leaf biomass
      real                               :: pat_lai_max   ! Maximum patch LAI
      logical                            :: check_lai     ! Local version of lai_criterion
      !----- External functions. ----------------------------------------------------------!
      real                 , external    :: sngloff       ! Safe dble-sngl conversion
      !------------------------------------------------------------------------------------!


      !----- Check whether we should check area (default) or LAI. -------------------------!
      if (present(lai_criterion)) then
         check_lai = lai_criterion
      else
         check_lai = .false.
      end if
      !------------------------------------------------------------------------------------!


      allocate (remain_table(csite%npatches))
      remain_table(:) = .true.

      !------------------------------------------------------------------------------------!
      !     Loop through all the patches in this site and determine which of these patches !
      ! is too small in area to be valid. Remove these patches via the mask function.      !
      ! Realocate a new site with only the valid patches, and normalize their areas and    !
      ! plant densities to reflect the area loss.                                          !
      !------------------------------------------------------------------------------------!
      elim_area  = 0.d0
      total_area = 0.d0
      do ipa = 1,csite%npatches
         if (check_lai) then
            !----- Compute the maximum patch-level LAI (i.e. all cohorts fully flushed) ---!
            cpatch      => csite%patch(ipa)
            pat_lai_max = 0.0
            do ico=1,cpatch%ncohorts
               ipft        = cpatch%pft(ico)
               bleaf_max   = size2bl(cpatch%dbh(ico),cpatch%height(ico),cpatch%sla(ico)    &
                                    ,ipft)
               pat_lai_max = pat_lai_max + cpatch%nplant(ico) * SLA(ipft) * bleaf_max
            end do
            !------------------------------------------------------------------------------!

            !----- In case the patch is unreasonably leafy, get rid of it. ----------------!
            if (pat_lai_max > pat_laimax_fine) then
               elim_area         = elim_area + dble(csite%area(ipa))
               remain_table(ipa) = .false.
            end if
            !------------------------------------------------------------------------------!
         elseif (csite%area(ipa) < min_patch_area) then
            !----- In case this patch has a very small area, get rid of it. ---------------!
            elim_area         = elim_area + dble(csite%area(ipa))
            remain_table(ipa) = .false.
            !------------------------------------------------------------------------------!
         end if
         total_area = total_area + dble(csite%area(ipa))
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
      new_area   = 0.d0
      area_scale = sngloff(1.d0 / (total_area - elim_area),tiny_offset)
      do ipa = 1,csite%npatches
         csite%area(ipa) = csite%area(ipa) * area_scale
         new_area        = new_area + dble(csite%area(ipa))
      end do

      if (abs(new_area-1.d0) > 1.d-5) then
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
   !  SUBROUTINE: RESCALE_PATCHES   
   !> \brief This subroutine will rescale the area of the patches.
   !> This is almost the same as the terminate_patches subroutine,
   !> except that no patch is removed.
   !---------------------------------------------------------------------------------------!
   subroutine rescale_patches(csite)
      use update_derived_utils, only : update_cohort_extensive_props ! ! sub-routine
      use ed_state_vars       , only : polygontype                   & ! Structure
                                     , sitetype                      & ! Structure
                                     , patchtype                     ! ! Structure
      use disturb_coms        , only : min_patch_area                ! ! intent(in)
      use allometry           , only : size2bl                       ! ! function
      use ed_max_dims         , only : n_dist_types                  & ! intent(in)
                                     , n_pft                         ! ! intent(in)
      
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

      !Manfredo: inefficient: (.not. onlyone) should be checked outside the do loop
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
                              * size2bl(cpatch%dbh(ico),cpatch%height(ico)                 &
                                       ,cpatch%sla(ico),ipft)
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
   !  SUBROUTINE: FUSE_COHORTS
   !> \brief This subroutine will perform cohort fusion based on various
   !> similarity criteria to determine whether they can be fused with no 
   !> significant loss of information. 
   !> \details The user is welcome to set up a benchmark, but should be 
   !> aware that no miracles will happen here. If there are more very distinct
   !> cohorts than maxcohort, then the user will need to live with that and
   !> accept life is not always fair with those with limited computational resources.           
   !---------------------------------------------------------------------------------------!
   subroutine new_fuse_cohorts(csite,ipa, lsl, fuse_initial)

      use ed_state_vars       , only : sitetype            & ! Structure
                                     , patchtype           ! ! Structure
      use pft_coms            , only : dbh_crit            & ! intent(in)
                                     , hgt_max             & ! intent(in)
                                     , is_grass            & ! intent(in)
                                     , veg_hcap_min        & ! intent(in)
                                     , qsw                 & ! intent(in)
                                     , qbark               & ! intent(in)
                                     , agf_bs              ! ! intent(in)
      use fusion_fission_coms , only : niter_cohfus        & ! intent(in)
                                     , coh_size_tol_min    & ! intent(in)
                                     , coh_size_tol_mult   & ! intent(in)
                                     , coh_size_tol_max    & ! intent(in)
                                     , lai_tol             ! ! intent(in)
      use ed_max_dims         , only : n_pft               ! ! intent(in)
      use mem_polygons        , only : maxcohort           ! ! intent(in)
      use allometry           , only : size2bl             ! ! function
      use ed_misc_coms        , only : igrass              ! ! intent(in)
      use ed_therm_lib        , only : calc_veg_hcap       ! ! subroutine
      use stable_cohorts      , only : is_resolvable       ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target      :: csite             ! Current site
      integer                , intent(in)  :: ipa               ! Current patch ID
      integer                , intent(in)  :: lsl               ! Lowest soil level
      logical                , intent(in)  :: fuse_initial      ! Initialisation step?
      !----- Local arrays -----------------------------------------------------------------!
      logical, dimension(:)  , allocatable :: fuse_table     ! Flag, remaining cohorts
      type(patchtype)        , pointer     :: cpatch         ! Current patch
      type(patchtype)        , pointer     :: temppatch      ! Scratch patch
      !----- Local scalars. ---------------------------------------------------------------!
      integer      :: donc            ! Index: donor cohort
      integer      :: recc            ! Index: receptor cohort
      logical      :: donc_resolv     ! Donor cohort is resolvable
      logical      :: dr_may_fuse     ! Donor and receptor may be fused.
      logical      :: dr_eqv_recruit  ! Donor and receptor have the same recruit status.
      logical      :: dr_eqv_phen     ! Donor and receptor have the same phenology status.
      logical      :: dr_eqv_small    ! Donor and receptor have the same hydro-size status.
      logical      :: dr_le_lai_max   ! Donor and receptordo not exceed maximum LAI.
      real         :: newn            ! New nplants of merged coh.
      real         :: donc_lai_max    ! Maximum LAI: donor cohort
      real         :: donc_bleaf_max  ! Maximum BLeaf: donor cohort
      real         :: donc_bsapa_max  ! Maximum BSapwood (AG): donor cohort
      real         :: donc_bbarka_max ! Maximum BBark (AG): donor cohort
      real         :: donc_lhcap_max  ! Maximum leaf heat capacity: donor cohort
      real         :: donc_whcap_max  ! Maximum wood heat capacity: donor cohort
      real         :: recc_lai_max    ! Maximum LAI: receptor cohort
      real         :: recc_bleaf_max  ! Maximum BLeaf: receptor cohort
      real         :: recc_bsapa_max  ! Maximum BSapwood (AG): receptor cohort
      real         :: recc_bbarka_max ! Maximum BBark (AG): receptor cohort
      real         :: recc_lhcap_max  ! Maximum leaf heat capacity: receptor cohort
      real         :: recc_whcap_max  ! Maximum wood heat capacity: receptor cohort
      real         :: total_size      ! Total size
      real         :: coh_size_tol    ! Relative size tolerance
      integer      :: ncohorts_old    ! # of coh. before fusion test
      real         :: diff_dbh        ! Absolute DBH difference
      real         :: diff_hgt        ! Absolute height difference
      real         :: new_size        ! New size
      integer      :: ifus            ! Counter: fusion iteractions
      integer      :: dpft            ! PFT of donor cohort
      integer      :: rpft            ! PFT of receptor cohort
      logical      :: any_fusion      ! Flag: was there any fusion?
      !------------------------------------------------------------------------------------!


      cpatch => csite%patch(ipa)

      !------------------------------------------------------------------------------------!
      !     Return if maxcohort is 0 (flag for no cohort fusion), or if the patch is empty !
      ! or has a single cohort.                                                            !
      !------------------------------------------------------------------------------------!
      if (maxcohort == 0 .or. cpatch%ncohorts < 2) return
      !------------------------------------------------------------------------------------!

      !----- Initialize table. In principle, all cohorts stay. ----------------------------!
      allocate(fuse_table(cpatch%ncohorts))
      fuse_table(:) = .true.
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Start with minimum tolerance, iterate and relax tolerance in case it still    !
      ! has too many cohorts.                                                              !
      !------------------------------------------------------------------------------------!
      coh_size_tol = coh_size_tol_min
      coh_fusion: do ifus=1,niter_cohfus
         ncohorts_old =  count(fuse_table) ! Save current number of cohorts ---------------!



         !---------------------------------------------------------------------------------!
         !     Outer loop is the receptor cohort, which is more efficient to collect       !
         ! multiple similar cohorts.                                                       !
         !---------------------------------------------------------------------------------!
         recloop: do recc = 1,cpatch%ncohorts-1
            !------------------------------------------------------------------------------!
            !      Make sure this cohort hasn't been fused yet, otherwise skip it (this    !
            ! should never happen by the way).                                             !
            !------------------------------------------------------------------------------!
            if (.not. fuse_table(recc)) cycle recloop
            !------------------------------------------------------------------------------!

            !----- Skip cohort in case it is inviable.  Terminate it instead. -------------!
            if (.not. cpatch%is_viable(recc)) cycle recloop
            !------------------------------------------------------------------------------!

            !----- Handy aliases. ---------------------------------------------------------!
            rpft          = cpatch%pft(recc)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Inner loop is the donor cohort loop.                                    !
            !------------------------------------------------------------------------------!
            donloop:do donc = recc+1,cpatch%ncohorts
               !---------------------------------------------------------------------------!
               !      Make sure this cohort hasn't been fused yet, otherwise skip it.      !
               !---------------------------------------------------------------------------!
               if (.not. fuse_table(donc)) cycle donloop
               !---------------------------------------------------------------------------!

               !----- Skip cohort in case it is inviable.  Terminate it instead. ----------!
               if (.not. cpatch%is_viable(donc)) cycle donloop
               !---------------------------------------------------------------------------!


               !----- Handy aliases. ------------------------------------------------------!
               dpft        = cpatch%pft(donc)
               donc_resolv = cpatch%leaf_resolvable(donc) .or. cpatch%wood_resolvable(donc)
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !     Initial tests, so we don't compare cohort sizes between cohorts that  !
               ! can never be fused.                                                       !
               !                                                                           !
               ! 1. Cohorts have different PFTs.                                           !
               !                                                                           !
               ! The following checks are only applied in case the donor cohort is         !
               !    resolvable or this is not called during initialisation.                !
               !    2. Cohorts have different recruit statuses.                            !
               !    3. Cohorts have different phenology statuses.                          !
               !                                                                           !
               !---------------------------------------------------------------------------!
               if (dpft /= rpft) cycle donloop
               if (donc_resolv .or. (.not. fuse_initial)) then
                  dr_eqv_recruit =                                                         &
                    (cpatch%first_census    (donc) == cpatch%first_census    (recc)) .and. &
                    (cpatch%new_recruit_flag(donc) == cpatch%new_recruit_flag(recc)) .and. &
                    (cpatch%recruit_dbh     (donc) == cpatch%recruit_dbh     (recc)) .and. &
                    (cpatch%census_status   (donc) == cpatch%census_status   (recc))
                  dr_eqv_phen    =                                                         &
                     cpatch%phenology_status(donc) == cpatch%phenology_status(recc)
                  dr_eqv_small   =                                                         &
                     cpatch%is_small        (donc) .eqv. cpatch%is_small     (recc)
                  if (.not. dr_eqv_recruit) cycle donloop
                  if (.not. dr_eqv_phen   ) cycle donloop
                  if (.not. dr_eqv_small  ) cycle donloop
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Find maximum LAI.  Here we must check life form first.                !
               !---------------------------------------------------------------------------!
               if (is_grass(rpft) .and. igrass == 1) then
                  !----- New grasses.  Use actual LAI. ------------------------------------!
                  donc_lai_max = cpatch%lai(donc)
                  recc_lai_max = cpatch%lai(recc)
                  !------------------------------------------------------------------------!
               else
                  !----- Trees or old grasses. Use on-allometry LAI. ----------------------!
                  donc_lai_max = cpatch%nplant(donc)                                       &
                               * size2bl(cpatch%dbh(donc),cpatch%height(donc)              &
                                        ,cpatch%sla(donc),dpft)                            &
                               * cpatch%sla(donc)
                  recc_lai_max = cpatch%nplant(recc)                                       &
                               * size2bl(cpatch%dbh(recc),cpatch%height(recc)              &
                                        ,cpatch%sla(recc),rpft)                            &
                               * cpatch%sla(recc)
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               ! 4. Combined LAI shall not exceed maximum LAI for any cohort.  This is to  !
               !    prevent bulky cohorts, which prevents self-thinning to work properly.  !
               !    This check is not carried out during initialisation so tiny cohorts    !
               !    can be fused with large ones, and thus improving carbon conservation.  !
               !---------------------------------------------------------------------------!
               if (.not. fuse_initial) then
                  !------------------------------------------------------------------------!
                  !      Prevent fusion in case the cohort would be too leafy.             !
                  !------------------------------------------------------------------------!
                  dr_le_lai_max = ( donc_lai_max + recc_lai_max ) <= lai_tol
                  if (.not. dr_le_lai_max) cycle donloop
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Test for similarity.                                                  !
               !---------------------------------------------------------------------------!
               diff_dbh    = abs ( cpatch%dbh   (donc) - cpatch%dbh   (recc) )
               diff_hgt    = abs ( cpatch%height(donc) - cpatch%height(recc) )
               dr_may_fuse = ( diff_dbh < (dbh_crit(dpft)  * coh_size_tol) ) .and.         &
                             ( diff_hgt < (hgt_max (dpft)  * coh_size_tol) )
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      In case this is the initialisation, we also check whether the donor  !
               ! LAI is so small that it would be turned off.  In case so, we relax the    !
               ! size requirement to the maximum acceptable and fuse them.  This reduces   !
               ! the number of initial cohorts that would be killed otherwise just because !
               ! of small population.  This situation frequently occurs when large plots   !
               ! or airborne lidar are used to initialise the model.                       !
               !---------------------------------------------------------------------------!
               if (fuse_initial .and. (.not. dr_may_fuse)) then
                  !------------------------------------------------------------------------!
                  !    Find potential heat capacity -- Receptor cohort.  This is done      !
                  ! inside the donor loop because the receptor may change when we fuse     !
                  ! cohorts.                                                               !
                  !------------------------------------------------------------------------!
                  recc_bleaf_max  = size2bl(cpatch%dbh(recc),cpatch%height(recc)           &
                                           ,cpatch%sla(recc),rpft)
                  recc_bsapa_max  = agf_bs(rpft)                                           &
                                  * recc_bleaf_max * qsw  (rpft) * cpatch%height(recc)
                  recc_bbarka_max = agf_bs(rpft)                                           &
                                  * recc_bleaf_max * qbark(rpft) * cpatch%height(recc)
                  call calc_veg_hcap(recc_bleaf_max,cpatch%bdeada(recc),recc_bsapa_max     &
                                    ,recc_bbarka_max,cpatch%nplant(recc),rpft              &
                                    ,recc_lhcap_max,recc_whcap_max)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    Find potential heat capacity -- Donor cohort.                       !
                  !------------------------------------------------------------------------!
                  donc_bleaf_max  = size2bl(cpatch%dbh(donc),cpatch%height(donc)           &
                                           ,cpatch%sla(donc),dpft)
                  donc_bsapa_max  = agf_bs(dpft)                                           &
                                  * donc_bleaf_max * qsw  (dpft) * cpatch%height(donc)
                  donc_bbarka_max = agf_bs(dpft)                                           &
                                  * donc_bleaf_max * qbark(dpft) * cpatch%height(donc)
                  call calc_veg_hcap(donc_bleaf_max,cpatch%bdeada(donc),donc_bsapa_max     &
                                    ,donc_bbarka_max,cpatch%nplant(donc),dpft              &
                                    ,donc_lhcap_max,donc_whcap_max)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     In case heat capacity is less than minimum, ignore the size        !
                  ! similarity and fuse the cohort.                                        !
                  !------------------------------------------------------------------------!
                  diff_dbh    = abs ( cpatch%dbh   (donc) - cpatch%dbh   (recc) )
                  diff_hgt    = abs ( cpatch%height(donc) - cpatch%height(recc) )
                  dr_may_fuse = ( ( recc_lhcap_max < veg_hcap_min(rpft) ) .or.             &
                                  ( donc_lhcap_max < veg_hcap_min(dpft) )         ) .and.  &
                                ( diff_dbh < (dbh_crit(dpft)  * coh_size_tol_max) ) .and.  &
                                ( diff_hgt < (hgt_max (dpft)  * coh_size_tol_max) )
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Fuse cohorts in case they are very similar.                          !
               !---------------------------------------------------------------------------!
               if (dr_may_fuse) then

                  !----- New cohort has the total number of plants ------------------------!
                  newn = cpatch%nplant(donc) + cpatch%nplant(recc)
                  !------------------------------------------------------------------------!

                  !----- Check the total size of this cohort before and after fusion. -----!
                  total_size = cpatch%nplant(donc) * ( cpatch%balive(donc)                 &
                                                     + cpatch%bdeada(donc)                 &
                                                     + cpatch%bdeadb(donc)                 &
                                                     + cpatch%bstorage(donc) )             &
                             + cpatch%nplant(recc) * ( cpatch%balive(recc)                 &
                                                     + cpatch%bdeada(recc)                 &
                                                     + cpatch%bdeadb(recc)                 &
                                                     + cpatch%bstorage(recc) )
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !    In case this is not initialisation, we must temporarily make both   !
                  ! cohorts "resolvable".  This is to avoid energy/water leaks when a non- !
                  ! -resolvable cohort is fused with a resolvable cohort.                  !
                  !------------------------------------------------------------------------!
                  if (.not. fuse_initial) then
                     call is_resolvable(csite,ipa,recc,fuse_initial,.true.                 &
                                       ,'new_fuse_cohorts (recc,before)')
                     call is_resolvable(csite,ipa,donc,fuse_initial,.true.                 &
                                       ,'new_fuse_cohorts (donc,before)')
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Proceed with fusion.                                               !
                  !------------------------------------------------------------------------!
                  call fuse_2_cohorts(cpatch,donc,recc,csite%can_prss(ipa)                 &
                                     ,csite%can_shv(ipa),lsl,fuse_initial)
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     In case this is not initialisation, we check whether the final     !
                  ! fused cohort is resolvable.  In case it is not, then we must subtract  !
                  ! the phenology effect that was temporarily added before fusing the      !
                  ! cohorts.  Unlike the "before" calls, here we do not force cohorts to   !
                  ! be resolvable.                                                         !
                  !------------------------------------------------------------------------!
                  if (.not. fuse_initial) then
                     call is_resolvable(csite,ipa,recc,fuse_initial,.false.                &
                                       ,'new_fuse_cohorts (recc,after)')
                  end if
                  !------------------------------------------------------------------------!


                  !----- Flag donating cohort as gone, so it won't be checked again. ------!
                  fuse_table(donc) = .false.
                  !------------------------------------------------------------------------!


                  !----- Check whether total size and LAI are conserved. ------------------!
                  new_size = cpatch%nplant(recc) * ( cpatch%balive(recc)                   &
                                                   + cpatch%bdeada(recc)                   &
                                                   + cpatch%bdeadb(recc)                   &
                                                   + cpatch%bstorage(recc) )
                  if (new_size < 0.99* total_size .or. new_size > 1.01* total_size )       &
                  then
                     write (unit=*,fmt='(a,1x,es14.7)') 'OLD SIZE: ',total_size
                     write (unit=*,fmt='(a,1x,es14.7)') 'NEW SIZE: ',new_size
                     call fatal_error('Cohort fusion didn''t conserve plant size!!!'       &
                                     &,'new_fuse_cohorts','fuse_fiss_utils.f90')
                  end if
                  !------------------------------------------------------------------------!
               else
                  !------------------------------------------------------------------------!
                  !      If we reach this point, it means that these two cohorts are too   !
                  ! dissimilar in size.  Cohorts are sorted by size, so there is no point  !
                  ! looking for other cohorts to be fused with this receptor.  We can move !
                  ! on to the next receptor.                                               !
                  !------------------------------------------------------------------------!
                  cycle recloop
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do donloop
            !------------------------------------------------------------------------------!
         end do recloop
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      In case we met maxcohort goals or in case tolerance exceeded the maximum,  !
         ! interrupt cohort fusion.                                                        !
         !---------------------------------------------------------------------------------!
         if ( count(fuse_table) <= abs(maxcohort)   ) exit coh_fusion
         !---------------------------------------------------------------------------------!

         coh_size_tol = coh_size_tol * coh_size_tol_mult
         ncohorts_old = count(fuse_table)
      end do coh_fusion
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      In case any fusion has happened, rearrange cohorts.                           !
      !------------------------------------------------------------------------------------!
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
         !---------------------------------------------------------------------------------!

         !----- Now I reallocate the current patch with its new reduced size. -------------!
         call deallocate_patchtype(cpatch)  
         call allocate_patchtype(cpatch,count(fuse_table))
         !---------------------------------------------------------------------------------!
 
         !----- Make fuse_table true to all remaining cohorts. ----------------------------!
         fuse_table(:)                 = .false.
         fuse_table(1:cpatch%ncohorts) = .true.
         call copy_patchtype_mask(temppatch,cpatch,fuse_table,size(fuse_table)             &
                                 ,count(fuse_table))
         !---------------------------------------------------------------------------------!

         !----- Discard the scratch patch. ------------------------------------------------!
         call deallocate_patchtype(temppatch)
         deallocate(temppatch)  
         !---------------------------------------------------------------------------------!

         !----- Sort cohorts by size again, and update the cohort census for this patch. --!
         call sort_cohorts(cpatch)
         csite%cohort_count(ipa) = count(fuse_table)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      !----- Deallocate the aux. table ----------------------------------------------------!
      deallocate(fuse_table)
     
      return
   end subroutine new_fuse_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   ! MLO - Cohort fusion may significantly affect the runs, so I am keeping both the old   !
   !       and new routines, but I think we should probably phase out this and stick with  !
   !       the new one.                                                                    !
   !                                                                                       !
   !   This subroutine will perform cohort fusion based on various similarity criteria to  !
   ! determine whether they can be fused with no significant loss of information. The user !
   ! is welcome to set up a benchmark, but should be aware that no miracles will happen    !
   ! here. If there are more very distinct cohorts than maxcohort, then the user will need !
   ! to live with that and accept life is not always fair with those with limited          !
   ! computational resources.                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine old_fuse_cohorts(csite,ipa, lsl, fuse_initial)

      use ed_state_vars       , only : sitetype            & ! Structure
                                     , patchtype           ! ! Structure
      use pft_coms            , only : hgt_max             & ! intent(in)
                                     , is_grass            ! ! intent(in)
      use fusion_fission_coms , only : fusetol_h           & ! intent(in)
                                     , fusetol             & ! intent(in)
                                     , lai_fuse_tol        & ! intent(in)
                                     , fuse_relax          & ! intent(in)
                                     , coh_tolerance_max   ! ! intent(in)
      use ed_max_dims         , only : n_pft               ! ! intent(in)
      use mem_polygons        , only : maxcohort           ! ! intent(in)
      use allometry           , only : size2bl             ! ! function
      use ed_misc_coms        , only : igrass              ! ! intent(in)
      use stable_cohorts      , only : is_resolvable       ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target      :: csite             ! Current site
      integer                , intent(in)  :: ipa               ! Current patch ID
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
      real                                 :: mean_height    ! Mean height        (???)
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
      !    Calculate mean DBH and HEIGHT to help normalise the differences. As of now,     !
      ! height is not being used, but it might be used  in the future if proved            !
      ! advantageous.                                                                      !
      !------------------------------------------------------------------------------------!
      mean_dbh    = 0.0
      mean_height = 0.0
      nshort      = 0
      ntall       = 0
      do ico3 = 1,cpatch%ncohorts
         !---------------------------------------------------------------------------------!
         !    Get fusion height threshold.  Height is a good predictor when plants are     !
         ! growing in height, but it approaches the maximum height DBH becomes the only    !
         ! possible predictor because height saturates.                                    !
         !---------------------------------------------------------------------------------!
         if (cpatch%height(ico3) < (0.95 * hgt_max(cpatch%pft(ico3))) ) then
            mean_height = mean_height + cpatch%height(ico3)
            nshort      = nshort + 1
         else
            mean_dbh  = mean_dbh + cpatch%dbh(ico3)
            ntall     = ntall + 1
         end if
      end do 
      !------------------------------------------------------------------------------------!
      if (ntall  > 0) mean_dbh    = mean_dbh    / real(ntall)
      if (nshort > 0) mean_height = mean_height / real(nshort)

      !----- Initialize table. In principle, all cohorts stay. ----------------------------!
      allocate(fuse_table(cpatch%ncohorts))
      fuse_table(:) = .true.

      force_fusion: do
         
         ncohorts_old =  count(fuse_table) ! Save current number of cohorts ---------------!
         
         donloop:do donc = 1,cpatch%ncohorts-1
            if (.not. fuse_table(donc)) cycle donloop ! This one is gone, move to next.
            
            if (.not. cpatch%is_viable(donc)) cycle donloop ! Inviable cohort, skip it.

            recloop: do recc = donc+1,cpatch%ncohorts
               if (.not. fuse_table(recc)) cycle recloop ! This one is gone, move to next.
                                                         ! Hope it never happens...

               if (.not. cpatch%is_viable(recc)) cycle recloop ! Inviable cohort, skip it.


               !---------------------------------------------------------------------------!
               !     Test for similarity.  Again, we use height to assess similarity only  !
               ! when the cohort is not approaching the maximum height.  If this is the    !
               ! case, then we use DBH to test.                                            !
               !---------------------------------------------------------------------------!
               if (cpatch%height(donc) >= (0.95 * hgt_max(cpatch%pft(donc))) ) then
                  mean_dbh=0.5*(cpatch%dbh(donc)+cpatch%dbh(recc))
                  fusion_test = ( abs(cpatch%dbh(donc) - cpatch%dbh(recc)))/mean_dbh       &
                              < fusetol * tolerance_mult
               elseif (fuse_relax) then
                  fusion_test = ( abs(cpatch%height(donc) - cpatch%height(recc))           &
                                     / (0.5*(cpatch%height(donc) + cpatch%height(recc))) < &
                                fusetol * tolerance_mult)  
               else
                  fusion_test = (abs(cpatch%height(donc) - cpatch%height(recc))  <         &
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
                                * size2bl(cpatch%dbh(recc),cpatch%height(recc)             &
                                         ,cpatch%sla(recc),cpatch%pft(recc))               &
                                + cpatch%nplant(donc)                                      &
                                * size2bl(cpatch%dbh(donc),cpatch%height(donc)             &
                                         ,cpatch%sla(donc),cpatch%pft(donc)))              &
                                * cpatch%sla(recc)
                  end if

                  !----- Checking the total size of this cohort before and after fusion. --!
                  total_size = cpatch%nplant(donc) * ( cpatch%balive  (donc)               &
                                                     + cpatch%bdeada  (donc)               &
                                                     + cpatch%bdeadb  (donc)               &
                                                     + cpatch%bstorage(donc) )             &
                             + cpatch%nplant(recc) * ( cpatch%balive  (recc)               &
                                                     + cpatch%bdeada  (recc)               &
                                                     + cpatch%bdeadb  (recc)               &
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
                  ! 8. Both cohorts must have the same small/large plant size status.      !
                  !------------------------------------------------------------------------!
                  if (    (cpatch%pft(donc)              == cpatch%pft             (recc)) &
                    .and. (lai_max                        < lai_fuse_tol*tolerance_mult  ) &
                    .and. (cpatch%first_census    (donc) == cpatch%first_census    (recc)) &
                    .and. (cpatch%new_recruit_flag(donc) == cpatch%new_recruit_flag(recc)) &
                    .and. (cpatch%recruit_dbh     (donc) == cpatch%recruit_dbh     (recc)) &
                    .and. (cpatch%census_status   (donc) == cpatch%census_status   (recc)) &
                    .and. (cpatch%phenology_status(donc) == cpatch%phenology_status(recc)) &
                    .and. (cpatch%is_small        (donc) .eqv. cpatch%is_small     (recc)) &
                     ) then

                     !---------------------------------------------------------------------!
                     !    In case this is not initialisation, we must temporarily make     !
                     ! both cohorts "resolvable".  This is to avoid energy/water leaks     !
                     ! when a non-resolvable cohort is fused with a resolvable cohort.     !
                     !---------------------------------------------------------------------!
                     if (.not. fuse_initial) then
                        call is_resolvable(csite,ipa,recc,fuse_initial,.true.              &
                                          ,'old_fuse_cohorts (recc,before)')
                        call is_resolvable(csite,ipa,donc,fuse_initial,.true.              &
                                          ,'old_fuse_cohorts (donc,before)')
                     end if
                     !---------------------------------------------------------------------!


                     !----- Proceed with fusion -------------------------------------------!
                     call fuse_2_cohorts(cpatch,donc,recc,csite%can_prss(ipa)              &
                                        ,csite%can_shv(ipa),lsl,fuse_initial)
                     !---------------------------------------------------------------------!



                     !---------------------------------------------------------------------!
                     !     In case this is not initialisation, we check whether the final  !
                     ! fused cohort is resolvable.  In case it is not, then we must        !
                     ! subtract the phenology effect that was temporarily added before     !
                     ! fusing the cohorts.  Unlike the "before" calls, here we do not      !
                     ! force cohorts to be resolvable.                                     !
                     !---------------------------------------------------------------------!
                     if (.not. fuse_initial) then
                        call is_resolvable(csite,ipa,recc,fuse_initial,.false.             &
                                          ,'old_fuse_cohorts (recc,after)')
                     end if
                     !---------------------------------------------------------------------!



                     !----- Flag donating cohort as gone, so it won't be checked again. ---!
                     fuse_table(donc) = .false.
                     
                     !----- Check whether total size and LAI are conserved. ---------------!
                     new_size = cpatch%nplant(recc) * ( cpatch%balive  (recc)              &
                                                      + cpatch%bdeada  (recc)              &
                                                      + cpatch%bdeadb  (recc)              &
                                                      + cpatch%bstorage(recc) )
                     if (new_size < 0.99* total_size .or. new_size > 1.01* total_size )    &
                     then
                        write (unit=*,fmt='(a,1x,es14.7)') 'OLD SIZE: ',total_size
                        write (unit=*,fmt='(a,1x,es14.7)') 'NEW SIZE: ',new_size
                        call fatal_error('Cohort fusion didn''t conserve plant size!!!'    &
                                        &,'old_fuse_cohorts','fuse_fiss_utils.f90')
                     end if
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    Recalculate the means                                            !
                     !---------------------------------------------------------------------!
                     mean_dbh    = 0.0
                     mean_height = 0.0
                     nshort      = 0
                     ntall       = 0
                     recalcloop: do ico3 = 1,cpatch%ncohorts
                        if (.not. fuse_table(ico3)) cycle recalcloop
                        !----- Get fusion height threshold --------------------------------!
                        if (cpatch%height(ico3) < (0.95 * hgt_max(cpatch%pft(ico3))) ) then
                           mean_height = mean_height + cpatch%height(ico3)
                           nshort      = nshort+1
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
   end subroutine old_fuse_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: SPLIT_COHORTS
   !> \brief This subroutine will split two cohorts if its LAI has become too large.
   !> \details This is only necessary when we solve radiation cohort by cohort rather
   !> than layer by layer.
   !---------------------------------------------------------------------------------------!
   subroutine split_cohorts(csite,ipa,green_leaf_factor,is_initial)
      use update_derived_utils , only : update_cohort_extensive_props ! ! sub-routine
      use ed_state_vars        , only : sitetype                      & ! structure
                                      , patchtype                     & ! structure
                                      , copy_patchtype                ! ! sub-routine
      use pft_coms             , only : is_grass                      ! ! intent(in)
      use fusion_fission_coms  , only : lai_tol                       ! ! intent(in)
      use ed_max_dims          , only : n_pft                         ! ! intent(in)
      use allometry            , only : dbh2h                         & ! function
                                      , bd2dbh                        & ! function
                                      , bl2dbh                        & ! function
                                      , bl2h                          & ! function
                                      , size2bl                       ! ! function
      use ed_misc_coms         , only : igrass                        ! ! intent(in)
      use ed_therm_lib         , only : calc_veg_hcap                 & ! function
                                      , update_veg_energy_cweh        ! ! sub-routine
      use stable_cohorts       , only : is_resolvable                 ! ! sub-routine
      use plant_hydro          , only : rwc2tw                        & ! sub-routine
                                      , twi2twe                       ! ! sub-routine
      implicit none
      !----- Constants --------------------------------------------------------------------!
      real                   , parameter   :: epsilon=0.0001     ! Tweak factor...
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target      :: csite              ! Current site
      integer                , intent(in)  :: ipa                ! Patch index
      real, dimension(n_pft) , intent(in)  :: green_leaf_factor  ! Elongation factor.
      logical                , intent(in)  :: is_initial         ! Call from initialisation?
      !----- Local variables --------------------------------------------------------------!
      type(patchtype)        , pointer     :: cpatch             ! Current patch
      type(patchtype)        , pointer     :: temppatch          ! Temporary patch
      logical, dimension(:)  , allocatable :: split_mask         ! Flag: split this cohort
      integer                              :: ico                ! Counter
      integer                              :: inew               ! Counter
      integer                              :: ncohorts_new       ! New # of cohorts
      integer                              :: tobesplit          ! # of cohorts to be split
      integer                              :: ipft               ! PFT type
      real                                 :: bleaf_mp           ! Maximum possible Bleaf
      real                                 :: tai_mp             ! Maximum possible TAI
      real                                 :: old_leaf_hcap      ! Old heat capacity (leaf)
      real                                 :: old_wood_hcap      ! Old heat capacity (wood)
      real                                 :: old_leaf_water     ! Old sfc. water (leaf)
      real                                 :: old_wood_water     ! Old sfc. water (wood)
      real                                 :: old_leaf_water_im2 ! Old int. water (leaf)
      real                                 :: old_wood_water_im2 ! Old int. water (wood)
      real                                 :: old_nplant         ! Old nplant
      real                                 :: new_nplant         ! New nplant
      real                                 :: old_size           ! Old size
      real                                 :: new_size           ! New size
      !------------------------------------------------------------------------------------!


      !----- Set patch. -------------------------------------------------------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    We add the iterative loop in case cohorts have too high LAI.  The routine will  !
      ! repeat the steps until all cohorts are below the maximum potential TAI.            !
      !------------------------------------------------------------------------------------!
      splitloop: do
         !----- Initialize the vector with splitting table --------------------------------!
         allocate(split_mask(cpatch%ncohorts))
         split_mask(:) = .false.
         old_nplant = 0.
         old_size   = 0.
         !---------------------------------------------------------------------------------! 



         !----- Loop through cohorts ------------------------------------------------------!
         do ico = 1,cpatch%ncohorts
            ipft = cpatch%pft(ico)

            !------------------------------------------------------------------------------!
            !     Bleaf_mp and tai_mp are the maximum potential leaf biomass and           !
            ! associated plant area index, given the seasonal constrain                    !
            ! (green_leaf_factor and SLA).                                                 !
            !------------------------------------------------------------------------------!
            bleaf_mp = green_leaf_factor(ipft)                                             &
                     * size2bl(cpatch%dbh(ico),cpatch%height(ico),cpatch%sla(ico),ipft)
            tai_mp   = cpatch%nplant(ico) * bleaf_mp * cpatch%sla(ico) + cpatch%wai(ico)
            !------------------------------------------------------------------------------! 

            !----- If the resulting TAI is too large, split this cohort. ------------------!
            split_mask(ico) = tai_mp > lai_tol
            
            old_nplant = old_nplant + cpatch%nplant(ico)
            old_size   = old_size   + cpatch%nplant(ico) * ( cpatch%balive  (ico)          &
                                                           + cpatch%bdeada  (ico)          &
                                                           + cpatch%bdeadb  (ico)          &
                                                           + cpatch%bstorage(ico) )
            !------------------------------------------------------------------------------! 
         end do
         !---------------------------------------------------------------------------------! 



         !----- Compute the new number of cohorts. ----------------------------------------!
         tobesplit    = count(split_mask)
         ncohorts_new = cpatch%ncohorts + tobesplit
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Check whether any cohort splitting is necessary.                            !
         !---------------------------------------------------------------------------------!
         if (tobesplit > 0) then

            !----- Allocate the temppatch. ------------------------------------------------!
            nullify(temppatch)
            allocate(temppatch)
            call allocate_patchtype(temppatch,cpatch%ncohorts)
            !------------------------------------------------------------------------------!


            !----- Fill the temp space with the current patches. --------------------------!
            call copy_patchtype(cpatch,temppatch,1,cpatch%ncohorts,1,cpatch%ncohorts)
            !------------------------------------------------------------------------------!

            !----- Deallocate the current patch. ------------------------------------------!
            call deallocate_patchtype(cpatch)
            !------------------------------------------------------------------------------!

            !----- Re-allocate the current patch. -----------------------------------------!
            call allocate_patchtype(cpatch,ncohorts_new)
            !------------------------------------------------------------------------------!

            !----- Transfer the temp values back in. --------------------------------------!
            call copy_patchtype(temppatch,cpatch,1,temppatch%ncohorts,1,temppatch%ncohorts)
            !------------------------------------------------------------------------------!

            !----- Remove the temporary patch. --------------------------------------------!
            call deallocate_patchtype(temppatch)
            deallocate(temppatch)
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Go through cohorts tat will be split, then apply minor changes in        !
            ! structural biomass so the cohorts are different in size.                     !
            !------------------------------------------------------------------------------!
            inew = size(split_mask)
            do ico = 1,size(split_mask)

               if (split_mask(ico)) then
                  !------------------------------------------------------------------------!
                  !   Half the densities of the original cohort.  All "extensive" vari-    !
                  ! ables must be rescaled.                                                !
                  !------------------------------------------------------------------------!
                  call update_cohort_extensive_props(cpatch,ico,ico,0.5)
                  !------------------------------------------------------------------------!


                  !----- Apply these values to the new cohort. ----------------------------!
                  inew = inew+1
                  call copy_patchtype(cpatch,cpatch,ico,ico,inew,inew)
                  !------------------------------------------------------------------------!

                  !----- Tweaking bdead, to ensure carbon is conserved. -------------------!
                  if (is_grass(cpatch%pft(ico)) .and. igrass==1) then 
                     !----- Use bleaf for grass. ------------------------------------------!
                     cpatch%bleaf(ico)  = cpatch%bleaf(ico) * (1.-epsilon)
                     cpatch%dbh  (ico)  = bl2dbh(cpatch%bleaf(ico), cpatch%sla(ico)        &
                                                ,cpatch%pft(ico))
                     cpatch%height(ico) = bl2h(cpatch%bleaf(ico), cpatch%sla(ico)          &
                                              ,cpatch%pft(ico))

                     cpatch%bleaf(inew)  = cpatch%bleaf(inew) * (1.+epsilon)
                     cpatch%dbh  (inew)  = bl2dbh(cpatch%bleaf(inew), cpatch%sla(inew)     &
                                                 ,cpatch%pft(inew))
                     cpatch%height(inew) = bl2h  (cpatch%bleaf(inew), cpatch%sla(inew)     &
                                                 ,cpatch%pft(inew))
                     !---------------------------------------------------------------------!
                  else
                     !-- use bdead for trees and old grasses. -----------------------------!
                     cpatch%bdeada(ico)  = cpatch%bdeada(ico) * (1.-epsilon)
                     cpatch%bdeadb(ico)  = cpatch%bdeadb(ico) * (1.-epsilon)
                     cpatch%dbh   (ico)  = bd2dbh(cpatch%pft(ico),cpatch%bdeada(ico)       &
                                                 ,cpatch%bdeadb(ico))
                     cpatch%height(ico)  = dbh2h(cpatch%pft(ico), cpatch%dbh(ico))

                     cpatch%bdeada(inew) = cpatch%bdeada(inew) * (1.+epsilon)
                     cpatch%bdeadb(inew) = cpatch%bdeadb(inew) * (1.+epsilon)
                     cpatch%dbh   (inew) = bd2dbh(cpatch%pft(inew),cpatch%bdeada(inew)     &
                                                 ,cpatch%bdeadb(inew))
                     cpatch%height(inew) = dbh2h(cpatch%pft(inew), cpatch%dbh(inew))
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Biomass has changed, we must modify the heat capacity, internal    !
                  ! water, and account for these changes in the heat capacity effect.      !
                  !------------------------------------------------------------------------!
                  !----- Original cohort. -------------------------------------------------!
                  old_leaf_hcap      = cpatch%leaf_hcap     (ico)
                  old_wood_hcap      = cpatch%wood_hcap     (ico)
                  old_leaf_water     = cpatch%leaf_water    (ico)
                  old_wood_water     = cpatch%wood_water    (ico)
                  old_leaf_water_im2 = cpatch%leaf_water_im2(ico)
                  old_wood_water_im2 = cpatch%wood_water_im2(ico)
                  call calc_veg_hcap(cpatch%bleaf(ico) ,cpatch%bdeada(ico)                 &
                                    ,cpatch%bsapwooda(ico),cpatch%bbarka(ico)              &
                                    ,cpatch%nplant(ico),cpatch%pft(ico)                    &
                                    ,cpatch%leaf_hcap(ico),cpatch%wood_hcap(ico))
                  call rwc2tw(cpatch%leaf_rwc(ico),cpatch%wood_rwc(ico)                    &
                             ,cpatch%bleaf(ico),cpatch%bsapwooda(ico)                      &
                             ,cpatch%bsapwoodb(ico),cpatch%bdeada(ico),cpatch%bdeadb(ico)  &
                             ,cpatch%broot(ico),cpatch%dbh(ico),cpatch%pft(ico)            &
                             ,cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico))
                  call twi2twe(cpatch%leaf_water_int(ico),cpatch%wood_water_int(ico)       &
                              ,cpatch%nplant(ico),cpatch%leaf_water_im2(ico)               &
                              ,cpatch%wood_water_im2(ico))
                  call update_veg_energy_cweh(csite,ipa,ico,old_leaf_hcap,old_wood_hcap    &
                                             ,old_leaf_water,old_wood_water                &
                                             ,old_leaf_water_im2,old_wood_water_im2        &
                                             ,.true.,is_initial)
                  !----- New cohort. ------------------------------------------------------!
                  old_leaf_hcap      = cpatch%leaf_hcap     (inew)
                  old_wood_hcap      = cpatch%wood_hcap     (inew)
                  old_leaf_water     = cpatch%leaf_water    (inew)
                  old_wood_water     = cpatch%wood_water    (inew)
                  old_leaf_water_im2 = cpatch%leaf_water_im2(inew)
                  old_wood_water_im2 = cpatch%wood_water_im2(inew)
                  call calc_veg_hcap(cpatch%bleaf(inew) ,cpatch%bdeada(inew)               &
                                    ,cpatch%bsapwooda(inew),cpatch%bbarka(inew)            &
                                    ,cpatch%nplant(inew),cpatch%pft(inew)                  &
                                    ,cpatch%leaf_hcap(inew),cpatch%wood_hcap(inew))
                  call rwc2tw(cpatch%leaf_rwc(inew),cpatch%wood_rwc(inew)                  &
                             ,cpatch%bleaf(inew),cpatch%bsapwooda(inew)                    &
                             ,cpatch%bsapwoodb(inew),cpatch%bdeada(inew)                   &
                             ,cpatch%bdeadb(inew),cpatch%broot(inew),cpatch%dbh(inew)      &
                             ,cpatch%pft(inew),cpatch%leaf_water_int(inew)                 &
                             ,cpatch%wood_water_int(inew))
                  call twi2twe(cpatch%leaf_water_int(inew),cpatch%wood_water_int(inew)     &
                              ,cpatch%nplant(inew),cpatch%leaf_water_im2(inew)             &
                              ,cpatch%wood_water_im2(inew))
                  call update_veg_energy_cweh(csite,ipa,inew,old_leaf_hcap,old_wood_hcap   &
                                             ,old_leaf_water,old_wood_water                &
                                             ,old_leaf_water_im2,old_wood_water_im2        &
                                             ,.true.,is_initial)
                  !----- Update the stability status. -------------------------------------!
                  call is_resolvable(csite,ipa,ico ,is_initial,.false.                     &
                                    ,'split_cohorts (old)')
                  call is_resolvable(csite,ipa,inew,is_initial,.false.                     &
                                    ,'split_cohorts (new)')
                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!


            !----- After splitting, cohorts may need to be sorted again... ----------------!
            call sort_cohorts(cpatch)
            !------------------------------------------------------------------------------!


            !----- Check whether the total # of plants is conserved... --------------------!
            new_nplant = 0.
            new_size   = 0.
            do ico=1,cpatch%ncohorts
               new_nplant = new_nplant + cpatch%nplant(ico)
               new_size   = new_size   + cpatch%nplant(ico) * ( cpatch%balive  (ico)       &
                                                              + cpatch%bdeada  (ico)       &
                                                              + cpatch%bdeadb  (ico)       &
                                                              + cpatch%bstorage(ico) )
            end do
            if (new_nplant < 0.99 * old_nplant .or. new_nplant > 1.01 * old_nplant .or.    &
                new_size   < 0.99 * old_size   .or. new_size   > 1.01 * old_size) then
               write (unit=*,fmt='(a,1x,es14.7)') 'OLD NPLANT: ',old_nplant
               write (unit=*,fmt='(a,1x,es14.7)') 'NEW NPLANT: ',new_nplant
               write (unit=*,fmt='(a,1x,es14.7)') 'OLD SIZE:   ',old_size
               write (unit=*,fmt='(a,1x,es14.7)') 'NEW SIZE:   ',new_size
               call fatal_error('Cohort splitting didn''t conserve plants!!!'              &
                                           &,'split_cohorts','fuse_fiss_utils.f90')
            end if
            !------------------------------------------------------------------------------!


            !----- Free memory. -----------------------------------------------------------!
            deallocate(split_mask)
            !------------------------------------------------------------------------------!
         else 
            !----- Free memory then exit loop. --------------------------------------------!
            deallocate(split_mask)
            exit splitloop
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do splitloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine split_cohorts
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: FUSE_2_COHORTS
   !> \brief This subroutine will merge two cohorts into 1.
   !> \details The donating cohort (donc) is the one that will be deallocated,
   !> while the receptor cohort (recc) will contain the information from both cohorts.
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_cohorts(cpatch,donc,recc,can_prss,can_shv,lsl,fuse_initial)
      use ed_state_vars      , only : patchtype              ! ! Structure
      use pft_coms           , only : is_grass               ! ! intent(in)
      use therm_lib          , only : uextcm2tl              & ! subroutine
                                    , vpdefil                & ! subroutine
                                    , qslif                  ! ! function
      use allometry          , only : size2krdepth           & ! function
                                    , bd2dbh                 & ! function
                                    , bl2dbh                 & ! function
                                    , bl2h                   & ! function
                                    , dbh2h                  & ! function
                                    , size2xb                & ! function
                                    , ed_balive              & ! function
                                    , distrib_root           ! ! subroutine
      use ed_max_dims        , only : n_mort                 ! ! intent(in)
      use ed_misc_coms       , only : writing_long           & ! intent(in)
                                    , writing_eorq           & ! intent(in)
                                    , writing_dcyc           & ! intent(in)
                                    , ndcycle                & ! intent(in)
                                    , igrass                 ! ! intent(in)
      use consts_coms        , only : t3ple                  & ! intent(in)
                                    , lnexp_min              & ! intent(in)
                                    , lnexp_max              & ! intent(in)
                                    , tiny_num               ! ! intent(in)
      use fusion_fission_coms, only : corr_cohort            ! ! intent(in)
      use grid_coms          , only : nzg                    ! ! intent(in)
      use plant_hydro        , only : rwc2psi                & ! subroutine
                                    , twe2twi                & ! subroutine
                                    , tw2rwc                 & ! subroutine
                                    , tw2psi                 ! ! subroutine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype) , target     :: cpatch            ! Current patch
      integer                      :: donc              ! Donating cohort.
      integer                      :: recc              ! Receptor cohort.
      real            , intent(in) :: can_prss          ! Canopy air pressure
      real            , intent(in) :: can_shv           ! Canopy air specific humidity
      integer         , intent(in) :: lsl               ! Lowest soil level
      logical         , intent(in) :: fuse_initial      ! Called from initialisation
      !----- Local variables --------------------------------------------------------------!
      integer                      :: imon              ! Month for cb loop
      integer                      :: t                 ! Time of day for dcycle loop
      integer                      :: imty              ! Mortality type
      integer                      :: isl               ! Soil layer
      real                         :: exp_mort_donc     ! exp(mortality) donor
      real                         :: exp_mort_recc     ! exp(mortality) receptor
      real                         :: m_exp_dlnndt_donc ! -exp(dlnndt) donor
      real                         :: m_exp_dlnndt_recc ! -exp(dlnndt) receptor
      real                         :: recc_basarea      ! BA of receiver
      real                         :: donc_basarea      ! BA of donor
      real                         :: recc_bleaf        ! Leaf biomass of receiver
      real                         :: donc_bleaf        ! Leaf biomass of donor
      real                         :: rlai              ! LAI weight of receiver
      real                         :: dlai              ! LAI weight of donor
      real                         :: rwai              ! WAI weight of receiver
      real                         :: dwai              ! WAI weight of donor
      real                         :: rnplant           ! nplant weight of receiver
      real                         :: dnplant           ! nplant weight of donor
      real                         :: rba               ! BA weight of receiver
      real                         :: dba               ! BA weight of donor
      real                         :: rbleaf            ! Leaf biomass weight of receiver
      real                         :: dbleaf            ! Leaf biomass weight of donor
      real                         :: total_transp      ! Transpiration to be conserved
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Find the scaling factor for variables that are not "extensive".                 !
      !  - If the unit is X/plant, then we scale by nplant.                                !
      !  - If the unit is X/m2_leaf, then we scale by LAI.                                 !
      !  - If the unit is X/m2_wood, then we scale by WAI.                                 !
      !  - If the unit is X and related to basal area, then we scale by BA.                !
      !  - If the unit is X/m2_gnd, then we add, since they are "extensive".               !
      !------------------------------------------------------------------------------------!
      rnplant = cpatch%nplant(recc) / (cpatch%nplant(recc) + cpatch%nplant(donc))
      dnplant = 1.0 - rnplant
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    This is a fix for when two cohorts with very very low LAI are fused             !
      ! it prevents numerical errors (RGK 8-18-2014).                                      !
      ! (MLO 11-24-2014): turned rlai and dlai relative weights, so it works in all cases. !
      ! Also, applied the same idea to WAI-dependent variables.                            !
      !------------------------------------------------------------------------------------!
      if ((cpatch%lai(recc) + cpatch%lai(donc)) > 0. ) then
         rlai    = cpatch%lai(recc) / ( cpatch%lai(recc) + cpatch%lai(donc) )
         dlai    = 1.0 - rlai
      else
         rlai    = 0.5
         dlai    = 0.5
      end if
      if ((cpatch%wai(recc) + cpatch%wai(donc)) > 0. ) then
         rwai    = cpatch%wai(recc) / ( cpatch%wai(recc) + cpatch%wai(donc) )
         dwai    = 1.0 - rwai
      else
         rwai    = 0.5
         dwai    = 0.5
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Scaling factor for basal area.                                                  !
      !------------------------------------------------------------------------------------!
      recc_basarea = cpatch%nplant(recc) * cpatch%basarea(recc)
      donc_basarea = cpatch%nplant(donc) * cpatch%basarea(donc)
      if ((recc_basarea + donc_basarea) > 0. ) then
         rba = recc_basarea / ( recc_basarea + donc_basarea )
         dba = 1.0 - rba
      else
         rba = 0.5
         dba = 0.5
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Scaling factor for leaf biomass.                                                !
      !------------------------------------------------------------------------------------!
      recc_bleaf = cpatch%nplant(recc) * cpatch%bleaf(recc)
      donc_bleaf = cpatch%nplant(donc) * cpatch%bleaf(donc)
      if ((recc_bleaf + donc_bleaf) > 0. ) then
         rbleaf = recc_bleaf / ( recc_bleaf + donc_bleaf )
         dbleaf = 1.0 - rbleaf
      else
         rbleaf = 0.5
         dbleaf = 0.5
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Update nplant, LAI, WAI, and crown area.   Many variables use them as weights, !
      ! but the weights are set above -- DO NOT USE CPATCH%NPLANT, CPATCH%LAI, CPATCH%WAI  !
      ! AND CPATCH%CROWN_AREA DIRECTLY AS WEIGHTS.  Bad things will happen.                !
      !------------------------------------------------------------------------------------!
      cpatch%nplant     (recc) = cpatch%nplant(recc) + cpatch%nplant(donc)
      cpatch%lai        (recc) = cpatch%lai   (recc) + cpatch%lai   (donc)
      cpatch%wai        (recc) = cpatch%wai   (recc) + cpatch%wai   (donc)
      !----- Make sure that crown area is bounded. ----------------------------------------!
      cpatch%crown_area (recc) = min(1.,cpatch%crown_area(recc)  + cpatch%crown_area(donc))
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     Fuse all carbon pools.  This is done before we find DBH and height, as these   !
      ! quantities depend on whether the PFT is (new) grass or not.                        !
      !------------------------------------------------------------------------------------!
      cpatch%bdeada   (recc) = cpatch%bdeada   (recc) * rnplant                            &
                             + cpatch%bdeada   (donc) * dnplant
      cpatch%bdeadb   (recc) = cpatch%bdeadb   (recc) * rnplant                            &
                             + cpatch%bdeadb   (donc) * dnplant
      cpatch%bleaf    (recc) = cpatch%bleaf    (recc) * rnplant                            &
                             + cpatch%bleaf    (donc) * dnplant
      cpatch%broot    (recc) = cpatch%broot    (recc) * rnplant                            &
                             + cpatch%broot    (donc) * dnplant
      cpatch%bsapwooda(recc) = cpatch%bsapwooda(recc) * rnplant                            &
                             + cpatch%bsapwooda(donc) * dnplant
      cpatch%bsapwoodb(recc) = cpatch%bsapwoodb(recc) * rnplant                            &
                             + cpatch%bsapwoodb(donc) * dnplant
      cpatch%bbarka   (recc) = cpatch%bbarka   (recc) * rnplant                            &
                             + cpatch%bbarka   (donc) * dnplant
      cpatch%bbarkb   (recc) = cpatch%bbarkb   (recc) * rnplant                            &
                             + cpatch%bbarkb   (donc) * dnplant
      cpatch%bstorage (recc) = cpatch%bstorage (recc) * rnplant                            &
                             + cpatch%bstorage (donc) * dnplant
      cpatch%btimber  (recc) = cpatch%btimber  (recc) * rnplant                            &
                             + cpatch%btimber  (donc) * dnplant
      cpatch%bseeds   (recc) = cpatch%bseeds   (recc) * rnplant                            &
                             + cpatch%bseeds   (donc) * dnplant
      cpatch%byield   (recc) = cpatch%byield   (recc) * rnplant                            &
                             + cpatch%byield   (donc) * dnplant
      !------------------------------------------------------------------------------------!



      !------ Find balive from the other pools. -------------------------------------------!
      cpatch%balive(recc) = ed_balive(cpatch,recc)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Find DBH and height.  Make sure that carbon is conserved.                     !
      !------------------------------------------------------------------------------------!
      if (is_grass(cpatch%pft(recc)) .and. igrass == 1) then
          !----- New grass scheme, use bleaf then find DBH and height. --------------------!
          cpatch%dbh   (recc) = bl2dbh(cpatch%bleaf(recc),cpatch%sla(recc)                 &
                                      ,cpatch%pft(recc))
          cpatch%height(recc) = bl2h  (cpatch%bleaf(recc),cpatch%sla(recc)                 &
                                      ,cpatch%pft(recc))
          !--------------------------------------------------------------------------------!
      else
          !----- Trees, or old grass scheme.  Use bdead then find DBH and height. ---------!
          cpatch%dbh   (recc) = bd2dbh(cpatch%pft(recc),cpatch%bdeada(recc)                &
                                                       ,cpatch%bdeadb(recc))
          cpatch%height(recc) = dbh2h(cpatch%pft(recc),cpatch%dbh(recc))
          !--------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !----- Rooting depth. ---------------------------------------------------------------!
      cpatch%krdepth(recc) = size2krdepth(cpatch%height(recc),cpatch%dbh(recc)             &
                                         ,cpatch%pft(recc),lsl)
      !------------------------------------------------------------------------------------!



      !----- Update the vertical distribution of roots. -----------------------------------!
      call distrib_root(cpatch%krdepth(recc),cpatch%pft(recc),cpatch%root_frac(:,recc))
      !------------------------------------------------------------------------------------!




      !----- Maintenance costs. -----------------------------------------------------------!
      cpatch%leaf_maintenance  (recc) = cpatch%leaf_maintenance (recc) * rnplant           &
                                      + cpatch%leaf_maintenance (donc) * dnplant
      cpatch%root_maintenance  (recc) = cpatch%root_maintenance (recc) * rnplant           &
                                      + cpatch%root_maintenance (donc) * dnplant
      cpatch%barka_maintenance (recc) = cpatch%barka_maintenance(recc) * rnplant           &
                                      + cpatch%barka_maintenance(donc) * dnplant
      cpatch%barkb_maintenance (recc) = cpatch%barkb_maintenance(recc) * rnplant           &
                                      + cpatch%barkb_maintenance(donc) * dnplant
      cpatch%leaf_drop         (recc) = cpatch%leaf_drop        (recc) * rnplant           &
                                      + cpatch%leaf_drop        (donc) * dnplant
      cpatch%root_drop         (recc) = cpatch%root_drop        (recc) * rnplant           &
                                      + cpatch%root_drop        (donc) * dnplant
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Bark thickness is calculated based on the fused size and biomass.              !
      !------------------------------------------------------------------------------------!
      cpatch%thbark(recc) = size2xb(cpatch%dbh(recc),cpatch%height(recc)                   &
                                   ,cpatch%bbarka(recc),cpatch%bbarkb(recc)                &
                                   ,cpatch%sla(recc),cpatch%pft(recc))
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
      !     CB and CB_Xmax are scaled by population, as they are in kgC/plant.             !
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
      !     Growth and PLC history should be scaled by biomass. Here, we use basal area    !
      !------------------------------------------------------------------------------------!
      do imon = 1,13
         cpatch%ddbh_monthly(imon,recc) = cpatch%ddbh_monthly (imon,recc) * rba            &
                                        + cpatch%ddbh_monthly (imon,donc) * dba
         cpatch%plc_monthly(imon,recc)  = cpatch%plc_monthly  (imon,recc) * rba            &
                                        + cpatch%plc_monthly  (imon,donc) * dba
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
                                  
      cpatch%today_nppbark      (recc) = cpatch%today_nppbark      (recc)                  &
                                       + cpatch%today_nppbark      (donc)
                                  
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

      cpatch%today_stem_resp    (recc) = cpatch%today_stem_resp    (recc)                  &
                                       + cpatch%today_stem_resp    (donc)
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
      cpatch%barka_growth_resp  (recc) = cpatch%barka_growth_resp  (recc) * rnplant        &
                                       + cpatch%barka_growth_resp  (donc) * dnplant
      cpatch%barkb_growth_resp  (recc) = cpatch%barkb_growth_resp  (recc) * rnplant        &
                                       + cpatch%barkb_growth_resp  (donc) * dnplant
      cpatch%leaf_storage_resp  (recc) = cpatch%leaf_storage_resp  (recc) * rnplant        &
                                       + cpatch%leaf_storage_resp  (donc) * dnplant
      cpatch%root_storage_resp  (recc) = cpatch%root_storage_resp  (recc) * rnplant        &
                                       + cpatch%root_storage_resp  (donc) * dnplant
      cpatch%sapa_storage_resp  (recc) = cpatch%sapa_storage_resp  (recc) * rnplant        &
                                       + cpatch%sapa_storage_resp  (donc) * dnplant
      cpatch%sapb_storage_resp  (recc) = cpatch%sapb_storage_resp  (recc) * rnplant        &
                                       + cpatch%sapb_storage_resp  (donc) * dnplant
      cpatch%barka_storage_resp (recc) = cpatch%barka_storage_resp (recc) * rnplant        &
                                       + cpatch%barka_storage_resp (donc) * dnplant
      cpatch%barkb_storage_resp (recc) = cpatch%barkb_storage_resp (recc) * rnplant        &
                                       + cpatch%barkb_storage_resp (donc) * dnplant
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Save total transpiration before updating psi_open, psi_closed, fs_open and     !
      ! fs_closed, and merge these variables in a way that conserves total transpiration.  !
      ! This requirement is essential for consistent plant hydraulics.                     !
      !------------------------------------------------------------------------------------!
      total_transp              = ( cpatch%psi_open(recc) * cpatch%fs_open(recc)           &
                                  + cpatch%psi_closed(recc) * (1. - cpatch%fs_open(recc))) &
                                  * rlai                                                   &
                                + ( cpatch%psi_open(donc) * cpatch%fs_open(donc)           &
                                  + cpatch%psi_closed(donc) * (1. - cpatch%fs_open(donc))) &
                                  * dlai
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
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !     XXT: The original scaling for fs_open can be problematic when photorespiration !
      ! is too high (see photosyn_driv.f90 for details).  In this case, fs_open is         !
      ! decoupled from fsw * fsn.  The current approach ensures that the transpiration of  !
      ! the combined cohort is preserved.                                                  !
      !------------------------------------------------------------------------------------!
      if ( abs(cpatch%psi_open(recc) - cpatch%psi_closed(recc)) < tiny_num) then
         !---------------------------------------------------------------------------------!
         !    Difference between psi_open and psi_closed is too small, likely because the  !
         ! stomata are closed.  Use the original scaling approach.                         !
         !---------------------------------------------------------------------------------!
         cpatch%fs_open(recc)  = cpatch%fsw(recc) * cpatch%fsn(recc)
         !---------------------------------------------------------------------------------!
      else
         !---------------------------------------------------------------------------------!
         !      Conserve the total leaf-level transpiration.                               !
         !---------------------------------------------------------------------------------!
         cpatch%fs_open(recc) = max(0.,min(1.,(total_transp - cpatch%psi_closed(recc))     &
                                   / (cpatch%psi_open(recc) - cpatch%psi_closed(recc))))
         !---------------------------------------------------------------------------------!
      end if
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
      !    Fuse population change rates (y).  By definition, this is the negative of       !
      ! mortality, and y is defined as:                                                    !
      !                                                                                    !
      !                                    ln(N) - ln(A)                                   !
      !                               y = ---------------                                  !
      !                                          dt                                        !
      !                                                                                    !
      ! where N is the current population, and A is the population before the rate was     !
      ! applied.  The cohorts represent a group of individuals with the same size and PFT, !
      ! so they don't mix new recruits and old plants. Therefore, we can assume that N is  !
      ! actually nplant.  We don't know A, but if the change rate is assumed constant      !
      ! during the interval dt, A = N * exp(-y dt).                                        !
      !                                                                                    !
      ! For fusion we don't really care about dt, so any number will do as long as it is   !
      ! the same for both cohorts.  With these assumptions, the change rate for the        !
      ! fused cohort yf is:                                                                !
      !                                                                                    !
      !  yf   =  ln (Nd+Nr) - ln(Ad+Ar) = ln[Nd+Nr] - ln[Nd*exp(-yd) + Nr*exp(-yr)]        !
      !                                                                                    !
      !               / Nd*exp(-yd) + Nr*exp(-yr) \                                        !
      !  yf   =  - ln |---------------------------|                                        !
      !               \         Nd + Nr           /                                        !
      !------------------------------------------------------------------------------------!
      m_exp_dlnndt_donc = exp(max(lnexp_min,min(lnexp_max,-cpatch%monthly_dlnndt(donc))))
      m_exp_dlnndt_recc = exp(max(lnexp_min,min(lnexp_max,-cpatch%monthly_dlnndt(recc))))

      cpatch%monthly_dlnndt(recc) = - log( rnplant * m_exp_dlnndt_recc                     &
                                         + dnplant * m_exp_dlnndt_donc )
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
      cpatch%stem_respiration(recc) = cpatch%stem_respiration(recc)                        &
                                    + cpatch%stem_respiration(donc)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Potential available water and elongation factor can be consider "intensive"    !
      ! variable water.                                                                    !
      !------------------------------------------------------------------------------------!
      cpatch%paw_avg(recc) = cpatch%paw_avg(recc) * rnplant + cpatch%paw_avg(donc) * dnplant
      cpatch%elongf (recc) = cpatch%elongf(recc)  * rnplant + cpatch%elongf (donc) * dnplant
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Light-phenology characteristics.  To conserve maintenance costs, leaf longevity !
      !  must be scaled with leaf biomass (but note that the scaling must be applied to    !
      !  turnover instead of leaf lifespan).  Likewise, SLA is in m2_leaf/kgC, so it must  !
      ! be scaled with leaf biomass.  Carboxylation rate (vm_bar) and dark respiration rate (rd_bar) are in umol/m2_leaf/s, so !
      ! They must be scaled with LAI.                                                        !
      !                                                                                    !
      !------------------------------------------------------------------------------------!
      cpatch%sla         (recc) = cpatch%sla         (recc) * rbleaf                       &
                                + cpatch%sla         (donc) * dbleaf
      cpatch%vm_bar      (recc) = cpatch%vm_bar      (recc) * rlai                         &
                                + cpatch%vm_bar      (donc) * dlai
      cpatch%rd_bar      (recc) = cpatch%rd_bar      (recc) * rlai                         &
                                + cpatch%rd_bar      (donc) * dlai
      !------ For Life span, we must check whether they are non-zero. ---------------------!
      if ( abs(cpatch%llspan(recc)*cpatch%llspan(donc)) > tiny_num ) then
         !---------------------------------------------------------------------------------!
         !     The denominator weights are not inadvertenly swapped.  This happens because !
         ! we are linearly scaling turnover, not leaf life span, to ensure that            !
         ! maintenance costs are consistent.                                               !
         !---------------------------------------------------------------------------------!
         cpatch%llspan   (recc) =   cpatch%llspan    (recc) * cpatch%llspan    (donc)      &
                                / ( cpatch%llspan    (recc) * dbleaf                       &
                                  + cpatch%llspan    (donc) * rbleaf )
         !---------------------------------------------------------------------------------!
      else
         !------ This should only happen when both are zero, so it really doesn't matter. -!
         cpatch%llspan   (recc) = cpatch%llspan      (recc) * rbleaf                       &
                                + cpatch%llspan      (donc) * dbleaf
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Plant hydrodynamics characteristics (XXT).  Extensive internal water (kg/m2)    !
      ! is updated here, and Intensive internal water content (kg/plant) is updated later. !
      ! Water fluxes are weighted by nplant (int) or added (im2).                          !
      !------------------------------------------------------------------------------------!
      cpatch%leaf_water_im2(recc) = cpatch%leaf_water_im2(recc)                            &
                                  + cpatch%leaf_water_im2(donc)
      cpatch%wood_water_im2(recc) = cpatch%wood_water_im2(recc)                            &
                                  + cpatch%wood_water_im2(donc)
      cpatch%wflux_gw      (recc) = cpatch%wflux_gw     (recc) * rnplant                   &
                                  + cpatch%wflux_gw     (donc) * dnplant
      cpatch%wflux_wl      (recc) = cpatch%wflux_wl     (recc) * rnplant                   &
                                  + cpatch%wflux_wl     (donc) * dnplant
      do isl = 1,nzg
         cpatch%wflux_gw_layer(isl,recc) = cpatch%wflux_gw_layer(isl,recc) * rnplant       &
                                         + cpatch%wflux_gw_layer(isl,donc) * dnplant
      enddo
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      We update temperature and liquid water fraction.  We check whether the heat   !
      ! capacity is non-zero.  If it is a normal number, use the standard thermodynamic    !
      ! library, otherwise, average temperature, this is probably a blend of tiny cohorts  !
      ! that couldn't be solved, or the wood is not solved.                                !
      !------------------------------------------------------------------------------------!
      if ( cpatch%leaf_hcap(recc) > 0. ) then
         !----- Update temperature using the standard thermodynamics. ---------------------!
         call uextcm2tl( cpatch%leaf_energy   (recc)                                       &
                       , cpatch%leaf_water    (recc)                                       &
                       + cpatch%leaf_water_im2(recc)                                       &
                       , cpatch%leaf_hcap     (recc)                                       &
                       , cpatch%leaf_temp     (recc)                                       &
                       , cpatch%leaf_fliq     (recc) )
         !---------------------------------------------------------------------------------!
      else 
         !----- Leaf temperature cannot be found using uextcm2tl, this is a singularity. --!
         cpatch%leaf_temp(recc)  = cpatch%leaf_temp(recc) * rnplant                        &
                                 + cpatch%leaf_temp(donc) * dnplant
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Always make liquid fraction consistent with temperature.                    !
         !---------------------------------------------------------------------------------!
         if (cpatch%leaf_temp(recc) == t3ple) then
            cpatch%leaf_fliq(recc) = 0.5
         elseif (cpatch%leaf_temp(recc) > t3ple) then
            cpatch%leaf_fliq(recc) = 1.0
         else
            cpatch%leaf_fliq(recc) = 0.0
         end if
         !---------------------------------------------------------------------------------!
      end if


      if ( cpatch%wood_hcap(recc) > 0. .or. cpatch%wood_water_im2(recc) > 0.) then
         !----- Update temperature using the standard thermodynamics. ---------------------!
         call uextcm2tl( cpatch%wood_energy   (recc)                                       &
                       , cpatch%wood_water    (recc)                                       &
                       + cpatch%wood_water_im2(recc)                                       &
                       , cpatch%wood_hcap     (recc)                                       &
                       , cpatch%wood_temp     (recc)                                       &
                       , cpatch%wood_fliq     (recc) )
         !---------------------------------------------------------------------------------!
      else 
         !----- Wood temperature cannot be found using uextcm2tl, this is a singularity. --!
         cpatch%wood_temp(recc)  = cpatch%wood_temp(recc) * rnplant                        &
                                 + cpatch%wood_temp(donc) * dnplant
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Always make liquid fraction consistent with temperature.                    !
         !---------------------------------------------------------------------------------!
         if (cpatch%wood_temp(recc) == t3ple) then
            cpatch%wood_fliq(recc) = 0.5
         elseif (cpatch%wood_temp(recc) > t3ple) then
            cpatch%wood_fliq(recc) = 1.0
         else
            cpatch%wood_fliq(recc) = 0.0
         end if
         !---------------------------------------------------------------------------------!
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
      !    Recalculate rwc and psi from water_int. This ensures that psi, rwc, and total   !
      ! water are consistent with each other.                                              !
      !------------------------------------------------------------------------------------!
      call twe2twi(cpatch%leaf_water_im2(recc),cpatch%wood_water_im2(recc)                 &
                  ,cpatch%nplant(recc),cpatch%leaf_water_int(recc)                         &
                  ,cpatch%wood_water_int(recc))
      call tw2rwc(cpatch%leaf_water_int(recc),cpatch%wood_water_int(recc)                  &
                 ,cpatch%is_small(recc),cpatch%bleaf(recc),cpatch%bsapwooda(recc)          &
                 ,cpatch%bsapwoodb(recc),cpatch%bdeada(recc),cpatch%bdeadb(recc)           &
                 ,cpatch%broot(recc),cpatch%dbh(recc),cpatch%pft(recc)                     &
                 ,cpatch%leaf_rwc(recc),cpatch%wood_rwc(recc))
      call rwc2psi(cpatch%leaf_rwc(recc),cpatch%wood_rwc(recc),cpatch%pft(recc)            &
                  ,cpatch%leaf_psi(recc),cpatch%wood_psi(recc))
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Phenology/Stomatal variabels associated with plant hydraulics.  These are also  !
      ! scaled by nplant.                                                                  !
      !------------------------------------------------------------------------------------!
      cpatch%high_leaf_psi_days(recc) = nint(                                              &
                                        cpatch%high_leaf_psi_days(recc) * rnplant          &
                                      + cpatch%high_leaf_psi_days(donc) * dnplant)
      cpatch%low_leaf_psi_days(recc)  = nint(                                              &
                                        cpatch%low_leaf_psi_days(recc) * rnplant           &
                                      + cpatch%low_leaf_psi_days(donc) * dnplant)
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
         cpatch%fmean_gpp               (recc) = cpatch%fmean_gpp               (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_gpp               (donc)     &
                                               * dnplant
         cpatch%fmean_npp               (recc) = cpatch%fmean_npp               (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_npp               (donc)     &
                                               * dnplant
         cpatch%fmean_leaf_resp         (recc) = cpatch%fmean_leaf_resp         (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_leaf_resp         (donc)     &
                                               * dnplant
         cpatch%fmean_root_resp         (recc) = cpatch%fmean_root_resp         (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_root_resp         (donc)     &
                                               * dnplant
         cpatch%fmean_stem_resp         (recc) = cpatch%fmean_stem_resp         (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_stem_resp         (donc)     &
                                               * dnplant
         cpatch%fmean_leaf_growth_resp  (recc) = cpatch%fmean_leaf_growth_resp  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_leaf_growth_resp  (donc)     &
                                               * dnplant
         cpatch%fmean_root_growth_resp  (recc) = cpatch%fmean_root_growth_resp  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_root_growth_resp  (donc)     &
                                               * dnplant
         cpatch%fmean_sapa_growth_resp  (recc) = cpatch%fmean_sapa_growth_resp  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_sapa_growth_resp  (donc)     &
                                               * dnplant
         cpatch%fmean_sapb_growth_resp  (recc) = cpatch%fmean_sapb_growth_resp  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_sapb_growth_resp  (donc)     &
                                               * dnplant
         cpatch%fmean_barka_growth_resp (recc) = cpatch%fmean_barka_growth_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_barka_growth_resp (donc)     &
                                               * dnplant
         cpatch%fmean_barkb_growth_resp (recc) = cpatch%fmean_barkb_growth_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_barkb_growth_resp (donc)     &
                                               * dnplant
         cpatch%fmean_leaf_storage_resp (recc) = cpatch%fmean_leaf_storage_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_leaf_storage_resp (donc)     &
                                               * dnplant
         cpatch%fmean_root_storage_resp (recc) = cpatch%fmean_root_storage_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_root_storage_resp (donc)     &
                                               * dnplant
         cpatch%fmean_sapa_storage_resp (recc) = cpatch%fmean_sapa_storage_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_sapa_storage_resp (donc)     &
                                               * dnplant
         cpatch%fmean_sapb_storage_resp (recc) = cpatch%fmean_sapb_storage_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_sapb_storage_resp (donc)     &
                                               * dnplant
         cpatch%fmean_barka_storage_resp(recc) = cpatch%fmean_barka_storage_resp(recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_barka_storage_resp(donc)     &
                                               * dnplant
         cpatch%fmean_barkb_storage_resp(recc) = cpatch%fmean_barkb_storage_resp(recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_barkb_storage_resp(donc)     &
                                               * dnplant
         cpatch%fmean_plresp            (recc) = cpatch%fmean_plresp            (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_plresp            (donc)     &
                                               * dnplant
         cpatch%fmean_light_level       (recc) = cpatch%fmean_light_level       (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_light_level       (donc)     &
                                               * dnplant
         cpatch%fmean_light_level_beam  (recc) = cpatch%fmean_light_level_beam  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_light_level_beam  (donc)     &
                                               * dnplant
         cpatch%fmean_light_level_diff  (recc) = cpatch%fmean_light_level_diff  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_light_level_diff  (donc)     &
                                               * dnplant
         cpatch%fmean_par_level_beam    (recc) = cpatch%fmean_par_level_beam    (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_par_level_beam    (donc)     &
                                               * dnplant
         cpatch%fmean_par_level_diffd   (recc) = cpatch%fmean_par_level_diffd   (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_par_level_diffd   (donc)     &
                                               * dnplant
         cpatch%fmean_par_level_diffu   (recc) = cpatch%fmean_par_level_diffu   (recc)     &
                                               * rnplant                                   &
                                               + cpatch%fmean_par_level_diffu   (donc)     &
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
         !    Plant hydrodynamics characteristics (XXT).  Internal water content and water !
         ! fluxes are weighted by nplant.  Area-based internal water is added.             !
         !---------------------------------------------------------------------------------!
         cpatch%fmean_leaf_water_int(recc) = cpatch%fmean_leaf_water_int(recc) * rnplant   &
                                           + cpatch%fmean_leaf_water_int(donc) * dnplant
         cpatch%fmean_wood_water_int(recc) = cpatch%fmean_wood_water_int(recc) * rnplant   &
                                           + cpatch%fmean_wood_water_int(donc) * dnplant
         cpatch%fmean_leaf_water_im2(recc) = cpatch%fmean_leaf_water_im2(recc)             &
                                           + cpatch%fmean_leaf_water_im2(donc)
         cpatch%fmean_wood_water_im2(recc) = cpatch%fmean_wood_water_im2(recc)             &
                                           + cpatch%fmean_wood_water_im2(donc)
         cpatch%fmean_wflux_wl      (recc) = cpatch%fmean_wflux_wl      (recc) * rnplant   &
                                           + cpatch%fmean_wflux_wl      (donc) * dnplant
         cpatch%fmean_wflux_gw      (recc) = cpatch%fmean_wflux_gw      (recc) * rnplant   &
                                           + cpatch%fmean_wflux_gw      (donc) * dnplant
         do isl = 1,nzg
            cpatch%fmean_wflux_gw_layer(isl,recc) =                                        &
                cpatch%fmean_wflux_gw_layer(isl,recc) * rnplant                            &
              + cpatch%fmean_wflux_gw_layer(isl,donc) * dnplant
         end do
         !---------------------------------------------------------------------------------!


         !----- For daily maximum and minimum psi, we simply use NPLANT as weight. --------!
         cpatch%dmax_leaf_psi(recc) = cpatch%dmax_leaf_psi(recc) * rnplant +               &
                                      cpatch%dmax_leaf_psi(donc) * dnplant
         cpatch%dmin_leaf_psi(recc) = cpatch%dmin_leaf_psi(recc) * rnplant +               &
                                      cpatch%dmin_leaf_psi(donc) * dnplant
         cpatch%dmax_wood_psi(recc) = cpatch%dmax_wood_psi(recc) * rnplant +               &
                                      cpatch%dmax_wood_psi(donc) * dnplant
         cpatch%dmin_wood_psi(recc) = cpatch%dmin_wood_psi(recc) * rnplant +               &
                                      cpatch%dmin_wood_psi(donc) * dnplant
         !---------------------------------------------------------------------------------!


         !----- Recalculate psi from water_int. -------------------------------------------!
         call tw2psi(cpatch%fmean_leaf_water_int(recc),cpatch%fmean_wood_water_int(recc)   &
                    ,cpatch%is_small(recc),cpatch%bleaf(recc),cpatch%bsapwooda(recc)       &
                    ,cpatch%bsapwoodb(recc),cpatch%bdeada(recc),cpatch%bdeadb(recc)        &
                    ,cpatch%broot(recc),cpatch%dbh(recc),cpatch%pft(recc)                  &
                    ,cpatch%fmean_leaf_psi(recc),cpatch%fmean_wood_psi(recc) )
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
            call uextcm2tl( cpatch%fmean_leaf_energy   (recc)                              &
                          , cpatch%fmean_leaf_water    (recc)                              &
                          + cpatch%fmean_leaf_water_im2(recc)                              &
                          , cpatch%fmean_leaf_hcap     (recc)                              &
                          , cpatch%fmean_leaf_temp     (recc)                              &
                          , cpatch%fmean_leaf_fliq     (recc) )
            !------------------------------------------------------------------------------!


            !----- Scale vapour pressure deficit using LAI. -------------------------------!
            cpatch%fmean_leaf_vpdef   (recc) = cpatch%fmean_leaf_vpdef      (recc) * rlai  &
                                             + cpatch%fmean_leaf_vpdef      (donc) * dlai
            !------------------------------------------------------------------------------!
         else
            !----- None of the cohorts has leaf biomass use nplant to scale them. ---------!
            cpatch%fmean_leaf_temp (recc) = cpatch%fmean_leaf_temp (recc) * rnplant        &
                                          + cpatch%fmean_leaf_temp (donc) * dnplant
            cpatch%fmean_leaf_vpdef(recc) = cpatch%fmean_leaf_vpdef(recc) * rnplant        &
                                          + cpatch%fmean_leaf_vpdef(donc) * dnplant
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Always make liquid fraction consistent with temperature.                 !
            !------------------------------------------------------------------------------!
            if (cpatch%fmean_leaf_temp(recc) == t3ple) then
               cpatch%fmean_leaf_fliq(recc) = 0.5
            elseif (cpatch%fmean_leaf_temp(recc) > t3ple) then
               cpatch%fmean_leaf_fliq(recc) = 1.0
            else
               cpatch%fmean_leaf_fliq(recc) = 0.0
            end if
            !------------------------------------------------------------------------------!
         end if
         !------ Wood. --------------------------------------------------------------------!
         if ( cpatch%fmean_wood_hcap     (recc) > 0. .or.                                  &
              cpatch%fmean_wood_water_im2(recc) > 0.      ) then
            !----- Update temperature using the standard thermodynamics. ------------------!
            call uextcm2tl( cpatch%fmean_wood_energy   (recc)                              &
                          , cpatch%fmean_wood_water    (recc)                              &
                          + cpatch%fmean_wood_water_im2(recc)                              &
                          , cpatch%fmean_wood_hcap     (recc)                              &
                          , cpatch%fmean_wood_temp     (recc)                              &
                          , cpatch%fmean_wood_fliq     (recc) )
            !------------------------------------------------------------------------------!
         else                                                                              
            !----- Wood temperature can't be found using uextcm2tl (singularity). ---------!
            cpatch%fmean_wood_temp(recc) = cpatch%fmean_wood_temp(recc) * rnplant          &
                                         + cpatch%fmean_wood_temp(donc) * dnplant
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Always make liquid fraction consistent with temperature.                 !
            !------------------------------------------------------------------------------!
            if (cpatch%fmean_wood_temp(recc) == t3ple) then
               cpatch%fmean_wood_fliq(recc) = 0.5
            elseif (cpatch%fmean_wood_temp(recc) > t3ple) then
               cpatch%fmean_wood_fliq(recc) = 1.0
            else
               cpatch%fmean_wood_fliq(recc) = 0.0
            end if
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
         cpatch%dmean_nppleaf           (recc) = cpatch%dmean_nppleaf           (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_nppleaf           (donc)     &
                                               * dnplant
         cpatch%dmean_nppfroot          (recc) = cpatch%dmean_nppfroot          (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_nppfroot          (donc)     &
                                               * dnplant
         cpatch%dmean_nppsapwood        (recc) = cpatch%dmean_nppsapwood        (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_nppsapwood        (donc)     &
                                               * dnplant
         cpatch%dmean_nppbark           (recc) = cpatch%dmean_nppbark           (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_nppbark           (donc)     &
                                               * dnplant
         cpatch%dmean_nppcroot          (recc) = cpatch%dmean_nppcroot          (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_nppcroot          (donc)     &
                                               * dnplant
         cpatch%dmean_nppseeds          (recc) = cpatch%dmean_nppseeds          (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_nppseeds          (donc)     &
                                               * dnplant
         cpatch%dmean_nppwood           (recc) = cpatch%dmean_nppwood           (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_nppwood           (donc)     &
                                               * dnplant
         cpatch%dmean_nppdaily          (recc) = cpatch%dmean_nppdaily          (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_nppdaily          (donc)     &
                                               * dnplant
         cpatch%dmean_gpp               (recc) = cpatch%dmean_gpp               (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_gpp               (donc)     &
                                               * dnplant
         cpatch%dmean_npp               (recc) = cpatch%dmean_npp               (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_npp               (donc)     &
                                               * dnplant
         cpatch%dmean_leaf_resp         (recc) = cpatch%dmean_leaf_resp         (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_leaf_resp         (donc)     &
                                               * dnplant
         cpatch%dmean_root_resp         (recc) = cpatch%dmean_root_resp         (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_root_resp         (donc)     &
                                               * dnplant
         cpatch%dmean_stem_resp         (recc) = cpatch%dmean_stem_resp         (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_stem_resp         (donc)     &
                                               * dnplant
         cpatch%dmean_leaf_growth_resp  (recc) = cpatch%dmean_leaf_growth_resp  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_leaf_growth_resp  (donc)     &
                                               * dnplant
         cpatch%dmean_root_growth_resp  (recc) = cpatch%dmean_root_growth_resp  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_root_growth_resp  (donc)     &
                                               * dnplant
         cpatch%dmean_sapa_growth_resp  (recc) = cpatch%dmean_sapa_growth_resp  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_sapa_growth_resp  (donc)     &
                                               * dnplant
         cpatch%dmean_sapb_growth_resp  (recc) = cpatch%dmean_sapb_growth_resp  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_sapb_growth_resp  (donc)     &
                                               * dnplant
         cpatch%dmean_barka_growth_resp (recc) = cpatch%dmean_barka_growth_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_barka_growth_resp (donc)     &
                                               * dnplant
         cpatch%dmean_barkb_growth_resp (recc) = cpatch%dmean_barkb_growth_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_barkb_growth_resp (donc)     &
                                               * dnplant
         cpatch%dmean_leaf_storage_resp (recc) = cpatch%dmean_leaf_storage_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_leaf_storage_resp (donc)     &
                                               * dnplant
         cpatch%dmean_root_storage_resp (recc) = cpatch%dmean_root_storage_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_root_storage_resp (donc)     &
                                               * dnplant
         cpatch%dmean_sapa_storage_resp (recc) = cpatch%dmean_sapa_storage_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_sapa_storage_resp (donc)     &
                                               * dnplant
         cpatch%dmean_sapb_storage_resp (recc) = cpatch%dmean_sapb_storage_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_sapb_storage_resp (donc)     &
                                               * dnplant
         cpatch%dmean_barka_storage_resp(recc) = cpatch%dmean_barka_storage_resp(recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_barka_storage_resp(donc)     &
                                               * dnplant
         cpatch%dmean_barkb_storage_resp(recc) = cpatch%dmean_barkb_storage_resp(recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_barkb_storage_resp(donc)     &
                                               * dnplant
         cpatch%dmean_plresp            (recc) = cpatch%dmean_plresp            (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_plresp            (donc)     &
                                               * dnplant
         cpatch%dmean_light_level       (recc) = cpatch%dmean_light_level       (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_light_level       (donc)     &
                                               * dnplant
         cpatch%dmean_light_level_beam  (recc) = cpatch%dmean_light_level_beam  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_light_level_beam  (donc)     &
                                               * dnplant
         cpatch%dmean_light_level_diff  (recc) = cpatch%dmean_light_level_diff  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_light_level_diff  (donc)     &
                                               * dnplant
         cpatch%dmean_par_level_beam    (recc) = cpatch%dmean_par_level_beam    (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_par_level_beam    (donc)     &
                                               * dnplant
         cpatch%dmean_par_level_diffd   (recc) = cpatch%dmean_par_level_diffd   (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_par_level_diffd   (donc)     &
                                               * dnplant
         cpatch%dmean_par_level_diffu   (recc) = cpatch%dmean_par_level_diffu   (recc)     &
                                               * rnplant                                   &
                                               + cpatch%dmean_par_level_diffu   (donc)     &
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
         !    Plant hydrodynamics characteristics (XXT).  Internal water content and water !
         ! fluxes are weighted by nplant.  Area-based internal water is added.             !
         !---------------------------------------------------------------------------------!
         cpatch%dmean_leaf_water_int(recc) = cpatch%dmean_leaf_water_int(recc) * rnplant   &
                                           + cpatch%dmean_leaf_water_int(donc) * dnplant
         cpatch%dmean_wood_water_int(recc) = cpatch%dmean_wood_water_int(recc) * rnplant   &
                                           + cpatch%dmean_wood_water_int(donc) * dnplant
         cpatch%dmean_leaf_water_im2(recc) = cpatch%dmean_leaf_water_im2(recc)             &
                                           + cpatch%dmean_leaf_water_im2(donc)
         cpatch%dmean_wood_water_im2(recc) = cpatch%dmean_wood_water_im2(recc)             &
                                           + cpatch%dmean_wood_water_im2(donc)
         cpatch%dmean_wflux_wl      (recc) = cpatch%dmean_wflux_wl      (recc) * rnplant   &
                                           + cpatch%dmean_wflux_wl      (donc) * dnplant
         cpatch%dmean_wflux_gw      (recc) = cpatch%dmean_wflux_gw      (recc) * rnplant   &
                                           + cpatch%dmean_wflux_gw      (donc) * dnplant
         do isl = 1,nzg
            cpatch%dmean_wflux_gw_layer(isl,recc) =                                        &
                cpatch%dmean_wflux_gw_layer(isl,recc) * rnplant                            &
              + cpatch%dmean_wflux_gw_layer(isl,donc) * dnplant
         end do
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
            call uextcm2tl( cpatch%dmean_leaf_energy   (recc)                              &
                          , cpatch%dmean_leaf_water    (recc)                              &
                          + cpatch%dmean_leaf_water_im2(recc)                              &
                          , cpatch%dmean_leaf_hcap     (recc)                              &
                          , cpatch%dmean_leaf_temp     (recc)                              &
                          , cpatch%dmean_leaf_fliq     (recc) )
            !------------------------------------------------------------------------------!


            !----- Scale vapour pressure deficit using LAI. -------------------------------!
            cpatch%dmean_leaf_vpdef   (recc) = cpatch%dmean_leaf_vpdef      (recc) * rlai  &
                                             + cpatch%dmean_leaf_vpdef      (donc) * dlai
            !------------------------------------------------------------------------------!
         else
            !----- None of the cohorts has leaf biomass use nplant to scale them. ---------!
            cpatch%dmean_leaf_temp (recc) = cpatch%dmean_leaf_temp (recc) * rnplant        &
                                          + cpatch%dmean_leaf_temp (donc) * dnplant
            cpatch%dmean_leaf_vpdef(recc) = cpatch%dmean_leaf_vpdef(recc) * rnplant        &
                                          + cpatch%dmean_leaf_vpdef(donc) * dnplant
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Always make liquid fraction consistent with temperature.                 !
            !------------------------------------------------------------------------------!
            if (cpatch%dmean_leaf_temp(recc) == t3ple) then
               cpatch%dmean_leaf_fliq(recc) = 0.5
            elseif (cpatch%dmean_leaf_temp(recc) > t3ple) then
               cpatch%dmean_leaf_fliq(recc) = 1.0
            else
               cpatch%dmean_leaf_fliq(recc) = 0.0
            end if
            !------------------------------------------------------------------------------!
         end if
         !------ Wood. --------------------------------------------------------------------!
         if ( cpatch%dmean_wood_hcap     (recc) > 0. .or.                                  &
              cpatch%dmean_wood_water_im2(recc) > 0.      ) then
            !----- Update temperature using the standard thermodynamics. ------------------!
            call uextcm2tl( cpatch%dmean_wood_energy   (recc)                              &
                          , cpatch%dmean_wood_water    (recc)                              &
                          + cpatch%dmean_wood_water_im2(recc)                              &
                          , cpatch%dmean_wood_hcap     (recc)                              &
                          , cpatch%dmean_wood_temp     (recc)                              &
                          , cpatch%dmean_wood_fliq     (recc) )
            !------------------------------------------------------------------------------!
         else                                                                              
            !----- Wood temperature can't be found using uextcm2tl (singularity). ---------!
            cpatch%dmean_wood_temp(recc) = cpatch%dmean_wood_temp(recc) * rnplant          &
                                         + cpatch%dmean_wood_temp(donc) * dnplant
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Always make liquid fraction consistent with temperature.                 !
            !------------------------------------------------------------------------------!
            if (cpatch%dmean_wood_temp(recc) == t3ple) then
               cpatch%dmean_wood_fliq(recc) = 0.5
            elseif (cpatch%dmean_wood_temp(recc) > t3ple) then
               cpatch%dmean_wood_fliq(recc) = 1.0
            else
               cpatch%dmean_wood_fliq(recc) = 0.0
            end if
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
                                                   , dnplant                               &
                                                   , corr_cohort, .false.)
         cpatch%mmsqu_npp        (recc) = fuse_msqu( cpatch%mmean_npp        (recc)        &
                                                   , cpatch%mmsqu_npp        (recc)        &
                                                   , rnplant                               &
                                                   , cpatch%mmean_npp        (donc)        &
                                                   , cpatch%mmsqu_npp        (donc)        &
                                                   , dnplant                               &
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
         cpatch%mmsqu_wflux_wl   (recc) = fuse_msqu( cpatch%mmean_wflux_wl   (recc)        &
                                                   , cpatch%mmsqu_wflux_wl   (recc)        &
                                                   , rnplant                               &
                                                   , cpatch%mmean_wflux_wl   (donc)        &
                                                   , cpatch%mmsqu_wflux_wl   (donc)        &
                                                   , dnplant                               &
                                                   , corr_cohort, .false.) !kg/pl/s
         cpatch%mmsqu_wflux_gw   (recc) = fuse_msqu( cpatch%mmean_wflux_gw   (recc)        &
                                                   , cpatch%mmsqu_wflux_gw   (recc)        &
                                                   , rnplant                               &
                                                   , cpatch%mmean_wflux_gw   (donc)        &
                                                   , cpatch%mmsqu_wflux_gw   (donc)        &
                                                   , dnplant                               &
                                                   , corr_cohort, .false.) !kg/pl/s
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
         cpatch%mmean_nppleaf           (recc) = cpatch%mmean_nppleaf           (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_nppleaf           (donc)     &
                                               * dnplant
         cpatch%mmean_nppfroot          (recc) = cpatch%mmean_nppfroot          (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_nppfroot          (donc)     &
                                               * dnplant
         cpatch%mmean_nppsapwood        (recc) = cpatch%mmean_nppsapwood        (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_nppsapwood        (donc)     &
                                               * dnplant
         cpatch%mmean_nppbark           (recc) = cpatch%mmean_nppbark           (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_nppbark           (donc)     &
                                               * dnplant
         cpatch%mmean_nppcroot          (recc) = cpatch%mmean_nppcroot          (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_nppcroot          (donc)     &
                                               * dnplant
         cpatch%mmean_nppseeds          (recc) = cpatch%mmean_nppseeds          (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_nppseeds          (donc)     &
                                               * dnplant
         cpatch%mmean_nppwood           (recc) = cpatch%mmean_nppwood           (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_nppwood           (donc)     &
                                               * dnplant
         cpatch%mmean_nppdaily          (recc) = cpatch%mmean_nppdaily          (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_nppdaily          (donc)     &
                                               * dnplant
         cpatch%mmean_gpp               (recc) = cpatch%mmean_gpp               (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_gpp               (donc)     &
                                               * dnplant
         cpatch%mmean_npp               (recc) = cpatch%mmean_npp               (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_npp               (donc)     &
                                               * dnplant
         cpatch%mmean_leaf_resp         (recc) = cpatch%mmean_leaf_resp         (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_leaf_resp         (donc)     &
                                               * dnplant
         cpatch%mmean_root_resp         (recc) = cpatch%mmean_root_resp         (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_root_resp         (donc)     &
                                               * dnplant
         cpatch%mmean_stem_resp         (recc) = cpatch%mmean_stem_resp         (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_stem_resp         (donc)     &
                                               * dnplant
         cpatch%mmean_leaf_growth_resp  (recc) = cpatch%mmean_leaf_growth_resp  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_leaf_growth_resp  (donc)     &
                                               * dnplant
         cpatch%mmean_root_growth_resp  (recc) = cpatch%mmean_root_growth_resp  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_root_growth_resp  (donc)     &
                                               * dnplant
         cpatch%mmean_sapa_growth_resp  (recc) = cpatch%mmean_sapa_growth_resp  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_sapa_growth_resp  (donc)     &
                                               * dnplant
         cpatch%mmean_sapb_growth_resp  (recc) = cpatch%mmean_sapb_growth_resp  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_sapb_growth_resp  (donc)     &
                                               * dnplant
         cpatch%mmean_barka_growth_resp (recc) = cpatch%mmean_barka_growth_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_barka_growth_resp (donc)     &
                                               * dnplant
         cpatch%mmean_barkb_growth_resp (recc) = cpatch%mmean_barkb_growth_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_barkb_growth_resp (donc)     &
                                               * dnplant
         cpatch%mmean_leaf_storage_resp (recc) = cpatch%mmean_leaf_storage_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_leaf_storage_resp (donc)     &
                                               * dnplant
         cpatch%mmean_root_storage_resp (recc) = cpatch%mmean_root_storage_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_root_storage_resp (donc)     &
                                               * dnplant
         cpatch%mmean_sapa_storage_resp (recc) = cpatch%mmean_sapa_storage_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_sapa_storage_resp (donc)     &
                                               * dnplant
         cpatch%mmean_sapb_storage_resp (recc) = cpatch%mmean_sapb_storage_resp (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_sapb_storage_resp (donc)     &
                                               * dnplant
         cpatch%mmean_barka_storage_resp(recc) = cpatch%mmean_barka_storage_resp(recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_barka_storage_resp(donc)     &
                                               * dnplant
         cpatch%mmean_barkb_storage_resp(recc) = cpatch%mmean_barkb_storage_resp(recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_barkb_storage_resp(donc)     &
                                               * dnplant
         cpatch%mmean_plresp            (recc) = cpatch%mmean_plresp            (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_plresp            (donc)     &
                                               * dnplant
         cpatch%mmean_light_level       (recc) = cpatch%mmean_light_level       (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_light_level       (donc)     &
                                               * dnplant
         cpatch%mmean_light_level_beam  (recc) = cpatch%mmean_light_level_beam  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_light_level_beam  (donc)     &
                                               * dnplant
         cpatch%mmean_light_level_diff  (recc) = cpatch%mmean_light_level_diff  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_light_level_diff  (donc)     &
                                               * dnplant
         cpatch%mmean_bleaf             (recc) = cpatch%mmean_bleaf             (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_bleaf             (donc)     &
                                               * dnplant
         cpatch%mmean_broot             (recc) = cpatch%mmean_broot             (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_broot             (donc)     &
                                               * dnplant
         cpatch%mmean_bbarka            (recc) = cpatch%mmean_bbarka            (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_bbarka            (donc)     &
                                               * dnplant
         cpatch%mmean_bbarkb            (recc) = cpatch%mmean_bbarkb            (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_bbarkb            (donc)     &
                                               * dnplant
         cpatch%mmean_balive            (recc) = cpatch%mmean_balive            (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_balive            (donc)     &
                                               * dnplant
         cpatch%mmean_bstorage          (recc) = cpatch%mmean_bstorage          (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_bstorage          (donc)     &
                                               * dnplant
         cpatch%mmean_leaf_maintenance  (recc) = cpatch%mmean_leaf_maintenance  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_leaf_maintenance  (donc)     &
                                               * dnplant
         cpatch%mmean_root_maintenance  (recc) = cpatch%mmean_root_maintenance  (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_root_maintenance  (donc)     &
                                               * dnplant
         cpatch%mmean_barka_maintenance (recc) = cpatch%mmean_barka_maintenance (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_barka_maintenance (donc)     &
                                               * dnplant
         cpatch%mmean_barkb_maintenance (recc) = cpatch%mmean_barkb_maintenance (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_barkb_maintenance (donc)     &
                                               * dnplant
         cpatch%mmean_leaf_drop         (recc) = cpatch%mmean_leaf_drop         (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_leaf_drop         (donc)     &
                                               * dnplant
         cpatch%mmean_root_drop         (recc) = cpatch%mmean_root_drop         (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_root_drop         (donc)     &
                                               * dnplant
         cpatch%mmean_cb                (recc) = cpatch%mmean_cb                (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_cb                (donc)     &
                                               * dnplant
         cpatch%mmean_par_level_beam    (recc) = cpatch%mmean_par_level_beam    (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_par_level_beam    (donc)     &
                                               * dnplant
         cpatch%mmean_par_level_diffd   (recc) = cpatch%mmean_par_level_diffd   (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_par_level_diffd   (donc)     &
                                               * dnplant
         cpatch%mmean_par_level_diffu   (recc) = cpatch%mmean_par_level_diffu   (recc)     &
                                               * rnplant                                   &
                                               + cpatch%mmean_par_level_diffu   (donc)     &
                                               * dnplant
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Bark thickness is weighted by basal area.                                    !
         !---------------------------------------------------------------------------------!
         cpatch%mmean_thbark          (recc) = cpatch%mmean_thbark          (recc)         &
                                             * rba                                         &
                                             + cpatch%mmean_thbark          (donc)         &
                                             * dba
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
         !     Light-phenology characteristics.  To conserve maintenance costs, leaf       !
         !  longevity must be scaled with leaf biomass (but note that the scaling must be  !
         !  applied to turnover instead of leaf lifespan).  Likewise, SLA is in            !
         !  m2_leaf/kgC, so it must be scaled with leaf biomass.  Carboxylation rate       !
         !  (vm_bar) and dark respiration rate (rd_bar) are in umol/m2_leaf/s, so          !
         ! They must be scaled with LAI.                                                        !
         !---------------------------------------------------------------------------------!
         cpatch%mmean_sla             (recc) = cpatch%mmean_sla         (recc) * rbleaf    &
                                             + cpatch%mmean_sla         (donc) * dbleaf
         cpatch%mmean_vm_bar          (recc) = cpatch%mmean_vm_bar      (recc) * rlai      &
                                             + cpatch%mmean_vm_bar      (donc) * dlai
         cpatch%mmean_rd_bar          (recc) = cpatch%mmean_rd_bar      (recc) * rlai      &
                                             + cpatch%mmean_rd_bar      (donc) * dlai
         !------ For Life span, we must check whether they are non-zero. ------------------!
         if ( abs(cpatch%mmean_llspan(recc)*cpatch%mmean_llspan(donc)) > tiny_num ) then
            !------------------------------------------------------------------------------!
            !     The denominator weights are not inadvertenly swapped.  This happens      !
            ! because we are linearly scaling turnover, not leaf life span, to ensure that !
            ! maintenance costs are consistent.                                            !
            !------------------------------------------------------------------------------!
            cpatch%mmean_llspan(recc) =   cpatch%mmean_llspan(recc)                        &
                                      *   cpatch%mmean_llspan(donc)                        &
                                      / ( cpatch%mmean_llspan(recc) * dbleaf               &
                                        + cpatch%mmean_llspan(donc) * rbleaf )
            !------------------------------------------------------------------------------!
         else
            !------ This only happens when both are zero, so it really doesn't matter. ----!
            cpatch%mmean_llspan(recc) = cpatch%mmean_llspan(recc) * rbleaf                 &
                                      + cpatch%mmean_llspan(donc) * dbleaf
            !------------------------------------------------------------------------------!
         end if
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
         !    Plant hydrodynamics characteristics (XXT).  Internal water content and water !
         ! fluxes are weighted by nplant.  Area-based internal water is added.             !
         !---------------------------------------------------------------------------------!
         cpatch%mmean_leaf_water_int(recc) = cpatch%mmean_leaf_water_int(recc) * rnplant   &
                                           + cpatch%mmean_leaf_water_int(donc) * dnplant
         cpatch%mmean_wood_water_int(recc) = cpatch%mmean_wood_water_int(recc) * rnplant   &
                                           + cpatch%mmean_wood_water_int(donc) * dnplant
         cpatch%mmean_leaf_water_im2(recc) = cpatch%mmean_leaf_water_im2(recc)             &
                                           + cpatch%mmean_leaf_water_im2(donc)
         cpatch%mmean_wood_water_im2(recc) = cpatch%mmean_wood_water_im2(recc)             &
                                           + cpatch%mmean_wood_water_im2(donc)
         cpatch%mmean_wflux_gw      (recc) = cpatch%mmean_wflux_gw      (recc) * rnplant   &
                                           + cpatch%mmean_wflux_gw      (donc) * dnplant
         cpatch%mmean_wflux_wl      (recc) = cpatch%mmean_wflux_wl      (recc) * rnplant   &
                                           + cpatch%mmean_wflux_wl      (donc) * dnplant
         do isl = 1,nzg
            cpatch%mmean_wflux_gw_layer(isl,recc) =                                        &
                cpatch%mmean_wflux_gw_layer(isl,recc) * rnplant                            &
              + cpatch%mmean_wflux_gw_layer(isl,donc) * dnplant
         end do
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
            call uextcm2tl( cpatch%mmean_leaf_energy   (recc)                              &
                          , cpatch%mmean_leaf_water    (recc)                              &
                          + cpatch%mmean_leaf_water_im2(recc)                              &
                          , cpatch%mmean_leaf_hcap     (recc)                              &
                          , cpatch%mmean_leaf_temp     (recc)                              &
                          , cpatch%mmean_leaf_fliq     (recc) )
            !------------------------------------------------------------------------------!


            !----- Scale vapour pressure deficit using LAI. -------------------------------!
            cpatch%mmean_leaf_vpdef   (recc) = cpatch%mmean_leaf_vpdef      (recc) * rlai  &
                                             + cpatch%mmean_leaf_vpdef      (donc) * dlai
            !------------------------------------------------------------------------------!
         else
            !----- None of the cohorts has leaf biomass use nplant to scale them. ---------!
            cpatch%mmean_leaf_temp (recc) = cpatch%mmean_leaf_temp (recc) * rnplant        &
                                          + cpatch%mmean_leaf_temp (donc) * dnplant
            cpatch%mmean_leaf_vpdef(recc) = cpatch%mmean_leaf_vpdef(recc) * rnplant        &
                                          + cpatch%mmean_leaf_vpdef(donc) * dnplant
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Always make liquid fraction consistent with temperature.                 !
            !------------------------------------------------------------------------------!
            if (cpatch%mmean_leaf_temp(recc) == t3ple) then
               cpatch%mmean_leaf_fliq(recc) = 0.5
            elseif (cpatch%mmean_leaf_temp(recc) > t3ple) then
               cpatch%mmean_leaf_fliq(recc) = 1.0
            else
               cpatch%mmean_leaf_fliq(recc) = 0.0
            end if
            !------------------------------------------------------------------------------!
         end if
         !------ Wood. --------------------------------------------------------------------!
         if ( cpatch%mmean_wood_hcap     (recc) > 0. .or.                                  &
              cpatch%mmean_wood_water_im2(recc) > 0.      ) then
            !----- Update temperature using the standard thermodynamics. ------------------!
            call uextcm2tl( cpatch%mmean_wood_energy   (recc)                              &
                          , cpatch%mmean_wood_water    (recc)                              &
                          + cpatch%mmean_wood_water_im2(recc)                              &
                          , cpatch%mmean_wood_hcap     (recc)                              &
                          , cpatch%mmean_wood_temp     (recc)                              &
                          , cpatch%mmean_wood_fliq     (recc) )
            !------------------------------------------------------------------------------!
         else                                                                              
            !----- Wood temperature can't be found using uextcm2tl (singularity). ---------!
            cpatch%mmean_wood_temp(recc) = cpatch%mmean_wood_temp(recc) * rnplant          &
                                         + cpatch%mmean_wood_temp(donc) * dnplant
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Always make liquid fraction consistent with temperature.                 !
            !------------------------------------------------------------------------------!
            if (cpatch%mmean_wood_temp(recc) == t3ple) then
               cpatch%mmean_wood_fliq(recc) = 0.5
            elseif (cpatch%mmean_wood_temp(recc) > t3ple) then
               cpatch%mmean_wood_fliq(recc) = 1.0
            else
               cpatch%mmean_wood_fliq(recc) = 0.0
            end if
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
            cpatch%qmsqu_wflux_wl   (t,recc) = fuse_msqu(cpatch%qmean_wflux_wl   (t,recc)  &
                                                        ,cpatch%qmsqu_wflux_wl   (t,recc)  &
                                                        ,rnplant                           &
                                                        ,cpatch%qmean_wflux_wl   (t,donc)  &
                                                        ,cpatch%qmsqu_wflux_wl   (t,donc)  &
                                                        ,dnplant                           &
                                                        ,corr_cohort, .false.)
            cpatch%qmsqu_wflux_gw   (t,recc) = fuse_msqu(cpatch%qmean_wflux_gw   (t,recc)  &
                                                        ,cpatch%qmsqu_wflux_gw   (t,recc)  &
                                                        ,rnplant                           &
                                                        ,cpatch%qmean_wflux_gw   (t,donc)  &
                                                        ,cpatch%qmsqu_wflux_gw   (t,donc)  &
                                                        ,dnplant                           &
                                                        ,corr_cohort, .false.)
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
         cpatch%qmean_gpp               (:,recc) = cpatch%qmean_gpp               (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_gpp               (:,donc) &
                                                 * dnplant
         cpatch%qmean_npp               (:,recc) = cpatch%qmean_npp               (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_npp               (:,donc) &
                                                 * dnplant
         cpatch%qmean_leaf_resp         (:,recc) = cpatch%qmean_leaf_resp         (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_leaf_resp         (:,donc) &
                                                 * dnplant
         cpatch%qmean_root_resp         (:,recc) = cpatch%qmean_root_resp         (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_root_resp         (:,donc) &
                                                 * dnplant
         cpatch%qmean_stem_resp         (:,recc) = cpatch%qmean_stem_resp         (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_stem_resp         (:,donc) &
                                                 * dnplant
         cpatch%qmean_leaf_growth_resp  (:,recc) = cpatch%qmean_leaf_growth_resp  (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_leaf_growth_resp  (:,donc) &
                                                 * dnplant
         cpatch%qmean_root_growth_resp  (:,recc) = cpatch%qmean_root_growth_resp  (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_root_growth_resp  (:,donc) &
                                                 * dnplant
         cpatch%qmean_sapa_growth_resp  (:,recc) = cpatch%qmean_sapa_growth_resp  (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_sapa_growth_resp  (:,donc) &
                                                 * dnplant
         cpatch%qmean_sapb_growth_resp  (:,recc) = cpatch%qmean_sapb_growth_resp  (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_sapb_growth_resp  (:,donc) &
                                                 * dnplant
         cpatch%qmean_barka_growth_resp (:,recc) = cpatch%qmean_barka_growth_resp (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_barka_growth_resp (:,donc) &
                                                 * dnplant
         cpatch%qmean_barkb_growth_resp (:,recc) = cpatch%qmean_barkb_growth_resp (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_barkb_growth_resp (:,donc) &
                                                 * dnplant
         cpatch%qmean_leaf_storage_resp (:,recc) = cpatch%qmean_leaf_storage_resp (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_leaf_storage_resp (:,donc) &
                                                 * dnplant
         cpatch%qmean_root_storage_resp (:,recc) = cpatch%qmean_root_storage_resp (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_root_storage_resp (:,donc) &
                                                 * dnplant
         cpatch%qmean_sapa_storage_resp (:,recc) = cpatch%qmean_sapa_storage_resp (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_sapa_storage_resp (:,donc) &
                                                 * dnplant
         cpatch%qmean_sapb_storage_resp (:,recc) = cpatch%qmean_sapb_storage_resp (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_sapb_storage_resp (:,donc) &
                                                 * dnplant
         cpatch%qmean_barka_storage_resp(:,recc) = cpatch%qmean_barka_storage_resp(:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_barka_storage_resp(:,donc) &
                                                 * dnplant
         cpatch%qmean_barkb_storage_resp(:,recc) = cpatch%qmean_barkb_storage_resp(:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_barkb_storage_resp(:,donc) &
                                                 * dnplant
         cpatch%qmean_plresp            (:,recc) = cpatch%qmean_plresp            (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_plresp            (:,donc) &
                                                 * dnplant
         cpatch%qmean_light_level       (:,recc) = cpatch%qmean_light_level       (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_light_level       (:,donc) &
                                                 * dnplant
         cpatch%qmean_light_level_beam  (:,recc) = cpatch%qmean_light_level_beam  (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_light_level_beam  (:,donc) &
                                                 * dnplant
         cpatch%qmean_light_level_diff  (:,recc) = cpatch%qmean_light_level_diff  (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_light_level_diff  (:,donc) &
                                                 * dnplant
         cpatch%qmean_par_level_beam    (:,recc) = cpatch%qmean_par_level_beam    (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_par_level_beam    (:,donc) &
                                                 * dnplant
         cpatch%qmean_par_level_diffd   (:,recc) = cpatch%qmean_par_level_diffd   (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_par_level_diffd   (:,donc) &
                                                 * dnplant
         cpatch%qmean_par_level_diffu   (:,recc) = cpatch%qmean_par_level_diffu   (:,recc) &
                                                 * rnplant                                 &
                                                 + cpatch%qmean_par_level_diffu   (:,donc) &
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
         !    Plant hydrodynamics characteristics (XXT).  Internal water content and water !
         ! fluxes are weighted by nplant.  Area-based internal water is added.             !
         !---------------------------------------------------------------------------------!
         cpatch%qmean_leaf_water_int(:,recc) = cpatch%qmean_leaf_water_int(:,recc)         &
                                             * rnplant                                     &
                                             + cpatch%qmean_leaf_water_int(:,donc)         &
                                             * dnplant
         cpatch%qmean_wood_water_int(:,recc) = cpatch%qmean_wood_water_int(:,recc)         &
                                             * rnplant                                     &
                                             + cpatch%qmean_wood_water_int(:,donc)         &
                                             * dnplant
         cpatch%qmean_leaf_water_im2(:,recc) = cpatch%qmean_leaf_water_im2(:,recc)         &
                                             + cpatch%qmean_leaf_water_im2(:,donc)
         cpatch%qmean_wood_water_im2(:,recc) = cpatch%qmean_wood_water_im2(:,recc)         &
                                             + cpatch%qmean_wood_water_im2(:,donc)
         !----- Water fluxes are also weighted by nplant since they are kg H2O/s/plant. ---!
         cpatch%qmean_wflux_gw      (:,recc) = cpatch%qmean_wflux_gw   (:,recc) * rnplant  &
                                             + cpatch%qmean_wflux_gw   (:,donc) * dnplant
         cpatch%qmean_wflux_wl      (:,recc) = cpatch%qmean_wflux_wl   (:,recc) * rnplant  &
                                             + cpatch%qmean_wflux_wl   (:,donc) * dnplant
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     For psi, we simply use lai-weighted average for leaves and bdead-weighted   !
         ! average for wood because at this time point we don't know qmean biomass values. !
         !                                                                                 !
         ! MLO -> XX wood_psi is currently scaled by LAI too, so not consistent with the   !
         !        comment above.  Should it be scaled by WAI instead?                      !
         !---------------------------------------------------------------------------------!
         cpatch%qmean_leaf_psi      (:,recc) = cpatch%qmean_leaf_psi   (:,recc) * rlai     &
                                             + cpatch%qmean_leaf_psi   (:,donc) * dlai
         cpatch%qmean_wood_psi      (:,recc) = cpatch%qmean_wood_psi   (:,recc) * rlai     &
                                             + cpatch%qmean_wood_psi   (:,donc) * dlai
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
               !---------------------------------------------------------------------------!
               !    Update temperature and liquid fraction using standard thermodynamics.  !
               !---------------------------------------------------------------------------!
               call uextcm2tl( cpatch%qmean_leaf_energy   (t,recc)                         &
                             , cpatch%qmean_leaf_water    (t,recc)                         &
                             + cpatch%qmean_leaf_water_im2(t,recc)                         &
                             , cpatch%qmean_leaf_hcap     (t,recc)                         &
                             , cpatch%qmean_leaf_temp     (t,recc)                         &
                             , cpatch%qmean_leaf_fliq     (t,recc) )
               !---------------------------------------------------------------------------!


               !----- Scale vapour pressure deficit using LAI. ----------------------------!
               cpatch%qmean_leaf_vpdef   (t,recc) = cpatch%qmean_leaf_vpdef(t,recc) * rlai &
                                                  + cpatch%qmean_leaf_vpdef(t,donc) * dlai
               !------------------------------------------------------------------------------!
            else
               !----- None of the cohorts has leaf biomass use nplant to scale them. ------!
               cpatch%qmean_leaf_temp (t,recc) = cpatch%qmean_leaf_temp (t,recc) * rnplant &
                                               + cpatch%qmean_leaf_temp (t,donc) * dnplant
               cpatch%qmean_leaf_vpdef(t,recc) = cpatch%qmean_leaf_vpdef(t,recc) * rnplant &
                                               + cpatch%qmean_leaf_vpdef(t,donc) * dnplant
               !------------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Always make liquid fraction consistent with temperature.              !
               !---------------------------------------------------------------------------!
               if (cpatch%qmean_leaf_temp(t,recc) == t3ple) then
                  cpatch%qmean_leaf_fliq(t,recc) = 0.5
               elseif (cpatch%qmean_leaf_temp(t,recc) > t3ple) then
                  cpatch%qmean_leaf_fliq(t,recc) = 1.0
               else
                  cpatch%qmean_leaf_fliq(t,recc) = 0.0
               end if
               !---------------------------------------------------------------------------!
            end if
            !------ Wood. --------------------------------------------------------------------!
            if ( cpatch%qmean_wood_hcap     (t,recc) > 0. .or.                                  &
                 cpatch%qmean_wood_water_im2(t,recc) > 0.      ) then
               !----- Update temperature using the standard thermodynamics. ---------------!
               call uextcm2tl( cpatch%qmean_wood_energy   (t,recc)                         &
                             , cpatch%qmean_wood_water    (t,recc)                         &
                             + cpatch%qmean_wood_water_im2(t,recc)                         &
                             , cpatch%qmean_wood_hcap     (t,recc)                         &
                             , cpatch%qmean_wood_temp     (t,recc)                         &
                             , cpatch%qmean_wood_fliq     (t,recc) )
               !---------------------------------------------------------------------------!
            else
               !----- Wood temperature can't be found using uextcm2tl (singularity). ------!
               cpatch%qmean_wood_temp(t,recc) = cpatch%qmean_wood_temp(t,recc) * rnplant   &
                                              + cpatch%qmean_wood_temp(t,donc) * dnplant
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Always make liquid fraction consistent with temperature.              !
               !---------------------------------------------------------------------------!
               if (cpatch%qmean_wood_temp(t,recc) == t3ple) then
                  cpatch%qmean_wood_fliq(t,recc) = 0.5
               elseif (cpatch%qmean_wood_temp(t,recc) > t3ple) then
                  cpatch%qmean_wood_fliq(t,recc) = 1.0
               else
                  cpatch%qmean_wood_fliq(t,recc) = 0.0
               end if
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end if
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
   subroutine new_fuse_patches(cgrid,ifm,fuse_initial)
      use update_derived_utils, only : patch_pft_size_profile  ! ! sub-routine
      use ed_state_vars       , only : edtype                  & ! structure
                                     , polygontype             & ! structure
                                     , sitetype                & ! structure
                                     , patchtype               ! ! structure
      use fusion_fission_coms , only : ff_nhgt                 & ! intent(in)
                                     , niter_patfus            & ! intent(in)
                                     , print_fuse_details      & ! intent(in)
                                     , hgt_class               & ! intent(in)
                                     , pat_light_ext           & ! intent(in)
                                     , pat_light_tol_min       & ! intent(in)
                                     , pat_light_tol_max       & ! intent(in)
                                     , pat_light_tol_mult      & ! intent(in)
                                     , pat_light_mxd_fac       & ! intent(in)
                                     , pat_diff_age_tol        & ! intent(in)
                                     , pat_min_area_remain     & ! intent(in)
                                     , fuse_prefix             ! ! intent(in)
      use disturb_coms        , only : min_patch_area          & ! intent(in)
                                     , min_oldgrowth           ! ! intent(in)
      use ed_max_dims         , only : n_pft                   & ! intent(in)
                                     , str_len                 ! ! intent(in)
      use mem_polygons        , only : maxpatch                ! ! intent(in)
      use ed_node_coms        , only : mynum                   ! ! intent(in)
      use ed_misc_coms        , only : current_time            ! ! intent(in)
      use grid_coms           , only : nzg                     & ! intent(in)
                                     , nzs                     ! ! intent(in)
      use consts_coms         , only : lnexp_min               & ! intent(in)
                                     , lnexp_max               ! ! intent(in)
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
      logical, dimension(:) , allocatable :: fuse_table      ! Flag: still can be fused.
      logical, dimension(:) , allocatable :: remain          ! Flag: will remain
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
      integer                             :: npatches_remain ! # of patches that will remain
                                                             !   (i.e. excluding patches 
                                                             !   with small area that will
                                                             !   be terminated)
      integer                             :: npatches_new    ! New # of patches
      integer                             :: npatches_old    ! Old # of patches
      integer                             :: npatches_orig   ! Original # of patches
      integer                             :: ncohorts_remain ! # of cohorts that will remain
                                                             !   (i.e. excluding patches 
                                                             !   with small area that will
                                                             !   be terminated)
      integer                             :: npop            ! # of populated layers
      logical                             :: rec_old         ! Receptor patch is old
      logical                             :: don_old         ! Donor patch is old
      logical                             :: rec_pop         ! Receptor patch is not bare
      logical                             :: don_pop         ! Donor patch is not bare
      logical                             :: fuse_flag       ! Flag: fusion will happen
      logical                             :: recp_found      ! Found a receptor candidate
      logical                             :: same_age        ! Patches with same age
      logical                             :: same_lu         ! Patches with same age
      logical                             :: old_or_same_lu  ! Old patches or the same LU.
      real                                :: area_remain     ! Area of patches that will 
                                                             !    remain (i.e. excluding 
                                                             !    patches  with small area 
                                                             !    that will be terminated)
      real                                :: llevel_diff     ! Absolute difference in prof.
      real                                :: llevel_diff_max ! Maximum llevel_diff
      real                                :: llevel_diff_avg ! Maximum llevel_diff
      real                                :: hgt_diff_max    ! Height of maximum llevel_diff
      real                                :: lnexp_donp      ! LN Exp. of donor patch
      real                                :: lnexp_recp      ! LN Exp. of receptor  patch
      real                                :: llevel_donp     ! Light level of donor patch
      real                                :: llevel_recp     ! Light level of rec.  patch
      real                                :: dllev_diff_max  ! Donor LL @ max. diff.
      real                                :: rllev_diff_max  ! Receptor LL @ max. diff.
      real                                :: pat_light_mxd   ! Maximum deviation toler.
      real                                :: pat_light_tol   ! Light level Absolute toler.
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
      real                                :: dclai_diff_max  ! Donor CLAI @ max. diff.
      real                                :: rclai_diff_max  ! Receptor CLAI @ max. diff.
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
            allocate(remain(csite%npatches))

            !------------------------------------------------------------------------------!
            !     Allocate the fusion flag vector, and set all elements to .true., which   !
            ! means that every patch can be fused.  As soon as the patch is fused, we will !
            ! switch the flag to false.                                                    !
            !------------------------------------------------------------------------------!
            fuse_table(:) = .true.
            remain    (:) = fuse_table(:) .and. csite%area(:) >= min_patch_area
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
            same_age = .false.
            same_lu  = .false.
            donloop_check: do donp=csite%npatches,2,-1
               recloop_check: do recp=donp-1,1,-1
                  same_age = abs(csite%age(recp) - csite%age(donp)) <= pat_diff_age_tol
                  same_lu  = csite%dist_type(recp) == csite%dist_type(donp)
                  !----- At least two patches have the same age. --------------------------!
                  if (same_age .and. same_lu) exit donloop_check
                  !------------------------------------------------------------------------!
              end do recloop_check
            end do donloop_check
            !------------------------------------------------------------------------------!


            !----- Print same age status. -------------------------------------------------!
            if (print_fuse_details) then
                  open (unit=72,file=trim(fuse_fout),status='old',action='write'           &
                                        ,position='append')
                  write(unit=72,fmt='(a,1x,l1)') '     * same_age is ',same_age
                  write(unit=72,fmt='(a,1x,l1)') '     * same_lu is  ',same_lu
                  close(unit=72,status='keep')
            end if
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            donloope: do donp=csite%npatches,2,-1
               donpatch => csite%patch(donp)
               don_lu    = csite%dist_type(donp)
               don_old   = csite%age(donp) >= min_oldgrowth(don_lu)
               don_pop   = donpatch%ncohorts > 0
               
               !----- If patch is not empty, or has already been fused, move on. ----------!
               if ( (.not. fuse_table(donp)) .or. ( dont_force_fuse .and. don_pop) ) then
                  cycle donloope
               end if
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !     If we reach this point, it means that the donor patch is empty and    !
               ! hasn't been fused yet: look for an older empty patch and merge them.      !
               !---------------------------------------------------------------------------!
               recloope: do recp=donp-1,1,-1
                  recpatch => csite%patch(recp)
                  rec_lu  = csite%dist_type(recp)
                  rec_old = csite%age(recp) >= min_oldgrowth(rec_lu)
                  rec_pop = recpatch%ncohorts > 0

                  !------------------------------------------------------------------------!
                  !     Set this flag that checks whether the patches have the same        !
                  ! disturbance type or are too old so we don't need to distinguish them.  !
                  !------------------------------------------------------------------------!
                  old_or_same_lu = don_lu == rec_lu .or. ( don_old .and. rec_old)
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !     Skip the patch if it isn't empty, or it has already been fused, or !
                  ! if the donor and receptor have different disturbance types.            !
                  !------------------------------------------------------------------------!
                  if ( (.not. fuse_table(recp))                                       .or. &
                       (dont_force_fuse .and.  (rec_pop .or. (.not. old_or_same_lu))) ) then
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
                  call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)               &
                                     ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)               &
                                     ,cpoly%green_leaf_factor(:,isi),fuse_initial          &
                                     ,elim_nplant,elim_lai)
                  !------------------------------------------------------------------------!


                  !----- Record the fusion if requested by the user. ----------------------!
                  if (print_fuse_details) then
                     open (unit=72,file=trim(fuse_fout),status='old',action='write'        &
                                                        ,position='append')
                     write(unit=72,fmt='(2(a,i6),a)') '     * Patches ',donp,' and ',recp  &
                                                     ,' were fused (both were empty).'
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
            if (same_age .and. same_lu) then
               !----- Start with no multiplication factor for tolerance. ------------------!
               pat_light_tol  = pat_light_tol_min
               pat_light_mxd  = pat_light_tol * pat_light_mxd_fac
               !---------------------------------------------------------------------------!

               mainfuseloopa: do ifus=1,niter_patfus

                  npatches_old    = count(fuse_table)
                  npatches_new    = npatches_old

                  !------------------------------------------------------------------------!
                  !    Inform that the upcoming fusions are going to be with populated     !
                  ! patches, and record the tolerance used.                                !
                  !------------------------------------------------------------------------!
                  if (print_fuse_details) then
                     open (unit=72,file=trim(fuse_fout),status='old',action='write'        &
                                                        ,position='append')
                     write(unit=72,fmt='(a,1x,1(a,1x,es9.2,1x))')                          &
                         '   + Look for similar populated patches of same age and same LU' &
                        ,' - Rel. Tolerance  =',pat_light_tol
                     close(unit=72,status='keep')
                  end if
                  !------------------------------------------------------------------------!

                  donloopa: do donp=csite%npatches,2,-1
                     donpatch => csite%patch(donp)
                     don_lu = csite%dist_type(donp)
                     don_old = csite%age(donp) >= min_oldgrowth(don_lu)
                     don_pop = donpatch%ncohorts > 0
                     
                     !----- If patch is not empty, or has already been fused, move on. ----!
                     if ( (.not. fuse_table(donp))                    .or.                 &
                          ( dont_force_fuse .and. (.not. don_pop) ) ) then
                        cycle donloopa
                     end if
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     If we reach this point, it means that the donor patch is        !
                     ! populated and hasn't been fused yet: look for other patches with    !
                     ! the same age and merge them.                                        !
                     !---------------------------------------------------------------------!
                     recloopa: do recp=donp-1,1,-1
                        recpatch => csite%patch(recp)
                        rec_lu   = csite%dist_type(recp)
                        rec_old  = csite%age(recp) >= min_oldgrowth(rec_lu)
                        rec_pop  = recpatch%ncohorts > 0
                        same_age = abs(csite%age(recp)-csite%age(donp)) <= pat_diff_age_tol
                        same_lu  = csite%dist_type(recp) == csite%dist_type(donp)


                        !------------------------------------------------------------------!
                        !     Skip the patch if it isn't empty, or it has already been     !
                        ! fused, or if the donor and receptor have different disturbance   !
                        ! types.                                                           !
                        !------------------------------------------------------------------!
                        if ( (.not. fuse_table(recp))                                .or.  &
                             ( ( .not. (rec_pop .and. same_lu .and. same_age) )            &
                              .and. dont_force_fuse ) ) then
                           cycle recloopa
                        end if
                        !------------------------------------------------------------------!


                        !------------------------------------------------------------------!
                        !     Compare the light profile between two patches.               !
                        !------------------------------------------------------------------!
                        llevel_diff_max = 0.0
                        llevel_diff_avg = 0.0
                        npop            = 0
                        hgt_diff_max    = -999.0
                        dclai_diff_max  = -999.0
                        rclai_diff_max  = -999.0
                        dllev_diff_max  = -999.0
                        rllev_diff_max  = -999.0
                        hgtloopa: do ihgt=1,ff_nhgt
                           cumlai_recp = sum(csite%cumlai_profile(:,ihgt,recp))
                           cumlai_donp = sum(csite%cumlai_profile(:,ihgt,donp))

                           !---------------------------------------------------------------!
                           !      Exit loop in case cumulative LAI of both patches is      !
                           ! zero (above canopy at both sites).                            !
                           !---------------------------------------------------------------!
                           if ( (cumlai_recp + cumlai_donp) == 0) exit hgtloopa
                           !---------------------------------------------------------------!


                           !---------------------------------------------------------------!
                           !    Find the absolute difference in the light levels.          !
                           !---------------------------------------------------------------!
                           lnexp_donp  = max( lnexp_min                                    &
                                            , min(lnexp_max,pat_light_ext * cumlai_donp) )
                           lnexp_recp  = max( lnexp_min                                    &
                                            , min(lnexp_max,pat_light_ext * cumlai_recp) )
                           llevel_donp = exp(- lnexp_donp)
                           llevel_recp = exp(- lnexp_recp)
                           
                           llevel_diff = abs(llevel_donp - llevel_recp )
                           !---------------------------------------------------------------!

                           !---------------------------------------------------------------!
                           !    Update mean light level.                                   !
                           !---------------------------------------------------------------!
                           npop            = npop + 1 
                           llevel_diff_avg = ( llevel_diff_avg * real(npop-1)              &
                                             + llevel_diff ) / real(npop)
                           !---------------------------------------------------------------!



                           !---------------------------------------------------------------!
                           !    Update maximum light level difference.                     !
                           !---------------------------------------------------------------!
                           if (llevel_diff > llevel_diff_max) then
                              llevel_diff_max = llevel_diff
                              hgt_diff_max    = hgt_class(ihgt)
                              dclai_diff_max  = cumlai_donp
                              rclai_diff_max  = cumlai_recp
                              dllev_diff_max  = llevel_donp
                              rllev_diff_max  = llevel_recp
                           end if
                           !---------------------------------------------------------------!
                        end do hgtloopa
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Patch fusion criterion.  The average profile difference must !
                        ! be within (stricter) tolerance, and maximum difference must be   !
                        ! within (less strict) tolerance.                                  !
                        !------------------------------------------------------------------!
                        fuse_flag       = llevel_diff_avg <= pat_light_tol .and.           &
                                          llevel_diff_max <= pat_light_mxd
                        !------------------------------------------------------------------!


                        !------------------------------------------------------------------!
                        !     If fuse_flag is false, the patches aren't similar, move      !
                        ! to the next donor patch.                                         !
                        !------------------------------------------------------------------!
                        if (dont_force_fuse .and. (.not. fuse_flag)) cycle recloopa
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !    Given the potentially sheer number of combinations, we only   !
                        ! show details for successful cases.                               !
                        !------------------------------------------------------------------!
                        if (print_fuse_details) then
                           open  (unit=72,file=trim(fuse_fout),status='old',action='write' &
                                 ,position='append')
                           write (unit=72                                                  &
                                 ,fmt='(3(a,1x,i6,1x),7(a,1x,f9.3,1x),1(a,1x,l1,1x))')     &
                              '       * DONP =',donp,'RECP = ',recp,'NPOP = ',npop         &
                             ,'HGT_DIFF_MAX = ',hgt_diff_max                               &
                             ,'DCLAI_DIFF_MAX = ',dclai_diff_max                           &
                             ,'RCLAI_DIFF_MAX = ',rclai_diff_max                           &
                             ,'DLLEV_DIFF_MAX = ',dllev_diff_max                           &
                             ,'RLLEV_DIFF_MAX = ',rllev_diff_max                           &
                             ,'LLEVEL_DIFF_MAX = ',llevel_diff_max                         &
                             ,'LLEVEL_DIFF_AVG = ',llevel_diff_avg                         &
                             ,'FUSE_FLAG = ',fuse_flag
                           close (unit=72,status='keep')
                        end if
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Take an average of the patch properties of donpatch and      !
                        ! recpatch, and assign the average recpatch.                       !
                        !------------------------------------------------------------------!
                        call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)         &
                                           ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)         &
                                           ,cpoly%green_leaf_factor(:,isi),fuse_initial    &
                                           ,elim_nplant,elim_lai)
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
                        !------------------------------------------------------------------!
                     end do recloopa
                     !---------------------------------------------------------------------!
                  end do donloopa
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !      Check how many patches are valid and the original area that will  !
                  ! remain after patch termination.  We leave the loop when we have a      !
                  ! sufficient number of patches that will remain, and a sufficient area   !
                  ! that will be still represented.  Normally the second criterion is met, !
                  ! but we must check for when the initial conditions have a sheer number  !
                  ! of patches.                                                            !
                  !------------------------------------------------------------------------!
                  remain(:)       = fuse_table(:) .and. csite%area(:) >= min_patch_area
                  npatches_remain = count(remain)
                  area_remain     = sum(csite%area,mask=remain)
                  if ( npatches_remain <= abs(maxpatch) .and.                              &
                       area_remain     >= pat_min_area_remain ) exit mainfuseloopa
                  !------------------------------------------------------------------------!



                  !----- Increment tolerance ----------------------------------------------!
                  pat_light_tol = pat_light_tol * pat_light_tol_mult
                  pat_light_mxd = pat_light_tol * pat_light_mxd_fac
                  !------------------------------------------------------------------------!
               end do mainfuseloopa
               !---------------------------------------------------------------------------!
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
            pat_light_tol = pat_light_tol_min
            pat_light_mxd = pat_light_tol * pat_light_mxd_fac
            !------------------------------------------------------------------------------!

            mainfuseloop: do ifus=1,niter_patfus

               npatches_old = count(fuse_table)
               npatches_new = npatches_old

               !---------------------------------------------------------------------------!
               !    Inform that the upcoming fusions are going to be with populated        !
               ! patches, and record the tolerance used.                                   !
               !---------------------------------------------------------------------------!
               if (print_fuse_details) then
                  open (unit=72,file=trim(fuse_fout),status='old',action='write'           &
                                                     ,position='append')
                  write(unit=72,fmt='(a,1x,1(a,1x,es9.2,1x))')                             &
                                              '   + Looking for similar populated patches' &
                                             ,' - Rel. Tolerance  =',pat_light_tol
                  close(unit=72,status='keep')
               end if
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !     Loop from youngest to the second oldest patch.                        !
               !---------------------------------------------------------------------------!
               donloopp: do donp = csite%npatches,2,-1
                  donpatch => csite%patch(donp)
                  don_lu  = csite%dist_type(donp)
                  don_old = csite%age(donp) >= min_oldgrowth(don_lu)
                  don_pop = donpatch%ncohorts > 0 

                  !------------------------------------------------------------------------!
                  !     If this is an empty patch, or has already been merged, we skip it. !
                  !------------------------------------------------------------------------!
                  if ((.not. fuse_table(donp))) then
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
                     recpatch => csite%patch(recp)
                     rec_lu   = csite%dist_type(recp)
                     rec_old  = csite%age(recp) >= min_oldgrowth(rec_lu)
                     rec_pop  = recpatch%ncohorts > 0 

                     old_or_same_lu = don_lu == rec_lu .or. (don_old .and. rec_old)

                     recp_found     = old_or_same_lu .and. fuse_table(recp) .and.          &
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



                  !------------------------------------------------------------------------!
                  !     This should never happen because we have already fused all empty   !
                  ! patches, but, just in case... If both patches are empty they cannot be !
                  ! fused in this loop.                                                    !
                  !------------------------------------------------------------------------!
                  if ( dont_force_fuse .and. (.not. don_pop) .and. (.not. rec_pop)) then
                     cycle donloopp
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Compare the light profile between two patches.                     !
                  !------------------------------------------------------------------------!
                  llevel_diff_max = 0.0
                  llevel_diff_avg = 0.0
                  npop            = 0
                  hgt_diff_max    = -999.0
                  dclai_diff_max  = -999.0
                  rclai_diff_max  = -999.0
                  dllev_diff_max  = -999.0
                  rllev_diff_max  = -999.0
                  hgtloop: do ihgt=1,ff_nhgt

                     cumlai_recp = sum(csite%cumlai_profile(:,ihgt,recp))
                     cumlai_donp = sum(csite%cumlai_profile(:,ihgt,donp))

                     !---------------------------------------------------------------------!
                     !      Exit loop in case cumulative LAI of both patches is zero.      !
                     ! (above canopy at both sites).                                       !
                     !---------------------------------------------------------------------!
                     if ( (cumlai_recp + cumlai_donp) == 0) exit hgtloop
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    Find the absolute difference in the light levels.                !
                     !---------------------------------------------------------------------!
                     lnexp_donp  = max( lnexp_min                                          &
                                      , min(lnexp_max,pat_light_ext * cumlai_donp) )
                     lnexp_recp  = max( lnexp_min                                          &
                                      , min(lnexp_max,pat_light_ext * cumlai_recp) )
                     llevel_donp = exp(- lnexp_donp)
                     llevel_recp = exp(- lnexp_recp)
                     
                     llevel_diff = abs(llevel_donp - llevel_recp )
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    Update mean light level.                                         !
                     !---------------------------------------------------------------------!
                     npop            = npop + 1
                     llevel_diff_avg = ( llevel_diff_avg * real(npop-1) + llevel_diff )    &
                                     / real(npop)
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    In case we are printing the fusion steps, update maximum light   !
                     ! level difference.                                                   !
                     !---------------------------------------------------------------------!
                     if (llevel_diff > llevel_diff_max) then
                        llevel_diff_max = llevel_diff
                        hgt_diff_max    = hgt_class(ihgt)
                        dclai_diff_max  = cumlai_donp
                        rclai_diff_max  = cumlai_recp
                        dllev_diff_max  = llevel_donp
                        rllev_diff_max  = llevel_recp
                     end if
                     !---------------------------------------------------------------------!
                  end do hgtloop
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Patch fusion criterion.  The average profile difference must be    !
                  ! within (stricter) tolerance, and maximum difference must be within     !
                  ! (less strict) tolerance.                                               !
                  !------------------------------------------------------------------------!
                  fuse_flag       = llevel_diff_avg <= pat_light_tol .and.                 &
                                    llevel_diff_max <= pat_light_mxd
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     If fuse_flag is false, the patches aren't similar, move to the     !
                  ! next donor patch.                                                      !
                  !------------------------------------------------------------------------!
                  if (dont_force_fuse .and. (.not. fuse_flag)) cycle donloopp
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Given the potentially sheer number of combinations, we only         !
                  ! show details for successful cases.                                     !
                  !------------------------------------------------------------------------!
                  if (print_fuse_details) then
                     open  (unit=72,file=trim(fuse_fout),status='old',action='write'       &
                           ,position='append')
                     write (unit=72                                                        &
                           ,fmt='(3(a,1x,i6,1x),7(a,1x,f9.3,1x),1(a,1x,l1,1x))')           &
                        '       * DONP =',donp,'RECP = ',recp,'NPOP = ',npop               &
                       ,'HGT_DIFF_MAX = ',hgt_diff_max                                     &
                       ,'DCLAI_DIFF_MAX = ',dclai_diff_max                                 &
                       ,'RCLAI_DIFF_MAX = ',rclai_diff_max                                 &
                       ,'DLLEV_DIFF_MAX = ',dllev_diff_max                                 &
                       ,'RLLEV_DIFF_MAX = ',rllev_diff_max                                 &
                       ,'LLEVEL_DIFF_MAX = ',llevel_diff_max                               &
                       ,'LLEVEL_DIFF_AVG = ',llevel_diff_avg                               &
                       ,'FUSE_FLAG = ',fuse_flag
                     close (unit=72,status='keep')
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Reaching this point means that the patches are sufficiently       !
                  ! similar so they will be fused.   We take the average of the patch      !
                  ! properties of donpatch and recpatch, and leave the averaged values at  !
                  ! recpatch.                                                              !
                  !------------------------------------------------------------------------!
                  call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)               &
                                     ,cpoly%lsl(isi),cpoly%ntext_soil(:,isi)               &
                                     ,cpoly%green_leaf_factor(:,isi),fuse_initial          &
                                     ,elim_nplant,elim_lai)
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
               !      Check how many patches are valid and the original area that will     !
               ! remain after patch termination.  We leave the loop when we have a         !
               ! sufficient number of patches that will remain, and a sufficient area      !
               ! that will be still represented.  Normally the second criterion is met,    !
               ! but we must check for when the initial conditions have a sheer number     !
               ! of patches.                                                               !
               !---------------------------------------------------------------------------!
               remain(:)       = fuse_table(:) .and. csite%area(:) >= min_patch_area
               npatches_remain = count(remain)
               area_remain     = sum(csite%area,mask=remain)
               if ( npatches_remain <= abs(maxpatch) .and.                                 &
                    area_remain     >= pat_min_area_remain ) exit mainfuseloop
               !---------------------------------------------------------------------------!



               !----- Increment tolerance -------------------------------------------------!
               pat_light_tol = pat_light_tol * pat_light_tol_mult
               pat_light_mxd = pat_light_tol * pat_light_mxd_fac
               !---------------------------------------------------------------------------!
            end do mainfuseloop
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      There is a chance that the fusion process will still allow a            !
            ! significant information loss.  This could happen in case min_patch_area is   !
            ! very large and the current tolerance settings didn't merge enough patches    !
            ! that allowed most of the information to remain.  In case this happens, warn  !
            ! the user.  In the extreme case in which all patches would be gone, stop the  !
            ! run.                                                                         !
            !------------------------------------------------------------------------------!
            remain(:)       = fuse_table(:) .and. csite%area(:) >= min_patch_area
            npatches_remain = count(remain)
            area_remain     = sum(csite%area,mask=remain)
            if ( area_remain < pat_min_area_remain ) then
               write(unit=*,fmt='(a)') '--------------------------------------------------'
               if (npatches_remain > 0) then
                  write(unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING!'
                  write(unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING!'
                  write(unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING!'
                  write(unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING!'
                  write(unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING!'
                  write(unit=*,fmt='(a)') '   WARNING! WARNING! WARNING! WARNING! WARNING!'
               else
                  write(unit=*,fmt='(a)') '   PROBLEM! PROBLEM! PROBLEM! PROBLEM! PROBLEM!'
                  write(unit=*,fmt='(a)') '   PROBLEM! PROBLEM! PROBLEM! PROBLEM! PROBLEM!'
                  write(unit=*,fmt='(a)') '   PROBLEM! PROBLEM! PROBLEM! PROBLEM! PROBLEM!'
                  write(unit=*,fmt='(a)') '   PROBLEM! PROBLEM! PROBLEM! PROBLEM! PROBLEM!'
                  write(unit=*,fmt='(a)') '   PROBLEM! PROBLEM! PROBLEM! PROBLEM! PROBLEM!'
                  write(unit=*,fmt='(a)') '   PROBLEM! PROBLEM! PROBLEM! PROBLEM! PROBLEM!'
               end if
               write(unit=*,fmt='(a)') '--------------------------------------------------'
               write(unit=*,fmt='(a)') ' This simulation has too many tiny patches.'
               write(unit=*,fmt='(a)') ' Significant amount of information may be lost.'
               write(unit=*,fmt='(a)') '--------------------------------------------------'
               write(unit=*,fmt='(a,1x,i6)'  ) ' NPATCHES  (ORIG)   =',csite%npatches
               write(unit=*,fmt='(a,1x,i6)'  ) ' NPATCHES  (REMAIN) =',npatches_remain
               write(unit=*,fmt='(a,1x,f9.5)') ' SITE AREA (REMAIN) =',area_remain
               write(unit=*,fmt='(a,1x,f9.5)') ' MIN. SITE AREA OK  =',pat_min_area_remain
               write(unit=*,fmt='(a,1x,f9.5)') ' MIN. PATCH AREA    =',min_patch_area
               write(unit=*,fmt='(a,1x,f9.5)') ' MIN. LIGHT TOL.    =',pat_light_tol_min
               write(unit=*,fmt='(a,1x,f9.5)') ' MAX. LIGHT TOL.    =',pat_light_tol_max
               write(unit=*,fmt='(a,1x,f9.4)') ' DEVIATION FACTOR   =',pat_light_mxd_fac
               write(unit=*,fmt='(a)') '--------------------------------------------------'
               write(unit=*,fmt='(a)') ' Consider making min_patch_area smaller, or change'
               write(unit=*,fmt='(a)') ' patch fusion parameters.'
               write(unit=*,fmt='(a)') '--------------------------------------------------'
               
               if (npatches_remain == 0) then
                  call fatal_error('Patch fusion was going to terminate all patches'       &
                                  ,'fuse_patches','fuse_fiss_utils.f90')
               end if
            end if 
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
            deallocate(tempsite  )
            deallocate(fuse_table)
            deallocate(remain    )
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
      !     Print a banner to inform the user how many patches and cohorts exist.  Make a  !
      ! simple banner in case this is a multi-polygon simulation, or a detailed banner in  !
      ! case this is a single polygon simulation.                                          !
      !------------------------------------------------------------------------------------!
      tot_npolygons = cgrid%npolygons
      if (tot_npolygons > 1) then
         tot_ncohorts  = 0
         tot_npatches  = 0
         tot_nsites    = 0
                     npatches_remain = count(remain)
                     area_remain     = sum(csite%area,mask=remain)
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
         write (unit=*,fmt='(6(a,1x,i8,1x))')                                              &
              'After fusion, total count in node',mynum,'for grid',ifm                     &
             ,': POLYGONS=',tot_npolygons,'SITES=',tot_nsites,'PATCHES=',tot_npatches      &
             ,'COHORTS=',tot_ncohorts
      else
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         write(unit=*,fmt='(a)') ' Patch status for each site after fuse_patches.'
         write(unit=*,fmt='(a)') '   - >= MIN: patches with area >= min_patch_area.'
         write(unit=*,fmt='(a)') '   - TOTAL: all patches.'
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         ipy = 1
         cpoly => cgrid%polygon(ipy)
         do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            tot_npatches    = csite%npatches
            npatches_remain = count(csite%area >= min_patch_area)
            area_remain     = sum  (csite%area,mask=csite%area >= min_patch_area)
            tot_ncohorts    = 0
            ncohorts_remain = 0
            do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)
               tot_ncohorts = tot_ncohorts + cpatch%ncohorts
               if (csite%area(ipa) >= min_patch_area) then
                  ncohorts_remain = ncohorts_remain + cpatch%ncohorts
               end if
            end do
            
            write (unit=*,fmt='(a,i5,a,2(a,i5),2(a,f9.4),2(a,i5))') ' Site ',isi,' -- '    &
               ,'    PATCHES (>= MIN/TOTAL) = ',npatches_remain,'/',tot_npatches           &
               ,'    AREA (>= MIN/TOTAL) = '   ,area_remain,'/',1.0                        &
               ,'    COHORTS (>= MIN/TOTAL) = ',ncohorts_remain,'/',tot_ncohorts 
         end do
         write(unit=*,fmt='(a)') '--------------------------------------------------------'
         write(unit=*,fmt='(a)') ' '
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine new_fuse_patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: OLD_FUSE_PATCHES  
   !> \brief This subroutine will perform patch fusion based on some similarity
   !> criteria to determine whether they can be fused with no significant loss of
   !> information.
   !> \details The user is welcome to set up a benchmark, but they should be aware
   !> that no miracles will happen here. If there are more very distinct patches
   !> than maxpatch, then the user will need to live with that and accept life is
   !> not always fair with those with limited computational resources.
   !> \warning Cohort fusion may significantly affect the runs, so I am keeping both the 
   !> old and new routines, but I think we should probably phase out this and stick with
   !> the new one.
   !---------------------------------------------------------------------------------------!
   subroutine old_fuse_patches(cgrid,ifm,fuse_initial)
      use update_derived_utils, only : patch_pft_size_profile ! ! subroutine
      use ed_state_vars       , only : edtype                 & ! structure
                                     , polygontype            & ! structure
                                     , sitetype               & ! structure
                                     , patchtype              ! ! structure
      use fusion_fission_coms , only : ff_nhgt                & ! intent(in)
                                     , niter_patfus           & ! intent(in)
                                     , dark_cumlai_max        & ! intent(in)
                                     , dark_cumlai_mult       & ! intent(in)
                                     , sunny_cumlai_min       & ! intent(in)
                                     , sunny_cumlai_mult      & ! intent(in)
                                     , print_fuse_details     & ! intent(in)
                                     , light_toler_min        & ! intent(in)
                                     , light_toler_mult       & ! intent(in)
                                     , fuse_prefix            ! ! intent(in)
      use disturb_coms        , only : min_oldgrowth          ! ! intent(in)
      use ed_max_dims         , only : n_pft                  & ! intent(in)
                                     , str_len                ! ! intent(in)
      use mem_polygons        , only : maxpatch               ! ! intent(in)
      use ed_node_coms        , only : mynum                  ! ! intent(in)
      use ed_misc_coms        , only : current_time           ! ! intent(in)
      use grid_coms           , only : nzg                    & ! intent(in)
                                     , nzs                    ! ! intent(in)
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
                  call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)               &
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
                        call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)         &
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
                  call fuse_2_patches(csite,donp,recp,nzg,nzs,cpoly%met(isi)               &
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
   end subroutine old_fuse_patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: FUSE_2_PATCHES  
   !> \brief This subroutine will merge two patches into 1.
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_patches(csite,donp,recp,mzg,mzs,cmet,lsl,ntext_soil,green_leaf_factor &
                            ,fuse_initial,elim_nplant,elim_lai)
      use ed_state_vars       , only : sitetype                      & ! Structure 
                                     , patchtype                     & ! Structure
                                     , allocate_sitetype             & ! sub-routine
                                     , copy_sitetype                 & ! sub-routine
                                     , deallocate_sitetype           ! ! sub-routine
      use met_driver_coms     , only : met_driv_state                ! ! structure
      use soil_coms           , only : soil                          & ! intent(in)
                                     , tiny_sfcwater_mass            & ! intent(in)
                                     , matric_potential              ! ! intent(in)
      use ed_max_dims         , only : n_pft                         & ! intent(in)
                                     , n_dbh                         ! ! intent(in)
      use mem_polygons        , only : maxcohort                     ! ! intent(in)
      use therm_lib           , only : uextcm2tl                     & ! subroutine
                                     , uint2tl                       & ! subroutine
                                     , idealdenssh                   & ! function
                                     , idealdmolsh                   & ! function
                                     , press2exner                   & ! function
                                     , extemp2theta                  & ! function
                                     , extheta2temp                  & ! function
                                     , tq2enthalpy                   & ! function
                                     , hq2temp                       & ! function
                                     , virtt                         & ! function
                                     , thetaeiv                      & ! function
                                     , vpdefil                       ! ! function
      use ed_misc_coms        , only : writing_long                  & ! intent(in)
                                     , writing_eorq                  & ! intent(in)
                                     , writing_dcyc                  & ! intent(in)
                                     , ndcycle                       & ! intent(in)
                                     , frqsum                        & ! intent(in)
                                     , frqsumi                       ! ! intent(in)
      use consts_coms         , only : wdns                          & ! intent(in)
                                     , rdry                          ! ! intent(in)
      use fusion_fission_coms , only : ifusion                       & ! intent(in)
                                     , corr_patch                    ! ! intent(in)
      use rk4_coms            , only : checkbudget                   ! ! intent(in)
      use ed_type_init        , only : new_patch_sfc_props           ! ! sub-routine
      use update_derived_utils, only : update_patch_derived_props    & ! sub-routine
                                     , update_cohort_extensive_props & ! sub-routine
                                     , patch_pft_size_profile        ! ! sub-routine
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)         , target      :: csite             ! Current site
      integer                , intent(in)  :: donp              ! Donating patch
      integer                , intent(in)  :: recp              ! Receptor patch
      integer                , intent(in)  :: lsl               ! Lowest soil level
      integer                , intent(in)  :: mzg               ! # of soil layers
      integer                , intent(in)  :: mzs               ! # of sfc. water layers
      type(met_driv_state)   , target      :: cmet              ! Current met forcing
      integer, dimension(mzg), intent(in)  :: ntext_soil        ! Soil type
      real, dimension(n_pft) , intent(in)  :: green_leaf_factor ! Green leaf factor...
      logical                , intent(in)  :: fuse_initial      ! Initialisation?
      real                   , intent(out) :: elim_nplant       ! Eliminated nplant 
      real                   , intent(out) :: elim_lai          ! Eliminated lai
      !----- Local variables --------------------------------------------------------------!
      type(sitetype)         , pointer     :: osite             ! Orig. patches bef. fusion
      type(patchtype)        , pointer     :: cpatch            ! Current patch
      type(patchtype)        , pointer     :: temppatch         ! Temporary patch
      integer                              :: iii               ! Counters
      integer                              :: nsoil             ! Alias for soil texture
      integer                              :: t                 ! Counter for time of day
      integer                              :: ndc               ! # of cohorts - donp patch
      integer                              :: nrc               ! # of cohorts - recp patch
      real                                 :: xmean_can_exner   ! Exner function - CAS
      real                                 :: newarea           ! new patch area
      real                                 :: rawgt             ! Weight for receptor patch
      real                                 :: dawgt             ! Weight for donor patch
      !----- The following variables are for conserving canopy air space. -----------------!
      real  :: can_r_depth0_donp       !< Old can. depth, mass  (donor   ) [          J/kg]
      real  :: can_r_depth0_recp       !< Old can. depth, mass  (receptor) [          J/kg]
      real  :: can_d_depth0_donp       !< Old can. depth, mass  (donor   ) [          J/kg]
      real  :: can_d_depth0_recp       !< Old can. depth, mass  (receptor) [          J/kg]
      real  :: can_enthalpy_donp       !< Specific enthalpy     (donor   ) [          J/kg]
      real  :: can_enthalpy_recp       !< Specific enthalpy     (receptor) [          J/kg]
      real  :: can_rvap_recp           !< Water mixing ratio    (receptor) [         kg/kg]
      real  :: can_exner_recp          !< Exner function        (receptor) [        J/kg/K]
      real  :: cb_enthalpy_donp        !< Total enthalpy        (donor   ) [          J/m2]
      real  :: cb_enthalpy_recp        !< Total enthalpy        (receptor) [          J/m2]
      real  :: cb_mass_donp            !< Total air mass        (donor   ) [     kg_air/m2]
      real  :: cb_mass_recp            !< Total air mass        (receptor) [     kg_air/m2]
      real  :: cb_molar_donp           !< Total dry molar count (donor   ) [mol_dry_air/m2]
      real  :: cb_molar_recp           !< Total dry molar count (receptor) [mol_dry_air/m2]
      real  :: cb_water_donp           !< Total water mass      (donor   ) [     kg_h2o/m2]
      real  :: cb_water_recp           !< Total water mass      (receptor) [     kg_h2o/m2]
      real  :: cb_co2_donp             !< Total CO2 mass        (donor   ) [     kg_co2/m2]
      real  :: cb_co2_recp             !< Total CO2 mass        (receptor) [     kg_co2/m2]
      real  :: rbudget_zcaneffect_recp !< Dz effect, mass (rho) (receptor) [   kg_air/m2/s]
      real  :: rbudget_zcaneffect_donp !< Dz effect, mass (rho) (donor   ) [   kg_air/m2/s]
      real  :: dbudget_zcaneffect_recp !< Dz effect, molar cnt  (receptor) [mol_d_air/m2/s]
      real  :: dbudget_zcaneffect_donp !< Dz effect, molar cnt  (donor   ) [mol_d_air/m2/s]
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     In case we are checking budget, make a copy of both receptor and donor sites.  !
      !------------------------------------------------------------------------------------!
      if (checkbudget .and. (.not. fuse_initial)) then
         !----- Allocate site with 2 patches (1 = receptor; 2 = donor). -------------------!
         nullify (osite)
         allocate(osite)
         call allocate_sitetype(osite,2)
         call copy_sitetype(csite,osite,recp,recp,1,1)
         call copy_sitetype(csite,osite,donp,donp,2,2)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     This function fuses the two patches specified in the argument. It fuses the    !
      ! first patch in the argument (the "donor" = donp ) into the second patch in the     !
      ! argument (the "recipient" = recp ), and frees the memory associated with the donor !
      ! patch.                                                                             !
      !------------------------------------------------------------------------------------!

      !----- The new area is simply the sum of each patch area. ---------------------------!
      newarea  = csite%area(donp) + csite%area(recp)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Find the scaling factor based on area.  This should be more robust for when the !
      ! areas of the patches to be fused are very different.                               !
      !------------------------------------------------------------------------------------!
      dawgt   = csite%area(donp) / ( csite%area(donp) + csite%area(recp) )
      rawgt   = 1.0 - dawgt
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
      csite%age(recp)                = csite%age(donp)                * dawgt              &
                                     + csite%age(recp)                * rawgt

      csite%htry(recp)               = csite%htry(donp)               * dawgt              &
                                     + csite%htry(recp)               * rawgt

      csite%hprev(recp)              = csite%hprev(donp)              * dawgt              &
                                     + csite%hprev(recp)              * rawgt

      csite%fbeam(recp)              = csite%fbeam(donp)              * dawgt              &
                                     + csite%fbeam(recp)              * rawgt

      csite%fast_grnd_C(recp)        = csite%fast_grnd_C(donp)        * dawgt              &
                                     + csite%fast_grnd_C(recp)        * rawgt

      csite%fast_soil_C(recp)        = csite%fast_soil_C(donp)        * dawgt              &
                                     + csite%fast_soil_C(recp)        * rawgt

      csite%structural_grnd_C(recp)  = csite%structural_grnd_C(donp)  * dawgt              &
                                     + csite%structural_grnd_C(recp)  * rawgt

      csite%structural_soil_C(recp)  = csite%structural_soil_C(donp)  * dawgt              &
                                     + csite%structural_soil_C(recp)  * rawgt

      csite%structural_grnd_L(recp)  = csite%structural_grnd_L(donp)  * dawgt              &
                                     + csite%structural_grnd_L(recp)  * rawgt
                                     

      csite%structural_soil_L(recp)  = csite%structural_soil_L(donp)  * dawgt              &
                                     + csite%structural_soil_L(recp)  * rawgt
                                     

      csite%microbial_soil_C(recp)   = csite%microbial_soil_C(donp)   * dawgt              &
                                     + csite%microbial_soil_C(recp)   * rawgt

      csite%slow_soil_C(recp)        = csite%slow_soil_C(donp)        * dawgt              &
                                     + csite%slow_soil_C(recp)        * rawgt

      csite%passive_soil_C(recp)     = csite%passive_soil_C(donp)     * dawgt              &
                                     + csite%passive_soil_C(recp)     * rawgt

      csite%fast_grnd_N(recp)        = csite%fast_grnd_N(donp)        * dawgt              &
                                     + csite%fast_grnd_N(recp)        * rawgt

      csite%fast_soil_N(recp)        = csite%fast_soil_N(donp)        * dawgt              &
                                     + csite%fast_soil_N(recp)        * rawgt

      csite%structural_grnd_N(recp)  = csite%structural_grnd_N(donp)  * dawgt              &
                                     + csite%structural_grnd_N(recp)  * rawgt

      csite%structural_soil_N(recp)  = csite%structural_soil_N(donp)  * dawgt              &
                                     + csite%structural_soil_N(recp)  * rawgt

      csite%mineralized_soil_N(recp) = csite%mineralized_soil_N(donp) * dawgt              &
                                     + csite%mineralized_soil_N(recp) * rawgt

      csite%sum_dgd(recp)            = csite%sum_dgd(donp)            * dawgt              &
                                     + csite%sum_dgd(recp)            * rawgt

      csite%sum_chd(recp)            = csite%sum_chd(donp)            * dawgt              &
                                     + csite%sum_chd(recp)            * rawgt

      csite%ggbare(recp)             = csite%ggbare(donp)             * dawgt              &
                                     + csite%ggbare(recp)             * rawgt

      csite%ggnet(recp)              = csite%ggnet(donp)              * dawgt              &
                                     + csite%ggnet(recp)              * rawgt

      csite%ggsoil(recp)             = csite%ggsoil(donp)             * dawgt              &
                                     + csite%ggsoil(recp)             * rawgt

      csite%ustar (recp)             = csite%ustar (donp)             * dawgt              &
                                     + csite%ustar (recp)             * rawgt

      csite%tstar (recp)             = csite%tstar (donp)             * dawgt              &
                                     + csite%tstar (recp)             * rawgt

      csite%qstar (recp)             = csite%qstar (donp)             * dawgt              &
                                     + csite%qstar (recp)             * rawgt

      csite%cstar (recp)             = csite%cstar (donp)             * dawgt              &
                                     + csite%cstar (recp)             * rawgt
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Fusion of canopy air space must ensure mass and enthalpy conservation.        !
      ! Total mass [kg/m2] is rho*z.  We assume that the canopy depth of the fused patch   !
      ! is going to be the weighted average (in case it is slightly different, we correct  !
      ! it later in the subroutine).  Importantly, we must use the canopy depth from       !
      ! "before" the vegetation dynamics if we want to track conservation, and the         !
      ! zcaneffect variables have not yet been applied to storage, and also determine the  !
      ! zcaneffect on total mass ("mbudget_zcaneffect") and total dry-air molar count      !
      ! ("dbudget_zcaneffect"), which will be used after fusion.  In the case we are not   !
      ! checking the budget, we cannot determine the previous height so we use the current !
      ! depth instead.  This may lead to slightly different results if running with budget !
      ! check or not.                                                                      !
      !                                                                                    !
      !    For the time being this is only applied to the instantaneous (state) variables  !
      ! because these are the ones that are used in the budget assessment sub-routines.    !
      ! The ?mean variables are for diagnostics only; even though they could go with the   !
      ! same approach, we use the simpler fusion for the time being.                       !
      !------------------------------------------------------------------------------------!
      if (checkbudget) then
         !------ Previous canopy depth (mass). Use water budget to define it. -------------!
         can_r_depth0_donp       = csite%can_depth(donp)                                   &
                                 - frqsum * csite%wbudget_zcaneffect(donp)                 &
                                 / ( csite%can_rhos(donp) * csite%can_shv(donp) )
         can_r_depth0_recp       = csite%can_depth(recp)                                   &
                                 - frqsum * csite%wbudget_zcaneffect(recp)                 &
                                 / ( csite%can_rhos(recp) * csite%can_shv(recp) )
         !------ Previous canopy depth (molar). Use CO2 budget to define it. --------------!
         can_d_depth0_donp       = csite%can_depth(donp)                                   &
                                 - frqsum * csite%co2budget_zcaneffect(donp)               &
                                 / ( csite%can_dmol(donp) * csite%can_co2(donp) )
         can_d_depth0_recp       = csite%can_depth(recp)                                   &
                                 - frqsum * csite%co2budget_zcaneffect(recp)               &
                                 / ( csite%can_dmol(recp) * csite%can_co2(recp) )
         !------ Canopy change effect on total mass. --------------------------------------!
         rbudget_zcaneffect_donp = csite%can_rhos(donp)                                    &
                                 * (csite%can_depth(donp) - can_r_depth0_donp) * frqsumi
         rbudget_zcaneffect_recp = csite%can_rhos(recp)                                    &
                                 * (csite%can_depth(recp) - can_r_depth0_recp) * frqsumi
         !------ Canopy change effect on dry-air molar count. -----------------------------!
         dbudget_zcaneffect_donp = csite%can_dmol(donp)                                    &
                                 * (csite%can_depth(donp) - can_d_depth0_donp) * frqsumi
         dbudget_zcaneffect_recp = csite%can_dmol(recp)                                    &
                                 * (csite%can_depth(recp) - can_d_depth0_recp) * frqsumi
         !---------------------------------------------------------------------------------!
      else
         !------ Previous depth can't be determined. --------------------------------------!
         can_r_depth0_donp       = csite%can_depth(donp)
         can_r_depth0_recp       = csite%can_depth(recp)
         can_d_depth0_donp       = csite%can_depth(donp)
         can_d_depth0_recp       = csite%can_depth(recp)
         rbudget_zcaneffect_donp = 0.0
         rbudget_zcaneffect_recp = 0.0
         dbudget_zcaneffect_donp = 0.0
         dbudget_zcaneffect_recp = 0.0
         !---------------------------------------------------------------------------------!
      end if
      !------ Find the specific enthalpy of receptor and donor patch [J/kg]. --------------!
      can_enthalpy_donp = tq2enthalpy(csite%can_temp(donp),csite%can_shv (donp),.true.)
      can_enthalpy_recp = tq2enthalpy(csite%can_temp(recp),csite%can_shv (recp),.true.)
      !------ Find the total canopy air space mass [kg_air/m2]. ---------------------------!
      cb_mass_donp      = csite%can_rhos(donp) * can_r_depth0_donp
      cb_mass_recp      = csite%can_rhos(recp) * can_r_depth0_recp
      !------ Find the molar count of dry air in the canopy air space [mol_dry_air/m2]. ---!
      cb_molar_donp     = csite%can_dmol(donp) * can_d_depth0_donp
      cb_molar_recp     = csite%can_dmol(recp) * can_d_depth0_recp
      !------ Find the bulk enthalpy of receptor and donor patch [J/m2]. ------------------!
      cb_enthalpy_donp  = cb_mass_donp * can_enthalpy_donp
      cb_enthalpy_recp  = cb_mass_recp * can_enthalpy_recp
      !------ Find the total water mass [kg_h2o/m2]. --------------------------------------!
      cb_water_donp     = cb_mass_donp * csite%can_shv(donp)
      cb_water_recp     = cb_mass_recp * csite%can_shv(recp)
      !------ Find the total CO2 mass [kg_co2/m2]. ----------------------------------------!
      cb_co2_donp       = cb_molar_donp * csite%can_co2(donp)
      cb_co2_recp       = cb_molar_recp * csite%can_co2(recp)
      !------ Find the total properties (X/m2) of the fused patch. ------------------------!
      cb_enthalpy_recp  = cb_enthalpy_donp * dawgt + cb_enthalpy_recp * rawgt
      cb_mass_recp      = cb_mass_donp     * dawgt + cb_mass_recp     * rawgt
      cb_molar_recp     = cb_molar_donp    * dawgt + cb_molar_recp    * rawgt
      cb_water_recp     = cb_water_donp    * dawgt + cb_water_recp    * rawgt
      cb_co2_recp       = cb_co2_donp      * dawgt + cb_co2_recp      * rawgt
      !------------------------------------------------------------------------------------!





 


      !------------------------------------------------------------------------------------!
      !     Budget variables.                                                              !
      !------------------------------------------------------------------------------------!
      csite%co2budget_initialstorage(recp) = csite%co2budget_initialstorage(donp) * dawgt  &
                                           + csite%co2budget_initialstorage(recp) * rawgt
      csite%co2budget_residual      (recp) = csite%co2budget_residual      (donp) * dawgt  &
                                           + csite%co2budget_residual      (recp) * rawgt
      csite%co2budget_loss2atm      (recp) = csite%co2budget_loss2atm      (donp) * dawgt  &
                                           + csite%co2budget_loss2atm      (recp) * rawgt
      csite%co2budget_denseffect    (recp) = csite%co2budget_denseffect    (donp) * dawgt  &
                                           + csite%co2budget_denseffect    (recp) * rawgt
      csite%co2budget_zcaneffect    (recp) = csite%co2budget_zcaneffect    (donp) * dawgt  &
                                           + csite%co2budget_zcaneffect    (recp) * rawgt
      csite%co2budget_gpp           (recp) = csite%co2budget_gpp           (donp) * dawgt  &
                                           + csite%co2budget_gpp           (recp) * rawgt
      csite%co2budget_plresp        (recp) = csite%co2budget_plresp        (donp) * dawgt  &
                                           + csite%co2budget_plresp        (recp) * rawgt
      csite%co2budget_rh            (recp) = csite%co2budget_rh            (donp) * dawgt  &
                                           + csite%co2budget_rh            (recp) * rawgt
      csite%cbudget_initialstorage  (recp) = csite%cbudget_initialstorage  (donp) * dawgt  &
                                           + csite%cbudget_initialstorage  (recp) * rawgt
      csite%cbudget_residual        (recp) = csite%cbudget_residual        (donp) * dawgt  &
                                           + csite%cbudget_residual        (recp) * rawgt
      csite%cbudget_loss2atm        (recp) = csite%cbudget_loss2atm        (donp) * dawgt  &
                                           + csite%cbudget_loss2atm        (recp) * rawgt
      csite%cbudget_committed       (recp) = csite%cbudget_committed       (donp) * dawgt  &
                                           + csite%cbudget_committed       (recp) * rawgt
      csite%cbudget_denseffect      (recp) = csite%cbudget_denseffect      (donp) * dawgt  &
                                           + csite%cbudget_denseffect      (recp) * rawgt
      csite%cbudget_zcaneffect      (recp) = csite%cbudget_zcaneffect      (donp) * dawgt  &
                                           + csite%cbudget_zcaneffect      (recp) * rawgt
      csite%cbudget_loss2yield      (recp) = csite%cbudget_loss2yield      (donp) * dawgt  &
                                           + csite%cbudget_loss2yield      (recp) * rawgt
      csite%cbudget_seedrain        (recp) = csite%cbudget_seedrain        (donp) * dawgt  &
                                           + csite%cbudget_seedrain        (recp) * rawgt
      csite%ebudget_initialstorage  (recp) = csite%ebudget_initialstorage  (donp) * dawgt  &
                                           + csite%ebudget_initialstorage  (recp) * rawgt
      csite%ebudget_residual        (recp) = csite%ebudget_residual        (donp) * dawgt  &
                                           + csite%ebudget_residual        (recp) * rawgt
      csite%ebudget_netrad          (recp) = csite%ebudget_netrad          (donp) * dawgt  &
                                           + csite%ebudget_netrad          (recp) * rawgt
      csite%ebudget_loss2atm        (recp) = csite%ebudget_loss2atm        (donp) * dawgt  &
                                           + csite%ebudget_loss2atm        (recp) * rawgt
      csite%ebudget_denseffect      (recp) = csite%ebudget_denseffect      (donp) * dawgt  &
                                           + csite%ebudget_denseffect      (recp) * rawgt
      csite%ebudget_prsseffect      (recp) = csite%ebudget_prsseffect      (donp) * dawgt  &
                                           + csite%ebudget_prsseffect      (recp) * rawgt
      csite%ebudget_hcapeffect      (recp) = csite%ebudget_hcapeffect      (donp) * dawgt  &
                                           + csite%ebudget_hcapeffect      (recp) * rawgt
      csite%ebudget_wcapeffect      (recp) = csite%ebudget_wcapeffect      (donp) * dawgt  &
                                           + csite%ebudget_wcapeffect      (recp) * rawgt
      csite%ebudget_zcaneffect      (recp) = csite%ebudget_zcaneffect      (donp) * dawgt  &
                                           + csite%ebudget_zcaneffect      (recp) * rawgt
      csite%ebudget_pheneffect      (recp) = csite%ebudget_pheneffect      (donp) * dawgt  &
                                           + csite%ebudget_pheneffect      (recp) * rawgt
      csite%ebudget_loss2runoff     (recp) = csite%ebudget_loss2runoff     (donp) * dawgt  &
                                           + csite%ebudget_loss2runoff     (recp) * rawgt
      csite%ebudget_loss2drainage   (recp) = csite%ebudget_loss2drainage   (donp) * dawgt  &
                                           + csite%ebudget_loss2drainage   (recp) * rawgt
      csite%ebudget_precipgain      (recp) = csite%ebudget_precipgain      (donp) * dawgt  &
                                           + csite%ebudget_precipgain      (recp) * rawgt
      csite%wbudget_initialstorage  (recp) = csite%wbudget_initialstorage  (donp) * dawgt  &
                                           + csite%wbudget_initialstorage  (recp) * rawgt
      csite%wbudget_residual        (recp) = csite%wbudget_residual        (donp) * dawgt  &
                                           + csite%wbudget_residual        (recp) * rawgt
      csite%wbudget_loss2atm        (recp) = csite%wbudget_loss2atm        (donp) * dawgt  &
                                           + csite%wbudget_loss2atm        (recp) * rawgt
      csite%wbudget_denseffect      (recp) = csite%wbudget_denseffect      (donp) * dawgt  &
                                           + csite%wbudget_denseffect      (recp) * rawgt
      csite%wbudget_wcapeffect      (recp) = csite%wbudget_wcapeffect      (donp) * dawgt  &
                                           + csite%wbudget_wcapeffect      (recp) * rawgt
      csite%wbudget_zcaneffect      (recp) = csite%wbudget_zcaneffect      (donp) * dawgt  &
                                           + csite%wbudget_zcaneffect      (recp) * rawgt
      csite%wbudget_pheneffect      (recp) = csite%wbudget_pheneffect      (donp) * dawgt  &
                                           + csite%wbudget_pheneffect      (recp) * rawgt
      csite%wbudget_loss2runoff     (recp) = csite%wbudget_loss2runoff     (donp) * dawgt  &
                                           + csite%wbudget_loss2runoff     (recp) * rawgt
      csite%wbudget_loss2drainage   (recp) = csite%wbudget_loss2drainage   (donp) * dawgt  &
                                           + csite%wbudget_loss2drainage   (recp) * rawgt
      csite%wbudget_precipgain      (recp) = csite%wbudget_precipgain      (donp) * dawgt  &
                                           + csite%wbudget_precipgain      (recp) * rawgt
      !------ Additional budget variables for conserving canopy-air space. ----------------!
      if (checkbudget) then
         rbudget_zcaneffect_recp           = rbudget_zcaneffect_donp              * dawgt  &
                                           + rbudget_zcaneffect_recp              * rawgt
         dbudget_zcaneffect_recp           = dbudget_zcaneffect_donp              * dawgt  &
                                           + dbudget_zcaneffect_recp              * rawgt
      end if
      !------------------------------------------------------------------------------------!




      !------ Find the fused canopy depth using the regular weighting average. ------------!
      csite%can_depth  (recp)  = csite%can_depth(donp) * dawgt                             &
                               + csite%can_depth(recp) * rawgt
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Because the fused patch should already contain the zcan effect, we add the    !
      ! zcan effect back to the total masses.                                              !
      !------------------------------------------------------------------------------------!
      if (checkbudget) then
         cb_mass_recp     = cb_mass_recp     + frqsum * rbudget_zcaneffect_recp
         cb_molar_recp    = cb_molar_recp    + frqsum * dbudget_zcaneffect_recp
         cb_co2_recp      = cb_co2_recp      + frqsum * csite%co2budget_zcaneffect(recp)
         cb_water_recp    = cb_water_recp    + frqsum * csite%wbudget_zcaneffect  (recp)
         cb_enthalpy_recp = cb_enthalpy_recp + frqsum * csite%ebudget_zcaneffect  (recp)
      end if
      !------------------------------------------------------------------------------------!



      !------ Find air density and the specific properties. -------------------------------!
      csite%can_rhos   (recp)  = cb_mass_recp     / csite%can_depth(recp)
      csite%can_dmol   (recp)  = cb_molar_recp    / csite%can_depth(recp)
      can_enthalpy_recp        = cb_enthalpy_recp / cb_mass_recp
      csite%can_shv    (recp)  = cb_water_recp    / cb_mass_recp
      csite%can_co2    (recp)  = cb_co2_recp      / cb_molar_recp
      !------ Water mixing ratio (used by can_theiv). -------------------------------------!
      can_rvap_recp            = csite%can_shv(recp) / (1.0 - csite%can_shv(recp))
      !------ Find the pressure that is consistent with the ideal gas law. ----------------!
      csite%can_prss   (recp)  = rdry * csite%can_rhos(recp)                               &
                               * virtt(csite%can_temp(recp),can_rvap_recp)
      can_exner_recp           = press2exner(csite%can_prss(recp))
      !------ Find temperature from enthalpy, then find potential temperature. ------------!
      csite%can_temp    (recp) = hq2temp(can_enthalpy_recp,csite%can_shv(recp),.true.)
      csite%can_theta   (recp) = extemp2theta(can_exner_recp,csite%can_temp(recp))
      csite%can_theiv   (recp) = thetaeiv( csite%can_theta(recp), csite%can_prss (recp)    &
                                         , csite%can_temp (recp), can_rvap_recp            &
                                         , can_rvap_recp        )
      !------ Update vapour pressure deficit. ---------------------------------------------!
      csite%can_vpdef   (recp) = vpdefil( csite%can_prss(recp), csite%can_temp(recp)       &
                                        , csite%can_shv (recp), .true.               )
      !------------------------------------------------------------------------------------!
      !      Previous temperature likely needs to be done differently for a truly          !
      ! conserving model, however it does not directly affect the current state.  Leaving  !
      ! the simple approach for now.                                                       !
      !------------------------------------------------------------------------------------!
      csite%can_temp_pv (recp) = csite%can_temp_pv(donp) * dawgt                           &
                               + csite%can_temp_pv(recp) * rawgt
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Committed respiration.  These variables are used to release storage and growth !
      ! respiration to the canopy air space.  Cohort variables cannot be used because      !
      ! these are committed emissions that should be honoured even if the cohort is        !
      ! terminated, to ensure carbon conservation.                                         !
      !------------------------------------------------------------------------------------!
      csite%commit_storage_resp    (recp)  = csite%commit_storage_resp     (donp) * dawgt  &
                                           + csite%commit_storage_resp     (recp) * rawgt
      csite%commit_growth_resp     (recp)  = csite%commit_growth_resp      (donp) * dawgt  &
                                           + csite%commit_growth_resp      (recp) * rawgt
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !    There is no guarantee that there will be a minimum amount of mass in the tempo- !
      ! rary layer, nor is there any reason for both patches to have the same number of    !
      ! layers. In order to be safe, the fusion must happen in 5 stages.                   !
      !------------------------------------------------------------------------------------!
      !----- 1. Find the "extensive" sfcwater_energy (convert from J/kg to J/m2); ---------!
      do iii=1,csite%nlev_sfcwater(donp)
         csite%sfcwater_energy(iii,donp) = csite%sfcwater_energy(iii,donp)                 &
                                         * csite%sfcwater_mass  (iii,donp)
      end do
      do iii=1,csite%nlev_sfcwater(recp)
         csite%sfcwater_energy(iii,recp) = csite%sfcwater_energy(iii,recp)                 &
                                         * csite%sfcwater_mass  (iii,recp)
      end do
      !------------------------------------------------------------------------------------!
      ! 2. Merge all layers into one.  If needed, the layer will be split again next time  !
      !    the Runge-Kutta integrator is called.  After adding the value to the first      !
      !    layer, discard the value.                                                       !
      !------------------------------------------------------------------------------------!
      do iii=2,csite%nlev_sfcwater(donp)
         csite%sfcwater_energy(  1,donp) = csite%sfcwater_energy(  1,donp)                 &
                                         + csite%sfcwater_energy(iii,donp)
         csite%sfcwater_depth (  1,donp) = csite%sfcwater_depth (  1,donp)                 &
                                         + csite%sfcwater_depth (iii,donp)
         csite%sfcwater_mass  (  1,donp) = csite%sfcwater_mass  (  1,donp)                 &
                                         + csite%sfcwater_mass  (iii,donp)
         csite%sfcwater_energy(iii,donp) = 0.
         csite%sfcwater_depth (iii,donp) = 0.
         csite%sfcwater_mass  (iii,donp) = 0.
      end do
      do iii=2,csite%nlev_sfcwater(recp)
         csite%sfcwater_energy(  1,recp) = csite%sfcwater_energy(  1,recp)                 &
                                         + csite%sfcwater_energy(iii,recp)
         csite%sfcwater_depth (  1,recp) = csite%sfcwater_depth (  1,recp)                 &
                                         + csite%sfcwater_depth (iii,recp)
         csite%sfcwater_mass  (  1,recp) = csite%sfcwater_mass  (  1,recp)                 &
                                         + csite%sfcwater_mass  (iii,recp)
         csite%sfcwater_energy(iii,recp) = 0.
         csite%sfcwater_depth (iii,recp) = 0.
         csite%sfcwater_mass  (iii,recp) = 0.
      end do
      !----- 3. Merge the patches; --------------------------------------------------------!
      if (csite%nlev_sfcwater(recp) > 0 .or. csite%nlev_sfcwater(donp) > 0) then
         csite%sfcwater_mass  (1,recp) = csite%sfcwater_mass  (1,donp) * dawgt             &
                                       + csite%sfcwater_mass  (1,recp) * rawgt
         csite%sfcwater_depth (1,recp) = csite%sfcwater_depth (1,donp) * dawgt             &
                                       + csite%sfcwater_depth (1,recp) * rawgt
         csite%sfcwater_energy(1,recp) = csite%sfcwater_energy(1,donp) * dawgt             &
                                       + csite%sfcwater_energy(1,recp) * rawgt
      else
         csite%sfcwater_mass  (1,recp) = 0.
         csite%sfcwater_depth (1,recp) = 0.
         csite%sfcwater_energy(1,recp) = 0.
      end if
      !------------------------------------------------------------------------------------!
      ! 4. Convert energy back to J/kg;                                                    !
      ! 5. Find temperature and liquid water fraction;                                     !
      !    (Both are done in new_patch_sfc_props).                                         !
      !------------------------------------------------------------------------------------!
      !------------------------------------------------------------------------------------!
     

      !------------------------------------------------------------------------------------!
      !     Merge soil energy and water, then find temperature, liquid fraction and soil   !
      ! matric potential.                                                                  !
      !------------------------------------------------------------------------------------!
      do iii=1,mzg
         csite%soil_energy(iii,recp) = csite%soil_energy(iii,donp) * dawgt                 &
                                     + csite%soil_energy(iii,recp) * rawgt
         csite%soil_water (iii,recp) = csite%soil_water (iii,donp) * dawgt                 &
                                     + csite%soil_water (iii,recp) * rawgt
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     These variables shouldn't matter because they are reset every day/every month, !
      ! but just in case...                                                                !
      !------------------------------------------------------------------------------------!
      csite%avg_daily_temp       (recp) = csite%avg_daily_temp       (donp) * dawgt        &
                                        + csite%avg_daily_temp       (recp) * rawgt

      csite%avg_monthly_gndwater (recp) = csite%avg_monthly_gndwater (donp) * dawgt        &
                                        + csite%avg_monthly_gndwater (recp) * rawgt

      csite%avg_monthly_waterdef (recp) = csite%avg_monthly_waterdef (donp) * dawgt        &
                                        + csite%avg_monthly_waterdef (recp) * rawgt
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
      ! + csite%soil_mstpot(k,recp)                                                        !
      ! + csite%nlev_sfcwater(recp)                                                        !
      ! + csite%sfcwater_energy(k,recp) (Just converting back to J/kg)                     !
      ! + csite%csite%sfcwater_tempk(k,recp)                                               !
      ! + csite%sfcwater_fracliq(k,recp)                                                   !
      !------------------------------------------------------------------------------------!
      call new_patch_sfc_props(csite,recp,mzg,mzs,ntext_soil)
      !------------------------------------------------------------------------------------!





      !------ Seed bank. ------------------------------------------------------------------!
      do iii = 1,n_pft
         csite%repro(iii,recp)           = csite%repro(iii,donp)      * dawgt              &
                                         + csite%repro(iii,recp)      * rawgt
      end do
      !------------------------------------------------------------------------------------!






      !------------------------------------------------------------------------------------!
      !    Sub-daily means.  This should be skipped during initialisation because all      !
      ! averages will be zero.                                                             !
      !------------------------------------------------------------------------------------!
      if (.not. fuse_initial) then
         csite%fmean_rh               (recp) = csite%fmean_rh               (donp) * dawgt &
                                             + csite%fmean_rh               (recp) * rawgt
         csite%fmean_fgc_rh           (recp) = csite%fmean_fgc_rh           (donp) * dawgt &
                                             + csite%fmean_fgc_rh           (recp) * rawgt
         csite%fmean_fsc_rh           (recp) = csite%fmean_fsc_rh           (donp) * dawgt &
                                             + csite%fmean_fsc_rh           (recp) * rawgt
         csite%fmean_stgc_rh          (recp) = csite%fmean_stgc_rh          (donp) * dawgt &
                                             + csite%fmean_stgc_rh          (recp) * rawgt
         csite%fmean_stsc_rh          (recp) = csite%fmean_stsc_rh          (donp) * dawgt &
                                             + csite%fmean_stsc_rh          (recp) * rawgt
         csite%fmean_msc_rh           (recp) = csite%fmean_msc_rh           (donp) * dawgt &
                                             + csite%fmean_msc_rh           (recp) * rawgt
         csite%fmean_ssc_rh           (recp) = csite%fmean_ssc_rh           (donp) * dawgt &
                                             + csite%fmean_ssc_rh           (recp) * rawgt
         csite%fmean_psc_rh           (recp) = csite%fmean_psc_rh           (donp) * dawgt &
                                             + csite%fmean_psc_rh           (recp) * rawgt
         csite%fmean_nep              (recp) = csite%fmean_nep              (donp) * dawgt &
                                             + csite%fmean_nep              (recp) * rawgt
         csite%fmean_rk4step          (recp) = csite%fmean_rk4step          (donp) * dawgt &
                                             + csite%fmean_rk4step          (recp) * rawgt
         csite%fmean_available_water  (recp) = csite%fmean_available_water  (donp) * dawgt &
                                             + csite%fmean_available_water  (recp) * rawgt
         csite%fmean_veg_displace     (recp) = csite%fmean_veg_displace     (donp) * dawgt &
                                             + csite%fmean_veg_displace     (recp) * rawgt
         csite%fmean_rough            (recp) = csite%fmean_rough            (donp) * dawgt &
                                             + csite%fmean_rough            (recp) * rawgt
         csite%fmean_can_theiv        (recp) = csite%fmean_can_theiv        (donp) * dawgt &
                                             + csite%fmean_can_theiv        (recp) * rawgt
         csite%fmean_can_theta        (recp) = csite%fmean_can_theta        (donp) * dawgt &
                                             + csite%fmean_can_theta        (recp) * rawgt
         csite%fmean_can_vpdef        (recp) = csite%fmean_can_vpdef        (donp) * dawgt &
                                             + csite%fmean_can_vpdef        (recp) * rawgt
         csite%fmean_can_shv          (recp) = csite%fmean_can_shv          (donp) * dawgt &
                                             + csite%fmean_can_shv          (recp) * rawgt
         csite%fmean_can_co2          (recp) = csite%fmean_can_co2          (donp) * dawgt &
                                             + csite%fmean_can_co2          (recp) * rawgt
         csite%fmean_can_prss         (recp) = csite%fmean_can_prss         (donp) * dawgt &
                                             + csite%fmean_can_prss         (recp) * rawgt
         csite%fmean_gnd_temp         (recp) = csite%fmean_gnd_temp         (donp) * dawgt &
                                             + csite%fmean_gnd_temp         (recp) * rawgt
         csite%fmean_gnd_shv          (recp) = csite%fmean_gnd_shv          (donp) * dawgt &
                                             + csite%fmean_gnd_shv          (recp) * rawgt
         csite%fmean_can_ggnd         (recp) = csite%fmean_can_ggnd         (donp) * dawgt &
                                             + csite%fmean_can_ggnd         (recp) * rawgt
         csite%fmean_rshort_gnd       (recp) = csite%fmean_rshort_gnd       (donp) * dawgt &
                                             + csite%fmean_rshort_gnd       (recp) * rawgt
         csite%fmean_par_gnd          (recp) = csite%fmean_par_gnd          (donp) * dawgt &
                                             + csite%fmean_par_gnd          (recp) * rawgt
         csite%fmean_rlong_gnd        (recp) = csite%fmean_rlong_gnd        (donp) * dawgt &
                                             + csite%fmean_rlong_gnd        (recp) * rawgt
         csite%fmean_rlongup          (recp) = csite%fmean_rlongup          (donp) * dawgt &
                                             + csite%fmean_rlongup          (recp) * rawgt
         csite%fmean_parup            (recp) = csite%fmean_parup            (donp) * dawgt &
                                             + csite%fmean_parup            (recp) * rawgt
         csite%fmean_nirup            (recp) = csite%fmean_nirup            (donp) * dawgt &
                                             + csite%fmean_nirup            (recp) * rawgt
         csite%fmean_rshortup         (recp) = csite%fmean_rshortup         (donp) * dawgt &
                                             + csite%fmean_rshortup         (recp) * rawgt
         csite%fmean_rnet             (recp) = csite%fmean_rnet             (donp) * dawgt &
                                             + csite%fmean_rnet             (recp) * rawgt
         csite%fmean_albedo           (recp) = csite%fmean_albedo           (donp) * dawgt &
                                             + csite%fmean_albedo           (recp) * rawgt
         csite%fmean_albedo_par       (recp) = csite%fmean_albedo_par       (donp) * dawgt &
                                             + csite%fmean_albedo_par       (recp) * rawgt
         csite%fmean_albedo_nir       (recp) = csite%fmean_albedo_nir       (donp) * dawgt &
                                             + csite%fmean_albedo_nir       (recp) * rawgt
         csite%fmean_rlong_albedo     (recp) = csite%fmean_rlong_albedo     (donp) * dawgt &
                                             + csite%fmean_rlong_albedo     (recp) * rawgt
         csite%fmean_ustar            (recp) = csite%fmean_ustar            (donp) * dawgt &
                                             + csite%fmean_ustar            (recp) * rawgt
         csite%fmean_tstar            (recp) = csite%fmean_tstar            (donp) * dawgt &
                                             + csite%fmean_tstar            (recp) * rawgt
         csite%fmean_qstar            (recp) = csite%fmean_qstar            (donp) * dawgt &
                                             + csite%fmean_qstar            (recp) * rawgt
         csite%fmean_cstar            (recp) = csite%fmean_cstar            (donp) * dawgt &
                                             + csite%fmean_cstar            (recp) * rawgt
         csite%fmean_carbon_ac        (recp) = csite%fmean_carbon_ac        (donp) * dawgt &
                                             + csite%fmean_carbon_ac        (recp) * rawgt
         csite%fmean_carbon_st        (recp) = csite%fmean_carbon_st        (donp) * dawgt &
                                             + csite%fmean_carbon_st        (recp) * rawgt
         csite%fmean_vapor_gc         (recp) = csite%fmean_vapor_gc         (donp) * dawgt &
                                             + csite%fmean_vapor_gc         (recp) * rawgt
         csite%fmean_vapor_ac         (recp) = csite%fmean_vapor_ac         (donp) * dawgt &
                                             + csite%fmean_vapor_ac         (recp) * rawgt
         csite%fmean_throughfall      (recp) = csite%fmean_throughfall      (donp) * dawgt &
                                             + csite%fmean_throughfall      (recp) * rawgt
         csite%fmean_runoff           (recp) = csite%fmean_runoff           (donp) * dawgt &
                                             + csite%fmean_runoff           (recp) * rawgt
         csite%fmean_drainage         (recp) = csite%fmean_drainage         (donp) * dawgt &
                                             + csite%fmean_drainage         (recp) * rawgt
         csite%fmean_sensible_gc      (recp) = csite%fmean_sensible_gc      (donp) * dawgt &
                                             + csite%fmean_sensible_gc      (recp) * rawgt
         csite%fmean_sensible_ac      (recp) = csite%fmean_sensible_ac      (donp) * dawgt &
                                             + csite%fmean_sensible_ac      (recp) * rawgt
         csite%fmean_qthroughfall     (recp) = csite%fmean_qthroughfall     (donp) * dawgt &
                                             + csite%fmean_qthroughfall     (recp) * rawgt
         csite%fmean_qrunoff          (recp) = csite%fmean_qrunoff          (donp) * dawgt &
                                             + csite%fmean_qrunoff          (recp) * rawgt
         csite%fmean_qdrainage        (recp) = csite%fmean_qdrainage        (donp) * dawgt &
                                             + csite%fmean_qdrainage        (recp) * rawgt
         csite%fmean_soil_energy    (:,recp) = csite%fmean_soil_energy    (:,donp) * dawgt &
                                             + csite%fmean_soil_energy    (:,recp) * rawgt
         csite%fmean_soil_water     (:,recp) = csite%fmean_soil_water     (:,donp) * dawgt &
                                             + csite%fmean_soil_water     (:,recp) * rawgt
         csite%fmean_smoist_gg      (:,recp) = csite%fmean_smoist_gg      (:,donp) * dawgt &
                                             + csite%fmean_smoist_gg      (:,recp) * rawgt
         csite%fmean_transloss      (:,recp) = csite%fmean_transloss      (:,donp) * dawgt &
                                             + csite%fmean_transloss      (:,recp) * rawgt
         csite%fmean_sensible_gg    (:,recp) = csite%fmean_sensible_gg    (:,donp) * dawgt &
                                             + csite%fmean_sensible_gg    (:,recp) * rawgt
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Now we find the derived properties for the canopy air space.               !
         !---------------------------------------------------------------------------------!

         xmean_can_exner             = press2exner (csite%fmean_can_prss(recp))
         csite%fmean_can_temp (recp) = extheta2temp( xmean_can_exner                       &
                                                   , csite%fmean_can_theta (recp)          )
         csite%fmean_can_rhos (recp) = idealdenssh ( csite%fmean_can_prss  (recp)          &
                                                   , csite%fmean_can_temp  (recp)          &
                                                   , csite%fmean_can_shv   (recp)          )
         csite%fmean_can_dmol (recp) = idealdmolsh ( csite%fmean_can_prss  (recp)          &
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
         csite%fmean_sfcw_depth (recp) = csite%fmean_sfcw_depth        (donp) * dawgt      &
                                       + csite%fmean_sfcw_depth        (recp) * rawgt
         csite%fmean_sfcw_energy(recp) = csite%fmean_sfcw_energy       (donp)              &
                                       * csite%fmean_sfcw_mass         (donp) * dawgt      &
                                       + csite%fmean_sfcw_energy       (recp)              &
                                       * csite%fmean_sfcw_mass         (recp) * rawgt
         csite%fmean_sfcw_mass  (recp) = csite%fmean_sfcw_mass         (donp) * dawgt      &
                                       + csite%fmean_sfcw_mass         (recp) * rawgt
         csite%fmean_snowfac    (recp) = csite%fmean_snowfac           (donp) * dawgt      &
                                       + csite%fmean_snowfac           (recp) * rawgt
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
         if ( all(csite%dmean_can_prss > 10.0) ) then
            csite%dmean_A_decomp       (recp) = csite%dmean_A_decomp       (donp) * dawgt  &
                                              + csite%dmean_A_decomp       (recp) * rawgt
            csite%dmean_B_decomp       (recp) = csite%dmean_B_decomp       (donp) * dawgt  &
                                              + csite%dmean_B_decomp       (recp) * rawgt
            csite%dmean_Af_decomp      (recp) = csite%dmean_Af_decomp      (donp) * dawgt  &
                                              + csite%dmean_Af_decomp      (recp) * rawgt
            csite%dmean_Bf_decomp      (recp) = csite%dmean_Bf_decomp      (donp) * dawgt  &
                                              + csite%dmean_Bf_decomp      (recp) * rawgt
            csite%dmean_co2_residual   (recp) = csite%dmean_co2_residual   (donp) * dawgt  &
                                              + csite%dmean_co2_residual   (recp) * rawgt
            csite%dmean_energy_residual(recp) = csite%dmean_energy_residual(donp) * dawgt  &
                                              + csite%dmean_energy_residual(recp) * rawgt
            csite%dmean_water_residual (recp) = csite%dmean_water_residual (donp) * dawgt  &
                                              + csite%dmean_water_residual (recp) * rawgt
            csite%dmean_rh             (recp) = csite%dmean_rh             (donp) * dawgt  &
                                              + csite%dmean_rh             (recp) * rawgt
            csite%dmean_fgc_rh         (recp) = csite%dmean_fgc_rh         (donp) * dawgt  &
                                              + csite%dmean_fgc_rh         (recp) * rawgt
            csite%dmean_fsc_rh         (recp) = csite%dmean_fsc_rh         (donp) * dawgt  &
                                              + csite%dmean_fsc_rh         (recp) * rawgt
            csite%dmean_stgc_rh        (recp) = csite%dmean_stgc_rh        (donp) * dawgt  &
                                              + csite%dmean_stgc_rh        (recp) * rawgt
            csite%dmean_stsc_rh        (recp) = csite%dmean_stsc_rh        (donp) * dawgt  &
                                              + csite%dmean_stsc_rh        (recp) * rawgt
            csite%dmean_msc_rh         (recp) = csite%dmean_msc_rh         (donp) * dawgt  &
                                              + csite%dmean_msc_rh         (recp) * rawgt
            csite%dmean_ssc_rh         (recp) = csite%dmean_ssc_rh         (donp) * dawgt  &
                                              + csite%dmean_ssc_rh         (recp) * rawgt
            csite%dmean_psc_rh         (recp) = csite%dmean_psc_rh         (donp) * dawgt  &
                                              + csite%dmean_psc_rh         (recp) * rawgt
            csite%dmean_nep            (recp) = csite%dmean_nep            (donp) * dawgt  &
                                              + csite%dmean_nep            (recp) * rawgt
            csite%dmean_rk4step        (recp) = csite%dmean_rk4step        (donp) * dawgt  &
                                              + csite%dmean_rk4step        (recp) * rawgt
            csite%dmean_available_water(recp) = csite%dmean_available_water(donp) * dawgt  &
                                              + csite%dmean_available_water(recp) * rawgt
            csite%dmean_veg_displace   (recp) = csite%dmean_veg_displace   (donp) * dawgt  &
                                              + csite%dmean_veg_displace   (recp) * rawgt
            csite%dmean_rough          (recp) = csite%dmean_rough          (donp) * dawgt  &
                                              + csite%dmean_rough          (recp) * rawgt
            csite%dmean_can_theiv      (recp) = csite%dmean_can_theiv      (donp) * dawgt  &
                                              + csite%dmean_can_theiv      (recp) * rawgt
            csite%dmean_can_theta      (recp) = csite%dmean_can_theta      (donp) * dawgt  &
                                              + csite%dmean_can_theta      (recp) * rawgt
            csite%dmean_can_vpdef      (recp) = csite%dmean_can_vpdef      (donp) * dawgt  &
                                              + csite%dmean_can_vpdef      (recp) * rawgt
            csite%dmean_can_shv        (recp) = csite%dmean_can_shv        (donp) * dawgt  &
                                              + csite%dmean_can_shv        (recp) * rawgt
            csite%dmean_can_co2        (recp) = csite%dmean_can_co2        (donp) * dawgt  &
                                              + csite%dmean_can_co2        (recp) * rawgt
            csite%dmean_can_prss       (recp) = csite%dmean_can_prss       (donp) * dawgt  &
                                              + csite%dmean_can_prss       (recp) * rawgt
            csite%dmean_gnd_temp       (recp) = csite%dmean_gnd_temp       (donp) * dawgt  &
                                              + csite%dmean_gnd_temp       (recp) * rawgt
            csite%dmean_gnd_shv        (recp) = csite%dmean_gnd_shv        (donp) * dawgt  &
                                              + csite%dmean_gnd_shv        (recp) * rawgt
            csite%dmean_can_ggnd       (recp) = csite%dmean_can_ggnd       (donp) * dawgt  &
                                              + csite%dmean_can_ggnd       (recp) * rawgt
            csite%dmean_rshort_gnd     (recp) = csite%dmean_rshort_gnd     (donp) * dawgt  &
                                              + csite%dmean_rshort_gnd     (recp) * rawgt
            csite%dmean_par_gnd        (recp) = csite%dmean_par_gnd        (donp) * dawgt  &
                                              + csite%dmean_par_gnd        (recp) * rawgt
            csite%dmean_rlong_gnd      (recp) = csite%dmean_rlong_gnd      (donp) * dawgt  &
                                              + csite%dmean_rlong_gnd      (recp) * rawgt
            csite%dmean_rlongup        (recp) = csite%dmean_rlongup        (donp) * dawgt  &
                                              + csite%dmean_rlongup        (recp) * rawgt
            csite%dmean_parup          (recp) = csite%dmean_parup          (donp) * dawgt  &
                                              + csite%dmean_parup          (recp) * rawgt
            csite%dmean_nirup          (recp) = csite%dmean_nirup          (donp) * dawgt  &
                                              + csite%dmean_nirup          (recp) * rawgt
            csite%dmean_rshortup       (recp) = csite%dmean_rshortup       (donp) * dawgt  &
                                              + csite%dmean_rshortup       (recp) * rawgt
            csite%dmean_rnet           (recp) = csite%dmean_rnet           (donp) * dawgt  &
                                              + csite%dmean_rnet           (recp) * rawgt
            csite%dmean_albedo         (recp) = csite%dmean_albedo         (donp) * dawgt  &
                                              + csite%dmean_albedo         (recp) * rawgt
            csite%dmean_albedo_par     (recp) = csite%dmean_albedo_par     (donp) * dawgt  &
                                              + csite%dmean_albedo_par     (recp) * rawgt
            csite%dmean_albedo_nir     (recp) = csite%dmean_albedo_nir     (donp) * dawgt  &
                                              + csite%dmean_albedo_nir     (recp) * rawgt
            csite%dmean_rlong_albedo   (recp) = csite%dmean_rlong_albedo   (donp) * dawgt  &
                                              + csite%dmean_rlong_albedo   (recp) * rawgt
            csite%dmean_ustar          (recp) = csite%dmean_ustar          (donp) * dawgt  &
                                              + csite%dmean_ustar          (recp) * rawgt
            csite%dmean_tstar          (recp) = csite%dmean_tstar          (donp) * dawgt  &
                                              + csite%dmean_tstar          (recp) * rawgt
            csite%dmean_qstar          (recp) = csite%dmean_qstar          (donp) * dawgt  &
                                              + csite%dmean_qstar          (recp) * rawgt
            csite%dmean_cstar          (recp) = csite%dmean_cstar          (donp) * dawgt  &
                                              + csite%dmean_cstar          (recp) * rawgt
            csite%dmean_carbon_ac      (recp) = csite%dmean_carbon_ac      (donp) * dawgt  &
                                              + csite%dmean_carbon_ac      (recp) * rawgt
            csite%dmean_carbon_st      (recp) = csite%dmean_carbon_st      (donp) * dawgt  &
                                              + csite%dmean_carbon_st      (recp) * rawgt
            csite%dmean_vapor_gc       (recp) = csite%dmean_vapor_gc       (donp) * dawgt  &
                                              + csite%dmean_vapor_gc       (recp) * rawgt
            csite%dmean_vapor_ac       (recp) = csite%dmean_vapor_ac       (donp) * dawgt  &
                                              + csite%dmean_vapor_ac       (recp) * rawgt
            csite%dmean_throughfall    (recp) = csite%dmean_throughfall    (donp) * dawgt  &
                                              + csite%dmean_throughfall    (recp) * rawgt
            csite%dmean_runoff         (recp) = csite%dmean_runoff         (donp) * dawgt  &
                                              + csite%dmean_runoff         (recp) * rawgt
            csite%dmean_drainage       (recp) = csite%dmean_drainage       (donp) * dawgt  &
                                              + csite%dmean_drainage       (recp) * rawgt
            csite%dmean_sensible_gc    (recp) = csite%dmean_sensible_gc    (donp) * dawgt  &
                                              + csite%dmean_sensible_gc    (recp) * rawgt
            csite%dmean_sensible_ac    (recp) = csite%dmean_sensible_ac    (donp) * dawgt  &
                                              + csite%dmean_sensible_ac    (recp) * rawgt
            csite%dmean_qthroughfall   (recp) = csite%dmean_qthroughfall   (donp) * dawgt  &
                                              + csite%dmean_qthroughfall   (recp) * rawgt
            csite%dmean_qrunoff        (recp) = csite%dmean_qrunoff        (donp) * dawgt  &
                                              + csite%dmean_qrunoff        (recp) * rawgt
            csite%dmean_qdrainage      (recp) = csite%dmean_qdrainage      (donp) * dawgt  &
                                              + csite%dmean_qdrainage      (recp) * rawgt
            csite%dmean_smoist_gg    (:,recp) = csite%dmean_smoist_gg    (:,donp) * dawgt  &
                                              + csite%dmean_smoist_gg    (:,recp) * rawgt
            csite%dmean_transloss    (:,recp) = csite%dmean_transloss    (:,donp) * dawgt  &
                                              + csite%dmean_transloss    (:,recp) * rawgt
            csite%dmean_sensible_gg  (:,recp) = csite%dmean_sensible_gg  (:,donp) * dawgt  &
                                              + csite%dmean_sensible_gg  (:,recp) * rawgt
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Now we find the derived properties for the canopy air space.            !
            !------------------------------------------------------------------------------!
            xmean_can_exner             = press2exner (csite%dmean_can_prss(recp))
            csite%dmean_can_temp (recp) = extheta2temp( xmean_can_exner                    &
                                                      , csite%dmean_can_theta(recp)        )
            csite%dmean_can_rhos (recp) = idealdenssh ( csite%dmean_can_prss  (recp)       &
                                                      , csite%dmean_can_temp  (recp)       &
                                                      , csite%dmean_can_shv   (recp)       )
            csite%dmean_can_dmol (recp) = idealdmolsh ( csite%dmean_can_prss  (recp)       &
                                                      , csite%dmean_can_temp  (recp)       &
                                                      , csite%dmean_can_shv   (recp)       )
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Find the soil mean temperature, liquid water fraction, and matric       !
            ! potential.                                                                   !
            !------------------------------------------------------------------------------!
            do iii=lsl,mzg
               nsoil = ntext_soil(iii)
               call uextcm2tl( csite%dmean_soil_energy(iii,recp)                           &
                             , csite%dmean_soil_water (iii,recp) * wdns                    &
                             , soil(nsoil)%slcpd                                           &
                             , csite%dmean_soil_temp  (iii,recp)                           &
                             , csite%dmean_soil_fliq  (iii,recp))

               csite%dmean_soil_mstpot(iii,recp)  =                                        &
                                   matric_potential(nsoil,csite%dmean_soil_water(iii,recp))
            end do
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find the temporary surface water properties.  They may not be available  !
            ! at all times, so we must check.                                              !
            !------------------------------------------------------------------------------!
            !----- Temporarily make energy extensive [J/m2]. ------------------------------!
            csite%dmean_sfcw_depth (recp) = csite%dmean_sfcw_depth (donp) * dawgt          &
                                          + csite%dmean_sfcw_depth (recp) * rawgt
            csite%dmean_sfcw_energy(recp) = csite%dmean_sfcw_energy(donp)                  &
                                          * csite%dmean_sfcw_mass  (donp) * dawgt          &
                                          + csite%dmean_sfcw_energy(recp)                  &
                                          * csite%dmean_sfcw_mass  (recp) * rawgt
            csite%dmean_sfcw_mass  (recp) = csite%dmean_sfcw_mass  (donp) * dawgt          &
                                          + csite%dmean_sfcw_mass  (recp) * rawgt
            csite%dmean_snowfac    (recp) = csite%dmean_snowfac    (donp) * dawgt          &
                                          + csite%dmean_snowfac    (recp) * rawgt
            !----- Check whether there is enough surface water. ---------------------------!
            if (csite%dmean_sfcw_mass(recp) > tiny_sfcwater_mass) then
               csite%dmean_sfcw_energy     (recp) =   csite%dmean_sfcw_energy     (recp)   &
                                                  /   csite%dmean_sfcw_mass       (recp)
               call uint2tl( csite%dmean_sfcw_energy(recp)                                 &
                           , csite%dmean_sfcw_temp  (recp)                                 &
                           , csite%dmean_sfcw_fliq  (recp) )
            else
               csite%dmean_sfcw_mass  (recp)  = 0.
               csite%dmean_sfcw_depth (recp)  = 0.
               csite%dmean_sfcw_energy(recp)  = 0.
               csite%dmean_sfcw_temp  (recp)  = csite%dmean_soil_temp(mzg,recp)
               csite%dmean_sfcw_fliq  (recp)  = csite%dmean_soil_fliq(mzg,recp)
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !    Monthly means.                                                                  !
      !------------------------------------------------------------------------------------!
      if (writing_eorq .and. (.not. fuse_initial)) then
        if ( all(csite%mmean_can_prss > 10.0) ) then
            !------------------------------------------------------------------------------!
            !    First we find the mean sum of squares, because they depend on the means   !
            ! too, and the receptor values are lost after fusion.  All variables are       !
            ! intensive at the patch level.                                                !
            !------------------------------------------------------------------------------!
            csite%mmsqu_rh         (recp) = fuse_msqu( csite%mmean_rh         (donp)       &
                                                     , csite%mmsqu_rh         (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_rh         (recp)       &
                                                     , csite%mmsqu_rh         (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_fgc_rh     (recp) = fuse_msqu( csite%mmean_fgc_rh     (donp)       &
                                                     , csite%mmsqu_fgc_rh     (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_fgc_rh     (recp)       &
                                                     , csite%mmsqu_fgc_rh     (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_fsc_rh     (recp) = fuse_msqu( csite%mmean_fsc_rh     (donp)       &
                                                     , csite%mmsqu_fsc_rh     (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_fsc_rh     (recp)       &
                                                     , csite%mmsqu_fsc_rh     (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_stgc_rh    (recp) = fuse_msqu( csite%mmean_stgc_rh    (donp)       &
                                                     , csite%mmsqu_stgc_rh    (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_stgc_rh    (recp)       &
                                                     , csite%mmsqu_stgc_rh    (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_stsc_rh    (recp) = fuse_msqu( csite%mmean_stsc_rh    (donp)       &
                                                     , csite%mmsqu_stsc_rh    (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_stsc_rh    (recp)       &
                                                     , csite%mmsqu_stsc_rh    (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_msc_rh     (recp) = fuse_msqu( csite%mmean_msc_rh     (donp)       &
                                                     , csite%mmsqu_msc_rh     (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_msc_rh     (recp)       &
                                                     , csite%mmsqu_msc_rh     (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_ssc_rh     (recp) = fuse_msqu( csite%mmean_ssc_rh     (donp)       &
                                                     , csite%mmsqu_ssc_rh     (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_ssc_rh     (recp)       &
                                                     , csite%mmsqu_ssc_rh     (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_psc_rh     (recp) = fuse_msqu( csite%mmean_psc_rh     (donp)       &
                                                     , csite%mmsqu_psc_rh     (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_psc_rh     (recp)       &
                                                     , csite%mmsqu_psc_rh     (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_nep        (recp) = fuse_msqu( csite%mmean_nep        (donp)       &
                                                     , csite%mmsqu_nep        (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_nep        (recp)       &
                                                     , csite%mmsqu_nep        (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_rlongup    (recp) = fuse_msqu( csite%mmean_rlongup    (donp)       &
                                                     , csite%mmsqu_rlongup    (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_rlongup    (recp)       &
                                                     , csite%mmsqu_rlongup    (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_parup      (recp) = fuse_msqu( csite%mmean_parup      (donp)       &
                                                     , csite%mmsqu_parup      (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_parup      (recp)       &
                                                     , csite%mmsqu_parup      (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_nirup      (recp) = fuse_msqu( csite%mmean_nirup      (donp)       &
                                                     , csite%mmsqu_nirup      (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_nirup      (recp)       &
                                                     , csite%mmsqu_nirup      (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_rshortup   (recp) = fuse_msqu( csite%mmean_rshortup   (donp)       &
                                                     , csite%mmsqu_rshortup   (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_rshortup   (recp)       &
                                                     , csite%mmsqu_rshortup   (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_rnet       (recp) = fuse_msqu( csite%mmean_rnet       (donp)       &
                                                     , csite%mmsqu_rnet       (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_rnet       (recp)       &
                                                     , csite%mmsqu_rnet       (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_albedo     (recp) = fuse_msqu( csite%mmean_albedo     (donp)       &
                                                     , csite%mmsqu_albedo     (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_albedo     (recp)       &
                                                     , csite%mmsqu_albedo     (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_ustar      (recp) = fuse_msqu( csite%mmean_ustar      (donp)       &
                                                     , csite%mmsqu_ustar      (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_ustar      (recp)       &
                                                     , csite%mmsqu_ustar      (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_carbon_ac  (recp) = fuse_msqu( csite%mmean_carbon_ac  (donp)       &
                                                     , csite%mmsqu_carbon_ac  (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_carbon_ac  (recp)       &
                                                     , csite%mmsqu_carbon_ac  (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_carbon_st  (recp) = fuse_msqu( csite%mmean_carbon_st  (donp)       &
                                                     , csite%mmsqu_carbon_st  (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_carbon_st  (recp)       &
                                                     , csite%mmsqu_carbon_st  (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_vapor_gc   (recp) = fuse_msqu( csite%mmean_vapor_gc   (donp)       &
                                                     , csite%mmsqu_vapor_gc   (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_vapor_gc   (recp)       &
                                                     , csite%mmsqu_vapor_gc   (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_vapor_ac   (recp) = fuse_msqu( csite%mmean_vapor_ac   (donp)       &
                                                     , csite%mmsqu_vapor_ac   (donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_vapor_ac   (recp)       &
                                                     , csite%mmsqu_vapor_ac   (recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_sensible_gc(recp) = fuse_msqu( csite%mmean_sensible_gc(donp)       &
                                                     , csite%mmsqu_sensible_gc(donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_sensible_gc(recp)       &
                                                     , csite%mmsqu_sensible_gc(recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            csite%mmsqu_sensible_ac(recp) = fuse_msqu( csite%mmean_sensible_ac(donp)       &
                                                     , csite%mmsqu_sensible_ac(donp)       &
                                                     , dawgt                               &
                                                     , csite%mmean_sensible_ac(recp)       &
                                                     , csite%mmsqu_sensible_ac(recp)       &
                                                     , rawgt                               &
                                                     , corr_patch, .false.)
            !------------------------------------------------------------------------------!


            csite%mmean_rh             (recp) = csite%mmean_rh             (donp) * dawgt  &
                                              + csite%mmean_rh             (recp) * rawgt
            csite%mmean_fgc_rh         (recp) = csite%mmean_fgc_rh         (donp) * dawgt  &
                                              + csite%mmean_fgc_rh         (recp) * rawgt
            csite%mmean_fsc_rh         (recp) = csite%mmean_fsc_rh         (donp) * dawgt  &
                                              + csite%mmean_fsc_rh         (recp) * rawgt
            csite%mmean_stgc_rh        (recp) = csite%mmean_stgc_rh        (donp) * dawgt  &
                                              + csite%mmean_stgc_rh        (recp) * rawgt
            csite%mmean_stsc_rh        (recp) = csite%mmean_stsc_rh        (donp) * dawgt  &
                                              + csite%mmean_stsc_rh        (recp) * rawgt
            csite%mmean_msc_rh         (recp) = csite%mmean_msc_rh         (donp) * dawgt  &
                                              + csite%mmean_msc_rh         (recp) * rawgt
            csite%mmean_ssc_rh         (recp) = csite%mmean_ssc_rh         (donp) * dawgt  &
                                              + csite%mmean_ssc_rh         (recp) * rawgt
            csite%mmean_psc_rh         (recp) = csite%mmean_psc_rh         (donp) * dawgt  &
                                              + csite%mmean_psc_rh         (recp) * rawgt
            csite%mmean_nep            (recp) = csite%mmean_nep            (donp) * dawgt  &
                                              + csite%mmean_nep            (recp) * rawgt
            csite%mmean_rk4step        (recp) = csite%mmean_rk4step        (donp) * dawgt  &
                                              + csite%mmean_rk4step        (recp) * rawgt
            csite%mmean_available_water(recp) = csite%mmean_available_water(donp) * dawgt  &
                                              + csite%mmean_available_water(recp) * rawgt
            csite%mmean_veg_displace   (recp) = csite%mmean_veg_displace   (donp) * dawgt  &
                                              + csite%mmean_veg_displace   (recp) * rawgt
            csite%mmean_rough          (recp) = csite%mmean_rough          (donp) * dawgt  &
                                              + csite%mmean_rough          (recp) * rawgt
            csite%mmean_can_theiv      (recp) = csite%mmean_can_theiv      (donp) * dawgt  &
                                              + csite%mmean_can_theiv      (recp) * rawgt
            csite%mmean_can_theta      (recp) = csite%mmean_can_theta      (donp) * dawgt  &
                                              + csite%mmean_can_theta      (recp) * rawgt
            csite%mmean_can_vpdef      (recp) = csite%mmean_can_vpdef      (donp) * dawgt  &
                                              + csite%mmean_can_vpdef      (recp) * rawgt
            csite%mmean_can_shv        (recp) = csite%mmean_can_shv        (donp) * dawgt  &
                                              + csite%mmean_can_shv        (recp) * rawgt
            csite%mmean_can_co2        (recp) = csite%mmean_can_co2        (donp) * dawgt  &
                                              + csite%mmean_can_co2        (recp) * rawgt
            csite%mmean_can_prss       (recp) = csite%mmean_can_prss       (donp) * dawgt  &
                                              + csite%mmean_can_prss       (recp) * rawgt
            csite%mmean_gnd_temp       (recp) = csite%mmean_gnd_temp       (donp) * dawgt  &
                                              + csite%mmean_gnd_temp       (recp) * rawgt
            csite%mmean_gnd_shv        (recp) = csite%mmean_gnd_shv        (donp) * dawgt  &
                                              + csite%mmean_gnd_shv        (recp) * rawgt
            csite%mmean_can_ggnd       (recp) = csite%mmean_can_ggnd       (donp) * dawgt  &
                                              + csite%mmean_can_ggnd       (recp) * rawgt
            csite%mmean_rshort_gnd     (recp) = csite%mmean_rshort_gnd     (donp) * dawgt  &
                                              + csite%mmean_rshort_gnd     (recp) * rawgt
            csite%mmean_par_gnd        (recp) = csite%mmean_par_gnd        (donp) * dawgt  &
                                              + csite%mmean_par_gnd        (recp) * rawgt
            csite%mmean_rlong_gnd      (recp) = csite%mmean_rlong_gnd      (donp) * dawgt  &
                                              + csite%mmean_rlong_gnd      (recp) * rawgt
            csite%mmean_rlongup        (recp) = csite%mmean_rlongup        (donp) * dawgt  &
                                              + csite%mmean_rlongup        (recp) * rawgt
            csite%mmean_parup          (recp) = csite%mmean_parup          (donp) * dawgt  &
                                              + csite%mmean_parup          (recp) * rawgt
            csite%mmean_nirup          (recp) = csite%mmean_nirup          (donp) * dawgt  &
                                              + csite%mmean_nirup          (recp) * rawgt
            csite%mmean_rshortup       (recp) = csite%mmean_rshortup       (donp) * dawgt  &
                                              + csite%mmean_rshortup       (recp) * rawgt
            csite%mmean_rnet           (recp) = csite%mmean_rnet           (donp) * dawgt  &
                                              + csite%mmean_rnet           (recp) * rawgt
            csite%mmean_albedo         (recp) = csite%mmean_albedo         (donp) * dawgt  &
                                              + csite%mmean_albedo         (recp) * rawgt
            csite%mmean_albedo_par     (recp) = csite%mmean_albedo_par     (donp) * dawgt  &
                                              + csite%mmean_albedo_par     (recp) * rawgt
            csite%mmean_albedo_nir     (recp) = csite%mmean_albedo_nir     (donp) * dawgt  &
                                              + csite%mmean_albedo_nir     (recp) * rawgt
            csite%mmean_rlong_albedo   (recp) = csite%mmean_rlong_albedo   (donp) * dawgt  &
                                              + csite%mmean_rlong_albedo   (recp) * rawgt
            csite%mmean_ustar          (recp) = csite%mmean_ustar          (donp) * dawgt  &
                                              + csite%mmean_ustar          (recp) * rawgt
            csite%mmean_tstar          (recp) = csite%mmean_tstar          (donp) * dawgt  &
                                              + csite%mmean_tstar          (recp) * rawgt
            csite%mmean_qstar          (recp) = csite%mmean_qstar          (donp) * dawgt  &
                                              + csite%mmean_qstar          (recp) * rawgt
            csite%mmean_cstar          (recp) = csite%mmean_cstar          (donp) * dawgt  &
                                              + csite%mmean_cstar          (recp) * rawgt
            csite%mmean_carbon_ac      (recp) = csite%mmean_carbon_ac      (donp) * dawgt  &
                                              + csite%mmean_carbon_ac      (recp) * rawgt
            csite%mmean_carbon_st      (recp) = csite%mmean_carbon_st      (donp) * dawgt  &
                                              + csite%mmean_carbon_st      (recp) * rawgt
            csite%mmean_vapor_gc       (recp) = csite%mmean_vapor_gc       (donp) * dawgt  &
                                              + csite%mmean_vapor_gc       (recp) * rawgt
            csite%mmean_vapor_ac       (recp) = csite%mmean_vapor_ac       (donp) * dawgt  &
                                              + csite%mmean_vapor_ac       (recp) * rawgt
            csite%mmean_throughfall    (recp) = csite%mmean_throughfall    (donp) * dawgt  &
                                              + csite%mmean_throughfall    (recp) * rawgt
            csite%mmean_runoff         (recp) = csite%mmean_runoff         (donp) * dawgt  &
                                              + csite%mmean_runoff         (recp) * rawgt
            csite%mmean_drainage       (recp) = csite%mmean_drainage       (donp) * dawgt  &
                                              + csite%mmean_drainage       (recp) * rawgt
            csite%mmean_sensible_gc    (recp) = csite%mmean_sensible_gc    (donp) * dawgt  &
                                              + csite%mmean_sensible_gc    (recp) * rawgt
            csite%mmean_sensible_ac    (recp) = csite%mmean_sensible_ac    (donp) * dawgt  &
                                              + csite%mmean_sensible_ac    (recp) * rawgt
            csite%mmean_qthroughfall   (recp) = csite%mmean_qthroughfall   (donp) * dawgt  &
                                              + csite%mmean_qthroughfall   (recp) * rawgt
            csite%mmean_qrunoff        (recp) = csite%mmean_qrunoff        (donp) * dawgt  &
                                              + csite%mmean_qrunoff        (recp) * rawgt
            csite%mmean_qdrainage      (recp) = csite%mmean_qdrainage      (donp) * dawgt  &
                                              + csite%mmean_qdrainage      (recp) * rawgt
            csite%mmean_fast_grnd_c    (recp) = csite%mmean_fast_grnd_c    (donp) * dawgt  &
                                              + csite%mmean_fast_grnd_c    (recp) * rawgt
            csite%mmean_fast_soil_c    (recp) = csite%mmean_fast_soil_c    (donp) * dawgt  &
                                              + csite%mmean_fast_soil_c    (recp) * rawgt
            csite%mmean_struct_grnd_c  (recp) = csite%mmean_struct_grnd_c  (donp) * dawgt  &
                                              + csite%mmean_struct_grnd_c  (recp) * rawgt
            csite%mmean_struct_soil_c  (recp) = csite%mmean_struct_soil_c  (donp) * dawgt  &
                                              + csite%mmean_struct_soil_c  (recp) * rawgt
            csite%mmean_struct_grnd_l  (recp) = csite%mmean_struct_grnd_l  (donp) * dawgt  &
                                              + csite%mmean_struct_grnd_l  (recp) * rawgt
            csite%mmean_struct_soil_l  (recp) = csite%mmean_struct_soil_l  (donp) * dawgt  &
                                              + csite%mmean_struct_soil_l  (recp) * rawgt
            csite%mmean_microbe_soil_c (recp) = csite%mmean_microbe_soil_c (donp) * dawgt  &
                                              + csite%mmean_microbe_soil_c (recp) * rawgt
            csite%mmean_slow_soil_c    (recp) = csite%mmean_slow_soil_c    (donp) * dawgt  &
                                              + csite%mmean_slow_soil_c    (recp) * rawgt
            csite%mmean_passive_soil_c (recp) = csite%mmean_passive_soil_c (donp) * dawgt  &
                                              + csite%mmean_passive_soil_c (recp) * rawgt
            csite%mmean_fast_grnd_n    (recp) = csite%mmean_fast_grnd_n    (donp) * dawgt  &
                                              + csite%mmean_fast_grnd_n    (recp) * rawgt
            csite%mmean_fast_soil_n    (recp) = csite%mmean_fast_soil_n    (donp) * dawgt  &
                                              + csite%mmean_fast_soil_n    (recp) * rawgt
            csite%mmean_struct_grnd_n  (recp) = csite%mmean_struct_grnd_n  (donp) * dawgt  &
                                              + csite%mmean_struct_grnd_n  (recp) * rawgt
            csite%mmean_struct_soil_n  (recp) = csite%mmean_struct_soil_n  (donp) * dawgt  &
                                              + csite%mmean_struct_soil_n  (recp) * rawgt
            csite%mmean_mineral_soil_n (recp) = csite%mmean_mineral_soil_n (donp) * dawgt  &
                                              + csite%mmean_mineral_soil_n (recp) * rawgt
            csite%mmean_fgc_in         (recp) = csite%mmean_fgc_in         (donp) * dawgt  &
                                              + csite%mmean_fgc_in         (recp) * rawgt
            csite%mmean_fsc_in         (recp) = csite%mmean_fsc_in         (donp) * dawgt  &
                                              + csite%mmean_fsc_in         (recp) * rawgt
            csite%mmean_stgc_in        (recp) = csite%mmean_stgc_in        (donp) * dawgt  &
                                              + csite%mmean_stgc_in        (recp) * rawgt
            csite%mmean_stsc_in        (recp) = csite%mmean_stsc_in        (donp) * dawgt  &
                                              + csite%mmean_stsc_in        (recp) * rawgt
            csite%mmean_A_decomp       (recp) = csite%mmean_A_decomp       (donp) * dawgt  &
                                              + csite%mmean_A_decomp       (recp) * rawgt
            csite%mmean_B_decomp       (recp) = csite%mmean_B_decomp       (donp) * dawgt  &
                                              + csite%mmean_B_decomp       (recp) * rawgt
            csite%mmean_Af_decomp      (recp) = csite%mmean_Af_decomp      (donp) * dawgt  &
                                              + csite%mmean_Af_decomp      (recp) * rawgt
            csite%mmean_Bf_decomp      (recp) = csite%mmean_Bf_decomp      (donp) * dawgt  &
                                              + csite%mmean_Bf_decomp      (recp) * rawgt
            csite%mmean_co2_residual   (recp) = csite%mmean_co2_residual   (donp) * dawgt  &
                                              + csite%mmean_co2_residual   (recp) * rawgt
            csite%mmean_energy_residual(recp) = csite%mmean_energy_residual(donp) * dawgt  &
                                              + csite%mmean_energy_residual(recp) * rawgt
            csite%mmean_water_residual (recp) = csite%mmean_water_residual (donp) * dawgt  &
                                              + csite%mmean_water_residual (recp) * rawgt
            csite%mmean_smoist_gg    (:,recp) = csite%mmean_smoist_gg    (:,donp) * dawgt  &
                                              + csite%mmean_smoist_gg    (:,recp) * rawgt
            csite%mmean_transloss    (:,recp) = csite%mmean_transloss    (:,donp) * dawgt  &
                                              + csite%mmean_transloss    (:,recp) * rawgt
            csite%mmean_sensible_gg  (:,recp) = csite%mmean_sensible_gg  (:,donp) * dawgt  &
                                              + csite%mmean_sensible_gg  (:,recp) * rawgt
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Now we find the derived properties for the canopy air space.            !
            !------------------------------------------------------------------------------!
            xmean_can_exner             = press2exner (csite%mmean_can_prss(recp))
            csite%mmean_can_temp (recp) = extheta2temp( xmean_can_exner                    &
                                                      , csite%mmean_can_theta (recp)       )
            csite%mmean_can_rhos (recp) = idealdenssh ( csite%mmean_can_prss  (recp)       &
                                                      , csite%mmean_can_temp  (recp)       &
                                                      , csite%mmean_can_shv   (recp)       )
            csite%mmean_can_dmol (recp) = idealdmolsh ( csite%mmean_can_prss  (recp)       &
                                                      , csite%mmean_can_temp  (recp)       &
                                                      , csite%mmean_can_shv   (recp)       )
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Find the soil mean temperature, liquid water fraction, and matric       !
            ! potential.                                                                   !
            !------------------------------------------------------------------------------!
            do iii=lsl,mzg
               nsoil = ntext_soil(iii)
               call uextcm2tl( csite%mmean_soil_energy(iii,recp)                           &
                             , csite%mmean_soil_water (iii,recp) * wdns                    &
                             , soil(nsoil)%slcpd                                           &
                             , csite%mmean_soil_temp  (iii,recp)                           &
                             , csite%mmean_soil_fliq  (iii,recp))

               csite%mmean_soil_mstpot(iii,recp) =                                         &
                                   matric_potential(nsoil,csite%mmean_soil_water(iii,recp))
            end do
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Find the temporary surface water properties.  They may not be available  !
            ! at all times, so we must check.                                              !
            !------------------------------------------------------------------------------!
            !----- Temporarily make energy extensive [J/m2]. ------------------------------!
            csite%mmean_sfcw_depth (recp) = csite%mmean_sfcw_depth (donp) * dawgt          &
                                          + csite%mmean_sfcw_depth (recp) * rawgt
            csite%mmean_sfcw_energy(recp) = csite%mmean_sfcw_energy(donp)                  &
                                          * csite%mmean_sfcw_mass  (donp) * dawgt          &
                                          + csite%mmean_sfcw_energy(recp)                  &
                                          * csite%mmean_sfcw_mass  (recp) * rawgt
            csite%mmean_sfcw_mass  (recp) = csite%mmean_sfcw_mass  (donp) * dawgt          &
                                          + csite%mmean_sfcw_mass  (recp) * rawgt
            csite%mmean_snowfac    (recp) = csite%mmean_snowfac    (donp) * dawgt          &
                                          + csite%mmean_snowfac    (recp) * rawgt
            !----- Check whether there is enough surface water. ---------------------------!
            if (csite%mmean_sfcw_mass(recp) > tiny_sfcwater_mass) then
               csite%mmean_sfcw_energy     (recp) =   csite%mmean_sfcw_energy     (recp)   &
                                                  /   csite%mmean_sfcw_mass       (recp)
               call uint2tl( csite%mmean_sfcw_energy(recp)                                 &
                           , csite%mmean_sfcw_temp  (recp)                                 &
                           , csite%mmean_sfcw_fliq  (recp) )
            else
               csite%mmean_sfcw_mass  (recp)  = 0.
               csite%mmean_sfcw_depth (recp)  = 0.
               csite%mmean_sfcw_energy(recp)  = 0.
               csite%mmean_sfcw_temp  (recp)  = csite%mmean_soil_temp(mzg,recp)
               csite%mmean_sfcw_fliq  (recp)  = csite%mmean_soil_fliq(mzg,recp)
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Mean diel.                                                                      !
      !------------------------------------------------------------------------------------!
      if (writing_dcyc .and. (.not. fuse_initial)) then
         if ( all(csite%qmean_can_prss > 10.0) ) then
            !------------------------------------------------------------------------------!
            !    First we solve the mean sum of squares as they depend on the mean and     !
            ! the original receptor data is lost after fusion takes place.                 !
            !------------------------------------------------------------------------------!
            do t=1,ndcycle
               csite%qmsqu_rh        (t,recp) = fuse_msqu( csite%qmean_rh        (t,donp)  &
                                                         , csite%qmsqu_rh        (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_rh        (t,recp)  &
                                                         , csite%qmsqu_rh        (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_fgc_rh    (t,recp) = fuse_msqu( csite%qmean_fgc_rh    (t,donp)  &
                                                         , csite%qmsqu_fgc_rh    (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_fgc_rh    (t,recp)  &
                                                         , csite%qmsqu_fgc_rh    (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_fsc_rh    (t,recp) = fuse_msqu( csite%qmean_fsc_rh    (t,donp)  &
                                                         , csite%qmsqu_fsc_rh    (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_fsc_rh    (t,recp)  &
                                                         , csite%qmsqu_fsc_rh    (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_stgc_rh   (t,recp) = fuse_msqu( csite%qmean_stgc_rh   (t,donp)  &
                                                         , csite%qmsqu_stgc_rh   (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_stgc_rh   (t,recp)  &
                                                         , csite%qmsqu_stgc_rh   (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_stsc_rh   (t,recp) = fuse_msqu( csite%qmean_stsc_rh   (t,donp)  &
                                                         , csite%qmsqu_stsc_rh   (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_stsc_rh   (t,recp)  &
                                                         , csite%qmsqu_stsc_rh   (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_msc_rh    (t,recp) = fuse_msqu( csite%qmean_msc_rh    (t,donp)  &
                                                         , csite%qmsqu_msc_rh    (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_msc_rh    (t,recp)  &
                                                         , csite%qmsqu_msc_rh    (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_ssc_rh    (t,recp) = fuse_msqu( csite%qmean_ssc_rh    (t,donp)  &
                                                         , csite%qmsqu_ssc_rh    (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_ssc_rh    (t,recp)  &
                                                         , csite%qmsqu_ssc_rh    (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_psc_rh    (t,recp) = fuse_msqu( csite%qmean_psc_rh    (t,donp)  &
                                                         , csite%qmsqu_psc_rh    (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_psc_rh    (t,recp)  &
                                                         , csite%qmsqu_psc_rh    (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_nep       (t,recp) = fuse_msqu( csite%qmean_nep       (t,donp)  &
                                                         , csite%qmsqu_nep       (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_nep       (t,recp)  &
                                                         , csite%qmsqu_nep       (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_rlongup   (t,recp) = fuse_msqu( csite%qmean_rlongup   (t,donp)  &
                                                         , csite%qmsqu_rlongup   (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_rlongup   (t,recp)  &
                                                         , csite%qmsqu_rlongup   (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_parup     (t,recp) = fuse_msqu( csite%qmean_parup     (t,donp)  &
                                                         , csite%qmsqu_parup     (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_parup     (t,recp)  &
                                                         , csite%qmsqu_parup     (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_nirup     (t,recp) = fuse_msqu( csite%qmean_nirup     (t,donp)  &
                                                         , csite%qmsqu_nirup     (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_nirup     (t,recp)  &
                                                         , csite%qmsqu_nirup     (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_rshortup  (t,recp) = fuse_msqu( csite%qmean_rshortup  (t,donp)  &
                                                         , csite%qmsqu_rshortup  (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_rshortup  (t,recp)  &
                                                         , csite%qmsqu_rshortup  (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_rnet      (t,recp) = fuse_msqu( csite%qmean_rnet      (t,donp)  &
                                                         , csite%qmsqu_rnet      (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_rnet      (t,recp)  &
                                                         , csite%qmsqu_rnet      (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_albedo    (t,recp) = fuse_msqu( csite%qmean_albedo    (t,donp)  &
                                                         , csite%qmsqu_albedo    (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_albedo    (t,recp)  &
                                                         , csite%qmsqu_albedo    (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_ustar     (t,recp) = fuse_msqu( csite%qmean_ustar     (t,donp)  &
                                                         , csite%qmsqu_ustar     (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_ustar     (t,recp)  &
                                                         , csite%qmsqu_ustar     (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_carbon_ac (t,recp) = fuse_msqu( csite%qmean_carbon_ac (t,donp)  &
                                                         , csite%qmsqu_carbon_ac (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_carbon_ac (t,recp)  &
                                                         , csite%qmsqu_carbon_ac (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_carbon_st (t,recp) = fuse_msqu( csite%qmean_carbon_st (t,donp)  &
                                                         , csite%qmsqu_carbon_st (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_carbon_st (t,recp)  &
                                                         , csite%qmsqu_carbon_st (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_vapor_gc  (t,recp) = fuse_msqu( csite%qmean_vapor_gc  (t,donp)  &
                                                         , csite%qmsqu_vapor_gc  (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_vapor_gc  (t,recp)  &
                                                         , csite%qmsqu_vapor_gc  (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_vapor_ac  (t,recp) = fuse_msqu( csite%qmean_vapor_ac  (t,donp)  &
                                                         , csite%qmsqu_vapor_ac  (t,donp)  &
                                                         , dawgt                           &
                                                         , csite%qmean_vapor_ac  (t,recp)  &
                                                         , csite%qmsqu_vapor_ac  (t,recp)  &
                                                         , rawgt                           &
                                                         , corr_patch, .false.)
               csite%qmsqu_sensible_gc(t,recp) =                                           &
                                               fuse_msqu( csite%qmean_sensible_gc(t,donp)  &
                                                        , csite%qmsqu_sensible_gc(t,donp)  &
                                                        , dawgt                            &
                                                        , csite%qmean_sensible_gc(t,recp)  &
                                                        , csite%qmsqu_sensible_gc(t,recp)  &
                                                        , rawgt                            &
                                                        , corr_patch, .false.)
               csite%qmsqu_sensible_ac(t,recp) =                                           &
                                               fuse_msqu( csite%qmean_sensible_ac(t,donp)  &
                                                        , csite%qmsqu_sensible_ac(t,donp)  &
                                                        , dawgt                            &
                                                        , csite%qmean_sensible_ac(t,recp)  &
                                                        , csite%qmsqu_sensible_ac(t,recp)  &
                                                        , rawgt                            &
                                                        , corr_patch, .false.)
            end do
            !------------------------------------------------------------------------------!


            csite%qmean_rh             (:,recp) = csite%qmean_rh          (:,donp) * dawgt &
                                                + csite%qmean_rh          (:,recp) * rawgt
            csite%qmean_fgc_rh         (:,recp) = csite%qmean_fgc_rh      (:,donp) * dawgt &
                                                + csite%qmean_fgc_rh      (:,recp) * rawgt
            csite%qmean_fsc_rh         (:,recp) = csite%qmean_fsc_rh      (:,donp) * dawgt &
                                                + csite%qmean_fsc_rh      (:,recp) * rawgt
            csite%qmean_stgc_rh        (:,recp) = csite%qmean_stgc_rh     (:,donp) * dawgt &
                                                + csite%qmean_stgc_rh     (:,recp) * rawgt
            csite%qmean_stsc_rh        (:,recp) = csite%qmean_stsc_rh     (:,donp) * dawgt &
                                                + csite%qmean_stsc_rh     (:,recp) * rawgt
            csite%qmean_msc_rh         (:,recp) = csite%qmean_msc_rh      (:,donp) * dawgt &
                                                + csite%qmean_msc_rh      (:,recp) * rawgt
            csite%qmean_ssc_rh         (:,recp) = csite%qmean_ssc_rh      (:,donp) * dawgt &
                                                + csite%qmean_ssc_rh      (:,recp) * rawgt
            csite%qmean_psc_rh         (:,recp) = csite%qmean_psc_rh      (:,donp) * dawgt &
                                                + csite%qmean_psc_rh      (:,recp) * rawgt
            csite%qmean_nep            (:,recp) = csite%qmean_nep         (:,donp) * dawgt &
                                                + csite%qmean_nep         (:,recp) * rawgt
            csite%qmean_rk4step        (:,recp) = csite%qmean_rk4step     (:,donp) * dawgt &
                                                + csite%qmean_rk4step     (:,recp) * rawgt
            csite%qmean_available_water(:,recp) = csite%qmean_available_water(:,donp)      &
                                                * dawgt                                    &
                                                + csite%qmean_available_water(:,recp)      &
                                                * rawgt
            csite%qmean_veg_displace   (:,recp) = csite%qmean_veg_displace(:,donp) * dawgt &
                                                + csite%qmean_veg_displace(:,recp) * rawgt
            csite%qmean_rough          (:,recp) = csite%qmean_rough       (:,donp) * dawgt &
                                                + csite%qmean_rough       (:,recp) * rawgt
            csite%qmean_can_theiv      (:,recp) = csite%qmean_can_theiv   (:,donp) * dawgt &
                                                + csite%qmean_can_theiv   (:,recp) * rawgt
            csite%qmean_can_theta      (:,recp) = csite%qmean_can_theta   (:,donp) * dawgt &
                                                + csite%qmean_can_theta   (:,recp) * rawgt
            csite%qmean_can_vpdef      (:,recp) = csite%qmean_can_vpdef   (:,donp) * dawgt &
                                                + csite%qmean_can_vpdef   (:,recp) * rawgt
            csite%qmean_can_shv        (:,recp) = csite%qmean_can_shv     (:,donp) * dawgt &
                                                + csite%qmean_can_shv     (:,recp) * rawgt
            csite%qmean_can_co2        (:,recp) = csite%qmean_can_co2     (:,donp) * dawgt &
                                                + csite%qmean_can_co2     (:,recp) * rawgt
            csite%qmean_can_prss       (:,recp) = csite%qmean_can_prss    (:,donp) * dawgt &
                                                + csite%qmean_can_prss    (:,recp) * rawgt
            csite%qmean_gnd_temp       (:,recp) = csite%qmean_gnd_temp    (:,donp) * dawgt &
                                                + csite%qmean_gnd_temp    (:,recp) * rawgt
            csite%qmean_gnd_shv        (:,recp) = csite%qmean_gnd_shv     (:,donp) * dawgt &
                                                + csite%qmean_gnd_shv     (:,recp) * rawgt
            csite%qmean_can_ggnd       (:,recp) = csite%qmean_can_ggnd    (:,donp) * dawgt &
                                                + csite%qmean_can_ggnd    (:,recp) * rawgt
            csite%qmean_rshort_gnd     (:,recp) = csite%qmean_rshort_gnd  (:,donp) * dawgt &
                                                + csite%qmean_rshort_gnd  (:,recp) * rawgt
            csite%qmean_par_gnd        (:,recp) = csite%qmean_par_gnd     (:,donp) * dawgt &
                                                + csite%qmean_par_gnd     (:,recp) * rawgt
            csite%qmean_rlong_gnd      (:,recp) = csite%qmean_rlong_gnd   (:,donp) * dawgt &
                                                + csite%qmean_rlong_gnd   (:,recp) * rawgt
            csite%qmean_rlongup        (:,recp) = csite%qmean_rlongup     (:,donp) * dawgt &
                                                + csite%qmean_rlongup     (:,recp) * rawgt
            csite%qmean_parup          (:,recp) = csite%qmean_parup       (:,donp) * dawgt &
                                                + csite%qmean_parup       (:,recp) * rawgt
            csite%qmean_nirup          (:,recp) = csite%qmean_nirup       (:,donp) * dawgt &
                                                + csite%qmean_nirup       (:,recp) * rawgt
            csite%qmean_rshortup       (:,recp) = csite%qmean_rshortup    (:,donp) * dawgt &
                                                + csite%qmean_rshortup    (:,recp) * rawgt
            csite%qmean_rnet           (:,recp) = csite%qmean_rnet        (:,donp) * dawgt &
                                                + csite%qmean_rnet        (:,recp) * rawgt
            csite%qmean_albedo         (:,recp) = csite%qmean_albedo      (:,donp) * dawgt &
                                                + csite%qmean_albedo      (:,recp) * rawgt
            csite%qmean_albedo_par     (:,recp) = csite%qmean_albedo_par  (:,donp) * dawgt &
                                                + csite%qmean_albedo_par  (:,recp) * rawgt
            csite%qmean_albedo_nir     (:,recp) = csite%qmean_albedo_nir  (:,donp) * dawgt &
                                                + csite%qmean_albedo_nir  (:,recp) * rawgt
            csite%qmean_rlong_albedo   (:,recp) = csite%qmean_rlong_albedo(:,donp) * dawgt &
                                                + csite%qmean_rlong_albedo(:,recp) * rawgt
            csite%qmean_ustar          (:,recp) = csite%qmean_ustar       (:,donp) * dawgt &
                                                + csite%qmean_ustar       (:,recp) * rawgt
            csite%qmean_tstar          (:,recp) = csite%qmean_tstar       (:,donp) * dawgt &
                                                + csite%qmean_tstar       (:,recp) * rawgt
            csite%qmean_qstar          (:,recp) = csite%qmean_qstar       (:,donp) * dawgt &
                                                + csite%qmean_qstar       (:,recp) * rawgt
            csite%qmean_cstar          (:,recp) = csite%qmean_cstar       (:,donp) * dawgt &
                                                + csite%qmean_cstar       (:,recp) * rawgt
            csite%qmean_carbon_ac      (:,recp) = csite%qmean_carbon_ac   (:,donp) * dawgt &
                                                + csite%qmean_carbon_ac   (:,recp) * rawgt
            csite%qmean_carbon_st      (:,recp) = csite%qmean_carbon_st   (:,donp) * dawgt &
                                                + csite%qmean_carbon_st   (:,recp) * rawgt
            csite%qmean_vapor_gc       (:,recp) = csite%qmean_vapor_gc    (:,donp) * dawgt &
                                                + csite%qmean_vapor_gc    (:,recp) * rawgt
            csite%qmean_vapor_ac       (:,recp) = csite%qmean_vapor_ac    (:,donp) * dawgt &
                                                + csite%qmean_vapor_ac    (:,recp) * rawgt
            csite%qmean_throughfall    (:,recp) = csite%qmean_throughfall (:,donp) * dawgt &
                                                + csite%qmean_throughfall (:,recp) * rawgt
            csite%qmean_runoff         (:,recp) = csite%qmean_runoff      (:,donp) * dawgt &
                                                + csite%qmean_runoff      (:,recp) * rawgt
            csite%qmean_drainage       (:,recp) = csite%qmean_drainage    (:,donp) * dawgt &
                                                + csite%qmean_drainage    (:,recp) * rawgt
            csite%qmean_sensible_gc    (:,recp) = csite%qmean_sensible_gc (:,donp) * dawgt &
                                                + csite%qmean_sensible_gc (:,recp) * rawgt
            csite%qmean_sensible_ac    (:,recp) = csite%qmean_sensible_ac (:,donp) * dawgt &
                                                + csite%qmean_sensible_ac (:,recp) * rawgt
            csite%qmean_qthroughfall   (:,recp) = csite%qmean_qthroughfall(:,donp) * dawgt &
                                                + csite%qmean_qthroughfall(:,recp) * rawgt
            csite%qmean_qrunoff        (:,recp) = csite%qmean_qrunoff     (:,donp) * dawgt &
                                                + csite%qmean_qrunoff     (:,recp) * rawgt
            csite%qmean_qdrainage      (:,recp) = csite%qmean_qdrainage   (:,donp) * dawgt &
                                                + csite%qmean_qdrainage   (:,recp) * rawgt
            csite%qmean_smoist_gg    (:,:,recp) = csite%qmean_smoist_gg  (:,:,donp)* dawgt &
                                                + csite%qmean_smoist_gg  (:,:,recp)* rawgt
            csite%qmean_transloss    (:,:,recp) = csite%qmean_transloss  (:,:,donp)* dawgt &
                                                + csite%qmean_transloss  (:,:,recp)* rawgt
            csite%qmean_sensible_gg  (:,:,recp) = csite%qmean_sensible_gg(:,:,donp)* dawgt &
                                                + csite%qmean_sensible_gg(:,:,recp)* rawgt


            do t=1,ndcycle
               !---------------------------------------------------------------------------!
               !      Now we find the derived properties for the canopy air space.         !
               !---------------------------------------------------------------------------!
               xmean_can_exner               = press2exner ( csite%qmean_can_prss (t,recp) )
               csite%qmean_can_temp (t,recp) = extheta2temp( xmean_can_exner               &
                                                           , csite%qmean_can_theta(t,recp) )
               csite%qmean_can_rhos (t,recp) = idealdenssh ( csite%qmean_can_prss (t,recp) &
                                                           , csite%qmean_can_temp (t,recp) &
                                                           , csite%qmean_can_shv  (t,recp) )
               csite%qmean_can_dmol (t,recp) = idealdmolsh ( csite%qmean_can_prss (t,recp) &
                                                           , csite%qmean_can_temp (t,recp) &
                                                           , csite%qmean_can_shv  (t,recp) )
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Find the soil mean temperature, liquid water fraction, and matric    !
               ! potential.                                                                !
               !---------------------------------------------------------------------------!
               do iii=lsl,mzg
                  nsoil = ntext_soil(iii)
                  call uextcm2tl( csite%qmean_soil_energy(iii,t,recp)                      &
                                , csite%qmean_soil_water (iii,t,recp) * wdns               &
                                , soil(nsoil)%slcpd                                        &
                                , csite%qmean_soil_temp  (iii,t,recp)                      &
                                , csite%qmean_soil_fliq  (iii,t,recp) )

                  csite%qmean_soil_mstpot(iii,t,recp) =                                    &
                                 matric_potential(nsoil,csite%qmean_soil_water(iii,t,recp))
               end do
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Find the temporary surface water properties.  They may not be         !
               ! available at all times, so we must check.                                 !
               !---------------------------------------------------------------------------!
               !----- Temporarily make energy extensive [J/m2]. ---------------------------!
               csite%qmean_sfcw_depth (t,recp) = csite%qmean_sfcw_depth (t,donp) * dawgt   &
                                               + csite%qmean_sfcw_depth (t,recp) * rawgt
               csite%qmean_sfcw_energy(t,recp) = csite%qmean_sfcw_energy(t,donp)           &
                                               * csite%qmean_sfcw_mass  (t,donp) * dawgt   &
                                               + csite%qmean_sfcw_energy(t,recp)           &
                                               * csite%qmean_sfcw_mass  (t,recp) * rawgt
               csite%qmean_sfcw_mass  (t,recp) = csite%qmean_sfcw_mass  (t,donp) * dawgt   &
                                               + csite%qmean_sfcw_mass  (t,recp) * rawgt
               csite%qmean_snowfac    (t,recp) = csite%qmean_snowfac    (t,donp) * dawgt   &
                                               + csite%qmean_snowfac    (t,recp) * rawgt
               !----- Check whether there is enough surface water. ------------------------!
               if (csite%qmean_sfcw_mass(t,recp) > tiny_sfcwater_mass) then
                  csite%qmean_sfcw_energy   (t,recp) =   csite%qmean_sfcw_energy(t,recp)   &
                                                     /   csite%qmean_sfcw_mass  (t,recp)
                  call uint2tl( csite%qmean_sfcw_energy(t,recp)                            &
                              , csite%qmean_sfcw_temp  (t,recp)                            &
                              , csite%qmean_sfcw_fliq  (t,recp) )
               else
                  csite%qmean_sfcw_mass  (t,recp)  = 0.
                  csite%qmean_sfcw_depth (t,recp)  = 0.
                  csite%qmean_sfcw_energy(t,recp)  = 0.
                  csite%qmean_sfcw_temp  (t,recp)  = csite%qmean_soil_temp(mzg,t,recp)
                  csite%qmean_sfcw_fliq  (t,recp)  = csite%qmean_soil_fliq(mzg,t,recp)
               end if
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Now we need to adjust the densities of cohorts. Because the patch area         !
      ! increased we want to retain the same total amount of mass and energy.              !
      !------------------------------------------------------------------------------------!
      !----- 1. Adjust densities of cohorts in donor patch --------------------------------!
      cpatch => csite%patch(donp)
      ndc = cpatch%ncohorts
      call update_cohort_extensive_props(cpatch,1,ndc,dawgt)
      !----- 2. Adjust densities of cohorts in recipient patch ----------------------------!
      cpatch => csite%patch(recp)
      nrc = cpatch%ncohorts
      call update_cohort_extensive_props(cpatch,1,nrc,rawgt)
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
            select case (ifusion)
            case (0)
               call old_fuse_cohorts(csite,recp,lsl,fuse_initial)
            case (1)
               call new_fuse_cohorts(csite,recp,lsl,fuse_initial)
            end select
            call terminate_cohorts(csite,recp,cmet,fuse_initial,elim_nplant,elim_lai)
            call split_cohorts(csite,recp,green_leaf_factor,fuse_initial)
         end if
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    Now we update some variables that depend on cohort statistics, namely:          !
      ! + csite%veg_height(recp)                                                           !
      ! + csite%veg_displace(recp)                                                         !
      ! + csite%disp_height(recp)                                                          !
      ! + csite%veg_rough(recp)                                                            !
      ! + csite%total_sfcw_depth(recp)                                                     !
      ! + csite%snowfac(recp)                                                              !
      ! + csite%opencan_frac(recp)                                                         !
      !                                                                                    !
      !     Even though the canopy air space depth effect is fused, there may be an        !
      ! additional residual effect after cohort fusion/fission/termination.  Just to be    !
      ! safe, we also account for these changes.  The only time we do not check it is      !
      ! during initialisation.                                                             !
      !------------------------------------------------------------------------------------!
      call update_patch_derived_props(csite,recp,.not. fuse_initial)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !    This subroutine will update the size profile within patch.                      !
      ! + csite%cumlai_profile(:,:,recp)                                                   !
      !------------------------------------------------------------------------------------!
      call patch_pft_size_profile(csite,recp)
      !------------------------------------------------------------------------------------!


      !----- Update the patch area. -------------------------------------------------------!
      csite%area(recp) = newarea
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Last, we make sure that fusion preserved energy, water, and carbon storage.    !
      !------------------------------------------------------------------------------------!
      if (checkbudget .and. (.not. fuse_initial)) then
         !----- Call a separate routine that performs the check. --------------------------!
         call check_bfusion_patch(csite,osite,lsl,recp)
         !---------------------------------------------------------------------------------!

         !----- In case everything looks good, proceed, but free memory first. ------------!
         call deallocate_sitetype(osite)
         deallocate(osite)
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!

      return

   end subroutine fuse_2_patches
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !  SUBROUTINE: FUSE_MSQU       
   !> \brief This subroutine combines the mean sum of squares of two quantities (x and y).
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






   !=======================================================================================!
   !=======================================================================================!
   !     This sub-routine checks that energy, water, carbon, and carbon dioxide were       !
   ! conserved after fusion.                                                               !
   !---------------------------------------------------------------------------------------!
   subroutine check_bfusion_patch(csite,osite,lsl,ipa)
      use ed_state_vars, only : sitetype                 ! ! structure
      use budget_utils , only : tol_subday_budget        & ! structure
                              , tol_carbon_budget        & ! structure
                              , compute_co2_storage      & ! function
                              , compute_carbon_storage   & ! function
                              , compute_water_storage    & ! function
                              , compute_enthalpy_storage ! ! function
      use ed_misc_coms , only : current_time             & ! intent(in)
                              , frqsum                   ! ! intent(in)
      use consts_coms  , only : tiny_num                 ! ! intent(in)
      use therm_lib    , only : tq2enthalpy              ! ! function
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype), target       :: csite
      type(sitetype), target       :: osite
      integer       , intent(in)   :: lsl
      integer       , intent(in)   :: ipa
      !----- Local variables. -------------------------------------------------------------!
      integer                      :: opa
      logical                      :: co2fused_fine
      logical                      :: cfused_fine
      logical                      :: wfused_fine
      logical                      :: efused_fine
      logical                      :: budget_fine
      real                         :: rawgt
      real                         :: dawgt
      real                         :: newareai
      real                         :: co2fused_initialstorage
      real                         :: cfused_initialstorage
      real                         :: wfused_initialstorage
      real                         :: efused_initialstorage
      real                         :: co2fused_residual
      real                         :: cfused_residual
      real                         :: wfused_residual
      real                         :: efused_residual
      real                         :: co2fused_tolerance
      real                         :: cfused_tolerance
      real                         :: wfused_tolerance
      real                         :: efused_tolerance
      real                         :: delta_co2_vegtdyn
      real                         :: delta_c_vegtdyn
      real                         :: delta_e_vegtdyn
      real                         :: delta_w_vegtdyn
      real, dimension(6)           :: co2_cas
      real, dimension(6)           :: water_soil
      real, dimension(6)           :: water_sfcwater
      real, dimension(6)           :: water_cas
      real, dimension(6)           :: water_leaf
      real, dimension(6)           :: water_wood
      real, dimension(6)           :: enthalpy_soil
      real, dimension(6)           :: enthalpy_sfcwater
      real, dimension(6)           :: enthalpy_cas
      real, dimension(6)           :: enthalpy_leaf
      real, dimension(6)           :: enthalpy_wood
      real, dimension(6)           :: carbon_necro
      real, dimension(6)           :: carbon_repro
      real, dimension(6)           :: carbon_balive
      real, dimension(6)           :: carbon_bdead
      real, dimension(6)           :: carbon_bstorage
      real, dimension(6)           :: carbon_cas
      real, dimension(6)           :: carbon_committed
      real, dimension(6)           :: co2_loss2atm
      real, dimension(6)           :: co2_denseffect
      real, dimension(6)           :: co2_zcaneffect
      real, dimension(6)           :: co2_gpp
      real, dimension(6)           :: co2_plresp
      real, dimension(6)           :: co2_rh
      real, dimension(6)           :: carbon_loss2atm
      real, dimension(6)           :: carbon_denseffect
      real, dimension(6)           :: carbon_zcaneffect
      real, dimension(6)           :: carbon_seedrain
      real, dimension(6)           :: carbon_loss2yield
      real, dimension(6)           :: water_loss2atm
      real, dimension(6)           :: water_denseffect
      real, dimension(6)           :: water_wcapeffect
      real, dimension(6)           :: water_zcaneffect
      real, dimension(6)           :: water_pheneffect
      real, dimension(6)           :: water_precipgain
      real, dimension(6)           :: water_loss2runoff
      real, dimension(6)           :: water_loss2drainage
      real, dimension(6)           :: enthalpy_loss2atm
      real, dimension(6)           :: enthalpy_denseffect
      real, dimension(6)           :: enthalpy_prsseffect
      real, dimension(6)           :: enthalpy_hcapeffect
      real, dimension(6)           :: enthalpy_wcapeffect
      real, dimension(6)           :: enthalpy_zcaneffect
      real, dimension(6)           :: enthalpy_pheneffect
      real, dimension(6)           :: enthalpy_loss2runoff
      real, dimension(6)           :: enthalpy_loss2drainage
      real, dimension(6)           :: enthalpy_netrad
      real, dimension(6)           :: enthalpy_precipgain
      real, dimension(3)           :: can_prss
      real, dimension(3)           :: can_temp
      real, dimension(3)           :: can_shv
      real, dimension(3)           :: can_co2
      real, dimension(3)           :: can_rhos
      real, dimension(3)           :: can_dmol
      real, dimension(3)           :: can_enthalpy
      real, dimension(3)           :: can_depth
      !----- Local constants. -------------------------------------------------------------!
      character(len=10), parameter :: fmtl='(a,13x,l1)'
      character(len=10), parameter :: fmti='(a,1x,i14)'
      character(len=11), parameter :: fmth='(a,6(1x,a))'
      character(len=11), parameter :: fmtx='(a,3(1x,a))'
      character(len=13), parameter :: fmtf='(a,1x,es14.7)'
      character(len=16), parameter :: fmtd='(a,6(1x,es14.7))'
      character(len=16), parameter :: fmts='(a,3(1x,es14.7))'
      character(len=27), parameter :: fmtt='(a,i4.4,2(1x,i2.2),1x,f6.0)'
      !------------------------------------------------------------------------------------!



      !----- Compute current storage terms. -----------------------------------------------!
      co2fused_initialstorage = compute_co2_storage     (csite,ipa)
      cfused_initialstorage   = compute_carbon_storage  (csite,ipa,0)
      efused_initialstorage   = compute_enthalpy_storage(csite,lsl,ipa,0)
      wfused_initialstorage   = compute_water_storage   (csite,lsl,ipa,0)
      !------------------------------------------------------------------------------------!



      !----- Define tolerance. ------------------------------------------------------------!
      co2fused_tolerance      = csite%co2budget_initialstorage(ipa) * tol_subday_budget
      cfused_tolerance        = csite%cbudget_initialstorage  (ipa) * tol_carbon_budget
      efused_tolerance        = csite%ebudget_initialstorage  (ipa) * tol_subday_budget
      wfused_tolerance        = csite%wbudget_initialstorage  (ipa) * tol_subday_budget
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Find the long-term dynamic effects that have not yet been incorporated to the !
      ! storage (they will after the first thermodynamic step).                            !
      !------------------------------------------------------------------------------------!
      delta_co2_vegtdyn = frqsum * csite%co2budget_zcaneffect(ipa)
      delta_c_vegtdyn   = frqsum                                                           &
                        * ( csite%cbudget_loss2yield(ipa) + csite%cbudget_seedrain  (ipa)  &
                          + csite%cbudget_zcaneffect(ipa) )
      delta_e_vegtdyn   = frqsum                                                           &
                        * ( csite%ebudget_hcapeffect(ipa) + csite%ebudget_wcapeffect(ipa)  &
                          + csite%ebudget_zcaneffect(ipa) + csite%ebudget_pheneffect(ipa) )
      delta_w_vegtdyn   = frqsum                                                           &
                        * ( csite%wbudget_wcapeffect(ipa) + csite%wbudget_zcaneffect(ipa)  &
                          + csite%wbudget_pheneffect(ipa) )
      !------------------------------------------------------------------------------------!



      !----- Find residuals. --------------------------------------------------------------!
      co2fused_residual = csite%co2budget_initialstorage(ipa) + delta_co2_vegtdyn          &
                        - co2fused_initialstorage
      cfused_residual   = csite%cbudget_initialstorage  (ipa) + delta_c_vegtdyn            &
                        - cfused_initialstorage
      efused_residual   = csite%ebudget_initialstorage  (ipa) + delta_e_vegtdyn            &
                        - efused_initialstorage
      wfused_residual   = csite%wbudget_initialstorage  (ipa) + delta_w_vegtdyn            &
                        - wfused_initialstorage
      !------------------------------------------------------------------------------------!



      !----- Check that everything is within tolerance. -----------------------------------!
      co2fused_fine = abs(co2fused_residual) <= co2fused_tolerance
      cfused_fine   = abs(cfused_residual  ) <= cfused_tolerance
      efused_fine   = abs(efused_residual  ) <= efused_tolerance
      wfused_fine   = abs(wfused_residual  ) <= wfused_tolerance
      budget_fine   = co2fused_fine .and. cfused_fine .and. wfused_fine .and. efused_fine
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     In case anything went wrong, crash the entire simulation, and report the       !
      ! information to the end-user.                                                       !
      !------------------------------------------------------------------------------------!
      if (.not. budget_fine) then
         !----- Find the relative area weights. -------------------------------------------!
         newareai = 1. / ( osite%area(1) + osite%area(2))
         rawgt    = osite%area(1) * newareai
         dawgt    = osite%area(2) * newareai
         !---------------------------------------------------------------------------------!


         !----- Find the storage and budget of particular pools (receptor and donor). -----!
         do opa = 1,2
            !----- Storage terms. ---------------------------------------------------------!
            co2_cas               (opa) = compute_co2_storage     (osite,opa)
            water_soil            (opa) = compute_water_storage   (osite,lsl,opa,1)
            water_sfcwater        (opa) = compute_water_storage   (osite,lsl,opa,2)
            water_cas             (opa) = compute_water_storage   (osite,lsl,opa,3)
            water_leaf            (opa) = compute_water_storage   (osite,lsl,opa,4)
            water_wood            (opa) = compute_water_storage   (osite,lsl,opa,5)
            enthalpy_soil         (opa) = compute_enthalpy_storage(osite,lsl,opa,1)
            enthalpy_sfcwater     (opa) = compute_enthalpy_storage(osite,lsl,opa,2)
            enthalpy_cas          (opa) = compute_enthalpy_storage(osite,lsl,opa,3)
            enthalpy_leaf         (opa) = compute_enthalpy_storage(osite,lsl,opa,4)
            enthalpy_wood         (opa) = compute_enthalpy_storage(osite,lsl,opa,5)
            carbon_necro          (opa) = compute_carbon_storage  (osite,opa,1)
            carbon_repro          (opa) = compute_carbon_storage  (osite,opa,2)
            carbon_balive         (opa) = compute_carbon_storage  (osite,opa,3)
            carbon_bdead          (opa) = compute_carbon_storage  (osite,opa,4)
            carbon_bstorage       (opa) = compute_carbon_storage  (osite,opa,5)
            carbon_cas            (opa) = compute_carbon_storage  (osite,opa,6)
            carbon_committed      (opa) = compute_carbon_storage  (osite,opa,7)
            !------ Budget terms. ---------------------------------------------------------!
            co2_loss2atm          (opa) = osite%co2budget_loss2atm   (opa)
            co2_denseffect        (opa) = osite%co2budget_denseffect (opa)
            co2_zcaneffect        (opa) = osite%co2budget_zcaneffect (opa)
            co2_gpp               (opa) = osite%co2budget_gpp        (opa)
            co2_plresp            (opa) = osite%co2budget_plresp     (opa)
            co2_rh                (opa) = osite%co2budget_rh         (opa)
            carbon_loss2atm       (opa) = osite%cbudget_loss2atm     (opa)
            carbon_denseffect     (opa) = osite%cbudget_denseffect   (opa)
            carbon_zcaneffect     (opa) = osite%cbudget_zcaneffect   (opa)
            carbon_seedrain       (opa) = osite%cbudget_seedrain     (opa)
            carbon_loss2yield     (opa) = osite%cbudget_loss2yield   (opa)
            water_loss2atm        (opa) = osite%wbudget_loss2atm     (opa)
            water_denseffect      (opa) = osite%wbudget_denseffect   (opa)
            water_wcapeffect      (opa) = osite%wbudget_wcapeffect   (opa)
            water_zcaneffect      (opa) = osite%wbudget_zcaneffect   (opa)
            water_pheneffect      (opa) = osite%wbudget_pheneffect   (opa)
            water_precipgain      (opa) = osite%wbudget_precipgain   (opa)
            water_loss2runoff     (opa) = osite%wbudget_loss2runoff  (opa)
            water_loss2drainage   (opa) = osite%wbudget_loss2drainage(opa)
            enthalpy_loss2atm     (opa) = osite%ebudget_loss2atm     (opa)
            enthalpy_denseffect   (opa) = osite%ebudget_denseffect   (opa)
            enthalpy_prsseffect   (opa) = osite%ebudget_prsseffect   (opa)
            enthalpy_hcapeffect   (opa) = osite%ebudget_hcapeffect   (opa)
            enthalpy_wcapeffect   (opa) = osite%ebudget_wcapeffect   (opa)
            enthalpy_zcaneffect   (opa) = osite%ebudget_zcaneffect   (opa)
            enthalpy_pheneffect   (opa) = osite%ebudget_pheneffect   (opa)
            enthalpy_loss2runoff  (opa) = osite%ebudget_loss2runoff  (opa)
            enthalpy_loss2drainage(opa) = osite%ebudget_loss2drainage(opa)
            enthalpy_netrad       (opa) = osite%ebudget_netrad       (opa)
            enthalpy_precipgain   (opa) = osite%ebudget_precipgain   (opa)
            !------ Canopy air space properties. ------------------------------------------!
            can_prss    (opa) = osite%can_prss    (opa)
            can_temp    (opa) = osite%can_temp    (opa)
            can_shv     (opa) = osite%can_shv     (opa)
            can_co2     (opa) = osite%can_co2     (opa)
            can_rhos    (opa) = osite%can_rhos    (opa)
            can_dmol    (opa) = osite%can_dmol    (opa)
            can_enthalpy(opa) = tq2enthalpy(osite%can_temp(opa),osite%can_shv(opa),.true.)
            can_depth   (opa) = osite%can_depth   (opa)
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Find the storage and budget of fused pools ("bottom-up").                   !
         !---------------------------------------------------------------------------------!
         !----- Storage terms. ------------------------------------------------------------!
         co2_cas               (3) = co2_cas               (1) * rawgt                     &
                                   + co2_cas               (2) * dawgt
         water_soil            (3) = water_soil            (1) * rawgt                     &
                                   + water_soil            (2) * dawgt
         water_sfcwater        (3) = water_sfcwater        (1) * rawgt                     &
                                   + water_sfcwater        (2) * dawgt
         water_cas             (3) = water_cas             (1) * rawgt                     &
                                   + water_cas             (2) * dawgt
         water_leaf            (3) = water_leaf            (1) * rawgt                     &
                                   + water_leaf            (2) * dawgt
         water_wood            (3) = water_wood            (1) * rawgt                     &
                                   + water_wood            (2) * dawgt
         enthalpy_soil         (3) = enthalpy_soil         (1) * rawgt                     &
                                   + enthalpy_soil         (2) * dawgt
         enthalpy_sfcwater     (3) = enthalpy_sfcwater     (1) * rawgt                     &
                                   + enthalpy_sfcwater     (2) * dawgt
         enthalpy_cas          (3) = enthalpy_cas          (1) * rawgt                     &
                                   + enthalpy_cas          (2) * dawgt
         enthalpy_leaf         (3) = enthalpy_leaf         (1) * rawgt                     &
                                   + enthalpy_leaf         (2) * dawgt
         enthalpy_wood         (3) = enthalpy_wood         (1) * rawgt                     &
                                   + enthalpy_wood         (2) * dawgt
         carbon_necro          (3) = carbon_necro          (1) * rawgt                     &
                                   + carbon_necro          (2) * dawgt
         carbon_repro          (3) = carbon_repro          (1) * rawgt                     &
                                   + carbon_repro          (2) * dawgt
         carbon_balive         (3) = carbon_balive         (1) * rawgt                     &
                                   + carbon_balive         (2) * dawgt
         carbon_bdead          (3) = carbon_bdead          (1) * rawgt                     &
                                   + carbon_bdead          (2) * dawgt
         carbon_bstorage       (3) = carbon_bstorage       (1) * rawgt                     &
                                   + carbon_bstorage       (2) * dawgt
         carbon_cas            (3) = carbon_cas            (1) * rawgt                     &
                                   + carbon_cas            (2) * dawgt
         carbon_committed      (3) = carbon_committed      (1) * rawgt                     &
                                   + carbon_committed      (2) * dawgt
         !------ Budget terms. ------------------------------------------------------------!
         co2_loss2atm          (3) = co2_loss2atm          (1) * rawgt                     &
                                   + co2_loss2atm          (2) * dawgt
         co2_denseffect        (3) = co2_denseffect        (1) * rawgt                     &
                                   + co2_denseffect        (2) * dawgt
         co2_zcaneffect        (3) = co2_zcaneffect        (1) * rawgt                     &
                                   + co2_zcaneffect        (2) * dawgt
         co2_gpp               (3) = co2_gpp               (1) * rawgt                     &
                                   + co2_gpp               (2) * dawgt
         co2_plresp            (3) = co2_plresp            (1) * rawgt                     &
                                   + co2_plresp            (2) * dawgt
         co2_rh                (3) = co2_rh                (1) * rawgt                     &
                                   + co2_rh                (2) * dawgt
         carbon_loss2atm       (3) = carbon_loss2atm       (1) * rawgt                     &
                                   + carbon_loss2atm       (2) * dawgt
         carbon_denseffect     (3) = carbon_denseffect     (1) * rawgt                     &
                                   + carbon_denseffect     (2) * dawgt
         carbon_zcaneffect     (3) = carbon_zcaneffect     (1) * rawgt                     &
                                   + carbon_zcaneffect     (2) * dawgt
         carbon_seedrain       (3) = carbon_seedrain       (1) * rawgt                     &
                                   + carbon_seedrain       (2) * dawgt
         carbon_loss2yield     (3) = carbon_loss2yield     (1) * rawgt                     &
                                   + carbon_loss2yield     (2) * dawgt
         water_loss2atm        (3) = water_loss2atm        (1) * rawgt                     &
                                   + water_loss2atm        (2) * dawgt
         water_denseffect      (3) = water_denseffect      (1) * rawgt                     &
                                   + water_denseffect      (2) * dawgt
         water_wcapeffect      (3) = water_wcapeffect      (1) * rawgt                     &
                                   + water_wcapeffect      (2) * dawgt
         water_zcaneffect      (3) = water_zcaneffect      (1) * rawgt                     &
                                   + water_zcaneffect      (2) * dawgt
         water_pheneffect      (3) = water_pheneffect      (1) * rawgt                     &
                                   + water_pheneffect      (2) * dawgt
         water_precipgain      (3) = water_precipgain      (1) * rawgt                     &
                                   + water_precipgain      (2) * dawgt
         water_loss2runoff     (3) = water_loss2runoff     (1) * rawgt                     &
                                   + water_loss2runoff     (2) * dawgt
         water_loss2drainage   (3) = water_loss2drainage   (1) * rawgt                     &
                                   + water_loss2drainage   (2) * dawgt
         enthalpy_loss2atm     (3) = enthalpy_loss2atm     (1) * rawgt                     &
                                   + enthalpy_loss2atm     (2) * dawgt
         enthalpy_denseffect   (3) = enthalpy_denseffect   (1) * rawgt                     &
                                   + enthalpy_denseffect   (2) * dawgt
         enthalpy_prsseffect   (3) = enthalpy_prsseffect   (1) * rawgt                     &
                                   + enthalpy_prsseffect   (2) * dawgt
         enthalpy_hcapeffect   (3) = enthalpy_hcapeffect   (1) * rawgt                     &
                                   + enthalpy_hcapeffect   (2) * dawgt
         enthalpy_wcapeffect   (3) = enthalpy_wcapeffect   (1) * rawgt                     &
                                   + enthalpy_wcapeffect   (2) * dawgt
         enthalpy_zcaneffect   (3) = enthalpy_zcaneffect   (1) * rawgt                     &
                                   + enthalpy_zcaneffect   (2) * dawgt
         enthalpy_pheneffect   (3) = enthalpy_pheneffect   (1) * rawgt                     &
                                   + enthalpy_pheneffect   (2) * dawgt
         enthalpy_loss2runoff  (3) = enthalpy_loss2runoff  (1) * rawgt                     &
                                   + enthalpy_loss2runoff  (2) * dawgt
         enthalpy_loss2drainage(3) = enthalpy_loss2drainage(1) * rawgt                     &
                                   + enthalpy_loss2drainage(2) * dawgt
         enthalpy_netrad       (3) = enthalpy_netrad       (1) * rawgt                     &
                                   + enthalpy_netrad       (2) * dawgt
         enthalpy_precipgain   (3) = enthalpy_precipgain   (1) * rawgt                     &
                                   + enthalpy_precipgain   (2) * dawgt
         !---------------------------------------------------------------------------------!


         !------ Canopy air space properties. ---------------------------------------------!
         can_prss    (3) = csite%can_prss    (ipa)
         can_temp    (3) = csite%can_temp    (ipa)
         can_shv     (3) = csite%can_shv     (ipa)
         can_co2     (3) = csite%can_co2     (ipa)
         can_rhos    (3) = csite%can_rhos    (ipa)
         can_dmol    (3) = csite%can_dmol    (ipa)
         can_enthalpy(3) = tq2enthalpy(csite%can_temp(ipa),csite%can_shv(ipa),.true.)
         can_depth   (3) = csite%can_depth   (ipa)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the storage and budget of fused pools ("top-down").                    !
         !---------------------------------------------------------------------------------!
         !----- Storage terms. ------------------------------------------------------------!
         co2_cas               (4) = compute_co2_storage     (csite,ipa)
         water_soil            (4) = compute_water_storage   (csite,lsl,ipa,1)
         water_sfcwater        (4) = compute_water_storage   (csite,lsl,ipa,2)
         water_cas             (4) = compute_water_storage   (csite,lsl,ipa,3)
         water_leaf            (4) = compute_water_storage   (csite,lsl,ipa,4)
         water_wood            (4) = compute_water_storage   (csite,lsl,ipa,5)
         enthalpy_soil         (4) = compute_enthalpy_storage(csite,lsl,ipa,1)
         enthalpy_sfcwater     (4) = compute_enthalpy_storage(csite,lsl,ipa,2)
         enthalpy_cas          (4) = compute_enthalpy_storage(csite,lsl,ipa,3)
         enthalpy_leaf         (4) = compute_enthalpy_storage(csite,lsl,ipa,4)
         enthalpy_wood         (4) = compute_enthalpy_storage(csite,lsl,ipa,5)
         carbon_necro          (4) = compute_carbon_storage  (csite,ipa,1)
         carbon_repro          (4) = compute_carbon_storage  (csite,ipa,2)
         carbon_balive         (4) = compute_carbon_storage  (csite,ipa,3)
         carbon_bdead          (4) = compute_carbon_storage  (csite,ipa,4)
         carbon_bstorage       (4) = compute_carbon_storage  (csite,ipa,5)
         carbon_cas            (4) = compute_carbon_storage  (csite,ipa,6)
         carbon_committed      (4) = compute_carbon_storage  (csite,ipa,7)
         !------ Budget terms. ------------------------------------------------------------!
         co2_loss2atm          (4) = csite%co2budget_loss2atm   (ipa)
         co2_denseffect        (4) = csite%co2budget_denseffect (ipa)
         co2_zcaneffect        (4) = csite%co2budget_zcaneffect (ipa)
         co2_gpp               (4) = csite%co2budget_gpp        (ipa)
         co2_plresp            (4) = csite%co2budget_plresp     (ipa)
         co2_rh                (4) = csite%co2budget_rh         (ipa)
         carbon_loss2atm       (4) = csite%cbudget_loss2atm     (ipa)
         carbon_denseffect     (4) = csite%cbudget_denseffect   (ipa)
         carbon_zcaneffect     (4) = csite%cbudget_zcaneffect   (ipa)
         carbon_seedrain       (4) = csite%cbudget_seedrain     (ipa)
         carbon_loss2yield     (4) = csite%cbudget_loss2yield   (ipa)
         water_loss2atm        (4) = csite%wbudget_loss2atm     (ipa)
         water_denseffect      (4) = csite%wbudget_denseffect   (ipa)
         water_wcapeffect      (4) = csite%wbudget_wcapeffect   (ipa)
         water_zcaneffect      (4) = csite%wbudget_zcaneffect   (ipa)
         water_pheneffect      (4) = csite%wbudget_pheneffect   (ipa)
         water_precipgain      (4) = csite%wbudget_precipgain   (ipa)
         water_loss2runoff     (4) = csite%wbudget_loss2runoff  (ipa)
         water_loss2drainage   (4) = csite%wbudget_loss2drainage(ipa)
         enthalpy_loss2atm     (4) = csite%ebudget_loss2atm     (ipa)
         enthalpy_denseffect   (4) = csite%ebudget_denseffect   (ipa)
         enthalpy_prsseffect   (4) = csite%ebudget_prsseffect   (ipa)
         enthalpy_hcapeffect   (4) = csite%ebudget_hcapeffect   (ipa)
         enthalpy_wcapeffect   (4) = csite%ebudget_wcapeffect   (ipa)
         enthalpy_zcaneffect   (4) = csite%ebudget_zcaneffect   (ipa)
         enthalpy_pheneffect   (4) = csite%ebudget_pheneffect   (ipa)
         enthalpy_loss2runoff  (4) = csite%ebudget_loss2runoff  (ipa)
         enthalpy_loss2drainage(4) = csite%ebudget_loss2drainage(ipa)
         enthalpy_netrad       (4) = csite%ebudget_netrad       (ipa)
         enthalpy_precipgain   (4) = csite%ebudget_precipgain   (ipa)
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Find the absolute difference between the two approaches.                    !
         !---------------------------------------------------------------------------------!
         !----- Storage terms. ------------------------------------------------------------!
         co2_cas               (5) = co2_cas               (3) - co2_cas               (4)
         water_soil            (5) = water_soil            (3) - water_soil            (4)
         water_sfcwater        (5) = water_sfcwater        (3) - water_sfcwater        (4)
         water_cas             (5) = water_cas             (3) - water_cas             (4)
         water_leaf            (5) = water_leaf            (3) - water_leaf            (4)
         water_wood            (5) = water_wood            (3) - water_wood            (4)
         enthalpy_soil         (5) = enthalpy_soil         (3) - enthalpy_soil         (4)
         enthalpy_sfcwater     (5) = enthalpy_sfcwater     (3) - enthalpy_sfcwater     (4)
         enthalpy_cas          (5) = enthalpy_cas          (3) - enthalpy_cas          (4)
         enthalpy_leaf         (5) = enthalpy_leaf         (3) - enthalpy_leaf         (4)
         enthalpy_wood         (5) = enthalpy_wood         (3) - enthalpy_wood         (4)
         carbon_necro          (5) = carbon_necro          (3) - carbon_necro          (4)
         carbon_repro          (5) = carbon_repro          (3) - carbon_repro          (4)
         carbon_balive         (5) = carbon_balive         (3) - carbon_balive         (4)
         carbon_bdead          (5) = carbon_bdead          (3) - carbon_bdead          (4)
         carbon_bstorage       (5) = carbon_bstorage       (3) - carbon_bstorage       (4)
         carbon_cas            (5) = carbon_cas            (3) - carbon_cas            (4)
         carbon_committed      (5) = carbon_committed      (3) - carbon_committed      (4)
         !------ Budget terms. ------------------------------------------------------------!
         co2_loss2atm          (5) = co2_loss2atm          (3) - co2_loss2atm          (4)
         co2_denseffect        (5) = co2_denseffect        (3) - co2_denseffect        (4)
         co2_zcaneffect        (5) = co2_zcaneffect        (3) - co2_zcaneffect        (4)
         co2_gpp               (5) = co2_gpp               (3) - co2_gpp               (4)
         co2_plresp            (5) = co2_plresp            (3) - co2_plresp            (4)
         co2_rh                (5) = co2_rh                (3) - co2_rh                (4)
         carbon_loss2atm       (5) = carbon_loss2atm       (3) - carbon_loss2atm       (4)
         carbon_denseffect     (5) = carbon_denseffect     (3) - carbon_denseffect     (4)
         carbon_zcaneffect     (5) = carbon_zcaneffect     (3) - carbon_zcaneffect     (4)
         carbon_seedrain       (5) = carbon_seedrain       (3) - carbon_seedrain       (4)
         carbon_loss2yield     (5) = carbon_loss2yield     (3) - carbon_loss2yield     (4)
         water_loss2atm        (5) = water_loss2atm        (3) - water_loss2atm        (4)
         water_denseffect      (5) = water_denseffect      (3) - water_denseffect      (4)
         water_wcapeffect      (5) = water_wcapeffect      (3) - water_wcapeffect      (4)
         water_zcaneffect      (5) = water_zcaneffect      (3) - water_zcaneffect      (4)
         water_pheneffect      (5) = water_pheneffect      (3) - water_pheneffect      (4)
         water_precipgain      (5) = water_precipgain      (3) - water_precipgain      (4)
         water_loss2runoff     (5) = water_loss2runoff     (3) - water_loss2runoff     (4)
         water_loss2drainage   (5) = water_loss2drainage   (3) - water_loss2drainage   (4)
         enthalpy_loss2atm     (5) = enthalpy_loss2atm     (3) - enthalpy_loss2atm     (4)
         enthalpy_denseffect   (5) = enthalpy_denseffect   (3) - enthalpy_denseffect   (4)
         enthalpy_prsseffect   (5) = enthalpy_prsseffect   (3) - enthalpy_prsseffect   (4)
         enthalpy_hcapeffect   (5) = enthalpy_hcapeffect   (3) - enthalpy_hcapeffect   (4)
         enthalpy_wcapeffect   (5) = enthalpy_wcapeffect   (3) - enthalpy_wcapeffect   (4)
         enthalpy_zcaneffect   (5) = enthalpy_zcaneffect   (3) - enthalpy_zcaneffect   (4)
         enthalpy_pheneffect   (5) = enthalpy_pheneffect   (3) - enthalpy_pheneffect   (4)
         enthalpy_loss2runoff  (5) = enthalpy_loss2runoff  (3) - enthalpy_loss2runoff  (4)
         enthalpy_loss2drainage(5) = enthalpy_loss2drainage(3) - enthalpy_loss2drainage(4)
         enthalpy_netrad       (5) = enthalpy_netrad       (3) - enthalpy_netrad       (4)
         enthalpy_precipgain   (5) = enthalpy_precipgain   (3) - enthalpy_precipgain   (4)
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Find the relative difference between the two approaches.                    !
         !---------------------------------------------------------------------------------!
         !----- Storage terms. ------------------------------------------------------------!
         if(co2_cas(4) > tiny_num) then
            co2_cas(6) = co2_cas(5) / co2_cas(4)
         else 
            co2_cas(6) = 0. 
         end if
         if(water_soil(4) > tiny_num) then
            water_soil(6) = water_soil(5) / water_soil(4)
         else 
            water_soil(6) = 0. 
         end if
         if(water_sfcwater(4) > tiny_num) then
            water_sfcwater(6) = water_sfcwater(5) / water_sfcwater(4)
         else 
            water_sfcwater(6) = 0. 
         end if
         if(water_cas(4) > tiny_num) then
            water_cas(6) = water_cas(5) / water_cas(4)
         else 
            water_cas(6) = 0. 
         end if
         if(water_leaf(4) > tiny_num) then
            water_leaf(6) = water_leaf(5) / water_leaf(4)
         else 
            water_leaf(6) = 0. 
         end if
         if(water_wood(4) > tiny_num) then
            water_wood(6) = water_wood(5) / water_wood(4)
         else 
            water_wood(6) = 0. 
         end if
         if(enthalpy_soil(4) > tiny_num) then
            enthalpy_soil(6) = enthalpy_soil(5) / enthalpy_soil(4)
         else 
            enthalpy_soil(6) = 0. 
         end if
         if(enthalpy_sfcwater(4) > tiny_num) then
            enthalpy_sfcwater(6) = enthalpy_sfcwater(5) / enthalpy_sfcwater(4)
         else 
            enthalpy_sfcwater(6) = 0. 
         end if
         if(enthalpy_cas(4) > tiny_num) then
            enthalpy_cas(6) = enthalpy_cas(5) / enthalpy_cas(4)
         else 
            enthalpy_cas(6) = 0. 
         end if
         if(enthalpy_leaf(4) > tiny_num) then
            enthalpy_leaf(6) = enthalpy_leaf(5) / enthalpy_leaf(4)
         else 
            enthalpy_leaf(6) = 0. 
         end if
         if(enthalpy_wood(4) > tiny_num) then
            enthalpy_wood(6) = enthalpy_wood(5) / enthalpy_wood(4)
         else 
            enthalpy_wood(6) = 0. 
         end if
         if(carbon_necro(4) > tiny_num) then
            carbon_necro(6) = carbon_necro(5) / carbon_necro(4)
         else 
            carbon_necro(6) = 0. 
         end if
         if(carbon_repro(4) > tiny_num) then
            carbon_repro(6) = carbon_repro(5) / carbon_repro(4)
         else 
            carbon_repro(6) = 0. 
         end if
         if(carbon_balive(4) > tiny_num) then
            carbon_balive(6) = carbon_balive(5) / carbon_balive(4)
         else 
            carbon_balive(6) = 0. 
         end if
         if(carbon_bdead(4) > tiny_num) then
            carbon_bdead(6) = carbon_bdead(5) / carbon_bdead(4)
         else 
            carbon_bdead(6) = 0. 
         end if
         if(carbon_bstorage(4) > tiny_num) then
            carbon_bstorage(6) = carbon_bstorage(5) / carbon_bstorage(4)
         else 
            carbon_bstorage(6) = 0. 
         end if
         if(carbon_cas(4) > tiny_num) then
            carbon_cas(6) = carbon_cas(5) / carbon_cas(4)
         else 
            carbon_cas(6) = 0. 
         end if
         if(carbon_committed(4) > tiny_num) then
            carbon_committed(6) = carbon_committed(5) / carbon_committed(4)
         else 
            carbon_committed(6) = 0. 
         end if
         !------ Budget terms. ------------------------------------------------------------!
         if(abs(co2_loss2atm(4)) > tiny_num) then
            co2_loss2atm(6) = co2_loss2atm(5) / abs(co2_loss2atm(4))
         else
            co2_loss2atm(6) = 0.
         end if
         if(abs(co2_denseffect(4)) > tiny_num) then
            co2_denseffect(6) = co2_denseffect(5) / abs(co2_denseffect(4))
         else
            co2_denseffect(6) = 0.
         end if
         if(abs(co2_zcaneffect(4)) > tiny_num) then
            co2_zcaneffect(6) = co2_zcaneffect(5) / abs(co2_zcaneffect(4))
         else
            co2_zcaneffect(6) = 0.
         end if
         if(abs(co2_gpp(4)) > tiny_num) then
            co2_gpp(6) = co2_gpp(5) / abs(co2_gpp(4))
         else
            co2_gpp(6) = 0.
         end if
         if(abs(co2_plresp(4)) > tiny_num) then
            co2_plresp(6) = co2_plresp(5) / abs(co2_plresp(4))
         else
            co2_plresp(6) = 0.
         end if
         if(abs(co2_rh(4)) > tiny_num) then
            co2_rh(6) = co2_rh(5) / abs(co2_rh(4))
         else
            co2_rh(6) = 0.
         end if
         if(abs(carbon_loss2atm(4)) > tiny_num) then
            carbon_loss2atm(6) = carbon_loss2atm(5) / abs(carbon_loss2atm(4))
         else
            carbon_loss2atm(6) = 0.
         end if
         if(abs(carbon_denseffect(4)) > tiny_num) then
            carbon_denseffect(6) = carbon_denseffect(5) / abs(carbon_denseffect(4))
         else
            carbon_denseffect(6) = 0.
         end if
         if(abs(carbon_zcaneffect(4)) > tiny_num) then
            carbon_zcaneffect(6) = carbon_zcaneffect(5) / abs(carbon_zcaneffect(4))
         else
            carbon_zcaneffect(6) = 0.
         end if
         if(abs(carbon_seedrain(4)) > tiny_num) then
            carbon_seedrain(6) = carbon_seedrain(5) / abs(carbon_seedrain(4))
         else
            carbon_seedrain(6) = 0.
         end if
         if(abs(carbon_loss2yield(4)) > tiny_num) then
            carbon_loss2yield(6) = carbon_loss2yield(5) / abs(carbon_loss2yield(4))
         else
            carbon_loss2yield(6) = 0.
         end if
         if(abs(water_loss2atm(4)) > tiny_num) then
            water_loss2atm(6) = water_loss2atm(5) / abs(water_loss2atm(4))
         else
            water_loss2atm(6) = 0.
         end if
         if(abs(water_denseffect(4)) > tiny_num) then
            water_denseffect(6) = water_denseffect(5) / abs(water_denseffect(4))
         else
            water_denseffect(6) = 0.
         end if
         if(abs(water_wcapeffect(4)) > tiny_num) then
            water_wcapeffect(6) = water_wcapeffect(5) / abs(water_wcapeffect(4))
         else
            water_wcapeffect(6) = 0.
         end if
         if(abs(water_zcaneffect(4)) > tiny_num) then
            water_zcaneffect(6) = water_zcaneffect(5) / abs(water_zcaneffect(4))
         else
            water_zcaneffect(6) = 0.
         end if
         if(abs(water_pheneffect(4)) > tiny_num) then
            water_pheneffect(6) = water_pheneffect(5) / abs(water_pheneffect(4))
         else
            water_pheneffect(6) = 0.
         end if
         if(abs(water_precipgain(4)) > tiny_num) then
            water_precipgain(6) = water_precipgain(5) / abs(water_precipgain(4))
         else
            water_precipgain(6) = 0.
         end if
         if(abs(water_loss2runoff(4)) > tiny_num) then
            water_loss2runoff(6) = water_loss2runoff(5) / abs(water_loss2runoff(4))
         else
            water_loss2runoff(6) = 0.
         end if
         if(abs(water_loss2drainage(4)) > tiny_num) then
            water_loss2drainage(6) = water_loss2drainage(5) / abs(water_loss2drainage(4))
         else
            water_loss2drainage(6) = 0.
         end if
         if(abs(enthalpy_loss2atm(4)) > tiny_num) then
            enthalpy_loss2atm(6) = enthalpy_loss2atm(5) / abs(enthalpy_loss2atm(4))
         else
            enthalpy_loss2atm(6) = 0.
         end if
         if(abs(enthalpy_denseffect(4)) > tiny_num) then
            enthalpy_denseffect(6) = enthalpy_denseffect(5) / abs(enthalpy_denseffect(4))
         else
            enthalpy_denseffect(6) = 0.
         end if
         if(abs(enthalpy_prsseffect(4)) > tiny_num) then
            enthalpy_prsseffect(6) = enthalpy_prsseffect(5) / abs(enthalpy_prsseffect(4))
         else
            enthalpy_prsseffect(6) = 0.
         end if
         if(abs(enthalpy_hcapeffect(4)) > tiny_num) then
            enthalpy_hcapeffect(6) = enthalpy_hcapeffect(5) / abs(enthalpy_hcapeffect(4))
         else
            enthalpy_hcapeffect(6) = 0.
         end if
         if(abs(enthalpy_wcapeffect(4)) > tiny_num) then
            enthalpy_wcapeffect(6) = enthalpy_wcapeffect(5) / abs(enthalpy_wcapeffect(4))
         else
            enthalpy_wcapeffect(6) = 0.
         end if
         if(abs(enthalpy_zcaneffect(4)) > tiny_num) then
            enthalpy_zcaneffect(6) = enthalpy_zcaneffect(5) / abs(enthalpy_zcaneffect(4))
         else
            enthalpy_zcaneffect(6) = 0.
         end if
         if(abs(enthalpy_pheneffect(4)) > tiny_num) then
            enthalpy_pheneffect(6) = enthalpy_pheneffect(5) / abs(enthalpy_pheneffect(4))
         else
            enthalpy_pheneffect(6) = 0.
         end if
         if(abs(enthalpy_loss2runoff(4)) > tiny_num) then
            enthalpy_loss2runoff(6) = enthalpy_loss2runoff(5) / abs(enthalpy_loss2runoff(4))
         else
            enthalpy_loss2runoff(6) = 0.
         end if
         if(abs(enthalpy_loss2drainage(4)) > tiny_num) then
            enthalpy_loss2drainage(6) = enthalpy_loss2drainage(5)                          &
                                      / abs(enthalpy_loss2drainage(4))
         else
            enthalpy_loss2drainage(6) = 0.
         end if
         if(abs(enthalpy_netrad(4)) > tiny_num) then
            enthalpy_netrad(6) = enthalpy_netrad(5) / abs(enthalpy_netrad(4))
         else
            enthalpy_netrad(6) = 0.
         end if
         if(abs(enthalpy_precipgain(4)) > tiny_num) then
            enthalpy_precipgain(6) = enthalpy_precipgain(5) / abs(enthalpy_precipgain(4))
         else
            enthalpy_precipgain(6) = 0.
         end if
         !---------------------------------------------------------------------------------!




         !------ Report the bad news... ---------------------------------------------------!
         write(unit=*,fmt='(a)') '|====================================================|'
         write(unit=*,fmt='(a)') '|====================================================|'
         write(unit=*,fmt='(a)') '|       !!!   Patch Bfusion budget failed   !!!      |'
         write(unit=*,fmt='(a)') '|----------------------------------------------------|'
         write(unit=*,fmt=fmtt ) ' TIME                : ',current_time%year               &
                                                          ,current_time%month              &
                                                          ,current_time%date               &
                                                          ,current_time%time
         write(unit=*,fmt=fmti ) ' IPA                 : ',ipa
         write(unit=*,fmt=fmti ) ' DIST_TYPE(RECP)     : ',osite%dist_type(1)
         write(unit=*,fmt=fmti ) ' DIST_TYPE(DONP)     : ',osite%dist_type(2)
         write(unit=*,fmt=fmtf ) ' AGE(RECP)           : ',osite%age      (1)
         write(unit=*,fmt=fmtf ) ' AGE(DONP)           : ',osite%age      (2)
         write(unit=*,fmt=fmtf ) ' AREA(RECP)          : ',osite%area     (1)
         write(unit=*,fmt=fmtf ) ' AREA(DONP)          : ',osite%area     (2)
         write(unit=*,fmt=fmtf ) ' WEIGHT(RECP)        : ',rawgt
         write(unit=*,fmt=fmtf ) ' WEIGHT(DONP)        : ',dawgt
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' Canopy air space details'
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt=fmtx ) ' Pool           ','      Receptor','         Donor'      &
                                                   ,'      Top-down'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt=fmts ) ' CAN_PRSS     : ', (can_prss    (opa),opa=1,3)
         write(unit=*,fmt=fmts ) ' CAN_TEMP     : ', (can_temp    (opa),opa=1,3)
         write(unit=*,fmt=fmts ) ' CAN_SHV      : ', (can_shv     (opa),opa=1,3)
         write(unit=*,fmt=fmts ) ' CAN_CO2      : ', (can_co2     (opa),opa=1,3)
         write(unit=*,fmt=fmts ) ' CAN_RHOS     : ', (can_rhos    (opa),opa=1,3)
         write(unit=*,fmt=fmts ) ' CAN_DMOL     : ', (can_dmol    (opa),opa=1,3)
         write(unit=*,fmt=fmts ) ' CAN_ENTHALPY : ', (can_enthalpy(opa),opa=1,3)
         write(unit=*,fmt=fmts ) ' CAN_DEPTH    : ', (can_depth   (opa),opa=1,3)
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' Relative tolerances'
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt=fmtf ) ' REL_TOLER_SUBDAY  : ',tol_subday_budget
         write(unit=*,fmt=fmtf ) ' REL_TOLER_CARBON  : ',tol_carbon_budget
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' CO2 Summary'
         write(unit=*,fmt='(a)') ' .................................................. '
         write(unit=*,fmt=fmtl ) ' FINE              : ',co2fused_fine
         write(unit=*,fmt=fmtf ) ' TOLERANCE         : ',co2fused_tolerance
         write(unit=*,fmt=fmtf ) ' RESIDUAL          : ',co2fused_residual
         write(unit=*,fmt=fmtf ) ' STORAGE (FUSED)   : ',co2fused_initialstorage
         write(unit=*,fmt=fmtf ) ' STORAGE (IPA)     : ',csite%co2budget_initialstorage(ipa)
         write(unit=*,fmt=fmtf ) ' DELTA_VEGTDYN     : ',delta_co2_vegtdyn
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' Carbon Summary'
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt=fmtl ) ' FINE              : ',cfused_fine
         write(unit=*,fmt=fmtf ) ' TOLERANCE         : ',cfused_tolerance
         write(unit=*,fmt=fmtf ) ' RESIDUAL          : ',cfused_residual
         write(unit=*,fmt=fmtf ) ' STORAGE (FUSED)   : ',cfused_initialstorage
         write(unit=*,fmt=fmtf ) ' STORAGE (IPA)     : ',csite%cbudget_initialstorage(ipa)
         write(unit=*,fmt=fmtf ) ' DELTA_VEGTDYN     : ',delta_c_vegtdyn
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' Enthalpy Summary'
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt=fmtl ) ' FINE              : ',efused_fine
         write(unit=*,fmt=fmtf ) ' TOLERANCE         : ',efused_tolerance
         write(unit=*,fmt=fmtf ) ' RESIDUAL          : ',efused_residual
         write(unit=*,fmt=fmtf ) ' STORAGE (FUSED)   : ',efused_initialstorage
         write(unit=*,fmt=fmtf ) ' DELTA_VEGTDYN     : ',delta_e_vegtdyn
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' Water Summary'
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt=fmtl ) ' FINE              : ',wfused_fine
         write(unit=*,fmt=fmtf ) ' TOLERANCE         : ',wfused_tolerance
         write(unit=*,fmt=fmtf ) ' RESIDUAL          : ',wfused_residual
         write(unit=*,fmt=fmtf ) ' STORAGE (FUSED)   : ',wfused_initialstorage
         write(unit=*,fmt=fmtf ) ' STORAGE (IPA)     : ',csite%wbudget_initialstorage(ipa)
         write(unit=*,fmt=fmtf ) ' DELTA_VEGTDYN     : ',delta_w_vegtdyn
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' CO2 Detail'
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt=fmth ) ' Pool           ','      Receptor','         Donor'      &
                                                   ,'     Bottom-up','      Top-down'      &
                                                   ,'Abs difference','Rel difference'
         write(unit=*,fmt=fmtd ) ' CANOPY AIR   : ',(co2_cas               (opa),opa=1,6)
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt=fmtd ) ' LOSS2ATM     : ',(co2_loss2atm          (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' DENSEFFECT   : ',(co2_denseffect        (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' ZCANEFFECT   : ',(co2_zcaneffect        (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' GPP          : ',(co2_gpp               (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' PLRESP       : ',(co2_plresp            (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' RH           : ',(co2_rh                (opa),opa=1,6)
         write(unit=*,fmt='(a)') ' .......................................................'


         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' Carbon Detail'
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt=fmth ) ' Pool           ','      Receptor','         Donor'      &
                                                   ,'     Bottom-up','      Top-down'      &
                                                   ,'Abs difference','Rel difference'
         write(unit=*,fmt=fmtd ) ' NECROMASS    : ',(carbon_necro          (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' REPRODUCTION : ',(carbon_repro          (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' BALIVE       : ',(carbon_balive         (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' BDEAD        : ',(carbon_bdead          (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' BSTORAGE     : ',(carbon_bstorage       (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' CANOPY AIR   : ',(carbon_cas            (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' COMMITTED    : ',(carbon_committed      (opa),opa=1,6)
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt=fmtd ) ' LOSS2ATM     : ',(carbon_loss2atm       (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' DENSEFFECT   : ',(carbon_denseffect     (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' ZCANEFFECT   : ',(carbon_zcaneffect     (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' SEEDRAIN     : ',(carbon_seedrain       (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' LOSS2YIELD   : ',(carbon_loss2yield     (opa),opa=1,6)
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' Enthalpy Detail'
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt=fmth ) ' Pool           ','      Receptor','         Donor'      &
                                                   ,'     Bottom-up','      Top-down'      &
                                                   ,'Abs difference','Rel difference'
         write(unit=*,fmt=fmtd ) ' SOIL         : ',(enthalpy_soil         (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' SFC WATER    : ',(enthalpy_sfcwater     (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' CANOPY AIR   : ',(enthalpy_cas          (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' LEAF         : ',(enthalpy_leaf         (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' WOOD         : ',(enthalpy_wood         (opa),opa=1,6)
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt=fmtd ) ' LOSS2ATM     : ',(enthalpy_loss2atm     (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' DENSEFFECT   : ',(enthalpy_denseffect   (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' PRSSEFFECT   : ',(enthalpy_prsseffect   (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' HCAPEFFECT   : ',(enthalpy_hcapeffect   (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' WCAPEFFECT   : ',(enthalpy_wcapeffect   (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' ZCANEFFECT   : ',(enthalpy_zcaneffect   (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' PHENEFFECT   : ',(enthalpy_pheneffect   (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' LOSS2RUNOFF  : ',(enthalpy_loss2runoff  (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' LOSS2DRAINAGE: ',(enthalpy_loss2drainage(opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' NETRAD       : ',(enthalpy_netrad       (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' PRECIPGAIN   : ',(enthalpy_precipgain   (opa),opa=1,6)
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' Water Detail'
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt=fmth ) ' Pool           ','      Receptor','         Donor'      &
                                                   ,'     Bottom-up','      Top-down'      &
                                                   ,'Abs difference','Rel difference'
         write(unit=*,fmt=fmtd ) ' SOIL         : ',(water_soil            (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' SFC WATER    : ',(water_sfcwater        (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' CANOPY AIR   : ',(water_cas             (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' LEAF         : ',(water_leaf            (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' WOOD         : ',(water_wood            (opa),opa=1,6)
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt=fmtd ) ' LOSS2ATM     : ',(water_loss2atm        (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' DENSEFFECT   : ',(water_denseffect      (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' WCAPEFFECT   : ',(water_wcapeffect      (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' ZCANEFFECT   : ',(water_zcaneffect      (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' PHENEFFECT   : ',(water_pheneffect      (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' LOSS2RUNOFF  : ',(water_precipgain      (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' LOSS2DRAINAGE: ',(water_loss2runoff     (opa),opa=1,6)
         write(unit=*,fmt=fmtd ) ' PRECIPGAIN   : ',(water_loss2drainage   (opa),opa=1,6)
         write(unit=*,fmt='(a)') ' .......................................................'
         write(unit=*,fmt='(a)') ' '
         write(unit=*,fmt='(a)') ' -------------------------------------------------------'
         !---------------------------------------------------------------------------------!

         call fatal_error('Budget check has failed, see message above.'                    &
                         ,'check_bfusion_patch','fuse_fiss_utils.f90')
      end if
      !------------------------------------------------------------------------------------!

      return
   end subroutine check_bfusion_patch
   !=======================================================================================!
   !=======================================================================================!

end module fuse_fiss_utils
!==========================================================================================!
!==========================================================================================!






