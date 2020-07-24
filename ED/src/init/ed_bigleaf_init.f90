!==========================================================================================!
!==========================================================================================!
!     This sub-routine converts a conventional SAS initialisation into a big leaf.         !
!------------------------------------------------------------------------------------------!
subroutine ed_bigleaf_init(cgrid)
   use ed_node_coms   , only : mynum                ! ! intent(in)
   use ed_state_vars  , only : edtype               & ! structure
                             , polygontype          & ! structure
                             , sitetype             & ! structure
                             , patchtype            & ! structure
                             , allocate_sitetype    & ! subroutine
                             , deallocate_sitetype  & ! subroutine
                             , allocate_patchtype   & ! subroutine
                             , deallocate_patchtype ! ! subroutine
   use ed_max_dims    , only : n_pft                & ! intent(in)
                             , n_dist_types         ! ! intent(in)
   use allometry      , only : size2bl              & ! function
                             , dbh2bd               & ! function
                             , area_indices         & ! function
                             , ed_biomass           ! ! function
   use pft_coms       , only : hgt_max              & ! intent(in)
                             , dbh_bigleaf          & ! intent(in)
                             , sla                  & ! intent(in)
                             , q                    & ! intent(in)
                             , qsw                  & ! intent(in)
                             , agf_bs               ! ! intent(in)
   use fuse_fiss_utils, only : sort_cohorts         & ! subroutine
                             , sort_patches         ! ! subroutine
   use consts_coms    , only : pio4                 ! ! intent(in)
   use grid_coms      , only : nzg                  ! ! intent(in)
   implicit none

   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target                        :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer                       :: cpoly
   type(sitetype)   , pointer                       :: csite
   type(patchtype)  , pointer                       :: cpatch
   integer          , dimension(n_dist_types)       :: lu_npatch
   integer                                          :: site_npatch
   integer                                          :: ipy
   integer                                          :: isi
   integer                                          :: ipa
   integer                                          :: ico
   integer                                          :: ipft
   integer                                          :: ilu
   integer                                          :: ncohorts
   integer                                          :: npatchco
   integer                                          :: nsitepat
   integer                                          :: k
   logical          , dimension(n_dist_types)       :: lu_desert
   real             , dimension(n_pft)              :: pft_area
   real             , dimension(n_pft,n_dist_types) :: lai
   real             , dimension(n_dist_types)       :: area
   real             , dimension(n_dist_types)       :: fsc
   real             , dimension(n_dist_types)       :: ssc
   real             , dimension(n_dist_types)       :: stsc
   real             , dimension(n_dist_types)       :: stsl
   real             , dimension(n_dist_types)       :: msn
   real             , dimension(n_dist_types)       :: fsn
   real             , dimension(n_dist_types)       :: sum_dgd
   real             , dimension(n_dist_types)       :: sum_chd
   real                                             :: area_sum
   real                                             :: area_tot
   real                                             :: patch_lai
   real                                             :: patch_wai
   real                                             :: patch_agb
   real                                             :: patch_bsa
   real                                             :: patch_pop
   real                                             :: site_lai
   real                                             :: site_agb
   real                                             :: site_bsa
   real                                             :: site_pop
   real                                             :: site_fsc
   real                                             :: site_ssc
   real                                             :: site_stsc
   real                                             :: poly_lai
   real                                             :: poly_agb
   real                                             :: poly_bsa
   real                                             :: poly_pop
   real                                             :: poly_fsc
   real                                             :: poly_ssc
   real                                             :: poly_stsc
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !       Loop over all sites.                                                            !
   !---------------------------------------------------------------------------------------!
   write(unit=*,fmt='(141a)'    ) ('-',k=1,141)
   write(unit=*,fmt='(a,1x,i5)' ) ' Mynum = ',mynum
   write(unit=*,fmt='(141a)'    ) ('-',k=1,141)
   write(unit=*,fmt='(11(a,1x))') '         IPY','      NSITES','    NPATCHES'             &
                                 ,'    NCOHORTS','      NPLANT','         LAI'             &
                                 ,'         AGB','  BASAL_AREA','     FAST_SC'             &
                                 ,'     SLOW_SC','   STRUCT_SC'
   write(unit=*,fmt='(141a)'   ) ('-',k=1,141)
   polyloop: do ipy=1,cgrid%npolygons

      cpoly => cgrid%polygon(ipy)
      firstsiteloop: do isi=1,cpoly%nsites

         !----- Reset all structures before loop over paches. -----------------------------!
         lai     (:,:) = 0.
         area      (:) = 0.
         fsc       (:) = 0.
         ssc       (:) = 0.
         stsc      (:) = 0.
         stsl      (:) = 0.
         msn       (:) = 0.
         fsn       (:) = 0.
         sum_dgd   (:) = 0.
         sum_chd   (:) = 0.

         !---------------------------------------------------------------------------------!
         !     Determine the site-level properties.                                        !
         !---------------------------------------------------------------------------------!
         csite => cpoly%site(isi)
         inpatchloop: do ipa = 1,csite%npatches

            !------------------------------------------------------------------------------!
            !    Find the land use type for this patch.  We will aggregate all patches     !
            ! with this disturbance type.                                                  !
            !------------------------------------------------------------------------------!
            ilu = csite%dist_type(ipa)
            area   (ilu) = area   (ilu) + csite%area              (ipa)
            fsc  (ilu) = fsc  (ilu) + csite%fast_soil_C       (nzg,ipa) * csite%area(ipa)
            ssc  (ilu) = ssc  (ilu) + csite%slow_soil_C       (nzg,ipa) * csite%area(ipa)
            stsc (ilu) = stsc (ilu) + csite%structural_soil_C (nzg,ipa) * csite%area(ipa)
            stsl (ilu) = stsl (ilu) + csite%structural_soil_L (nzg,ipa) * csite%area(ipa)
            msn  (ilu) = msn  (ilu) + csite%mineralized_soil_N(nzg,ipa) * csite%area(ipa)
            fsn  (ilu) = fsn  (ilu) + csite%fast_soil_N       (nzg,ipa) * csite%area(ipa)
            sum_dgd(ilu) = sum_dgd(ilu) + csite%sum_dgd           (ipa) * csite%area(ipa)
            sum_chd(ilu) = sum_chd(ilu) + csite%sum_chd           (ipa) * csite%area(ipa)
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !    Integrate the cohort-level properties.                                    !
            !------------------------------------------------------------------------------!
            cpatch => csite%patch(ipa)
            incohortloop: do ico=1,cpatch%ncohorts
               !---------------------------------------------------------------------------!
               !     Find the PFT for this cohort.  We aggregate all cohorts of this PFT.  !
               ! The only variable we need is the leaf area index, which is preserved at   !
               ! the site level.  Because this subroutine is called immediately after the  !
               ! initialisation and before the first phenology check, LAI is always        !
               ! positive as long as the cohort exists.                                    !
               !---------------------------------------------------------------------------!
               ipft = cpatch%pft(ico)
               lai(ipft,ilu) = lai(ipft,ilu) + cpatch%lai(ico) * csite%area(ipa)
               !---------------------------------------------------------------------------!
            end do incohortloop
            !------------------------------------------------------------------------------!


            !----- Remove all cohorts. ----------------------------------------------------!
            call deallocate_patchtype(cpatch)
            cpatch%ncohorts = 0
            !------------------------------------------------------------------------------!
         end do inpatchloop
         !---------------------------------------------------------------------------------!


         !----- Remove all patches. -------------------------------------------------------!
         call deallocate_sitetype(csite)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     The new number of patches is going to be the total number of PFTs that      !
         ! exist amongst all land use types.                                               !
         ! 1.  We add bdead and bleaf because the new grasses don't have structural        !
         !     biomass.                                                                    !
         ! 2.  In case the biomass is zero but the area of the type is not, we make one    !
         !     empty patch.  This can happen when the previous simulation evolved to a     !
         !     desert, for example.                                                        !
         !---------------------------------------------------------------------------------!
         do ilu=1,n_dist_types
            lu_npatch(ilu) = count(lai(:,ilu) > 0.)
            lu_desert(ilu) = lu_npatch(ilu) == 0 .and. area(ilu) > 0.
            if (lu_desert(ilu)) lu_npatch(ilu) = 1
         end do
         site_npatch = sum(lu_npatch)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Allocate the right number of patches.                                      !
         !---------------------------------------------------------------------------------!
         call allocate_sitetype(csite,site_npatch)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Populate th patches with either 0 or 1 cohort each.                        !
         !---------------------------------------------------------------------------------!
         ipa = 0
         addluloop: do ilu=1,n_dist_types
            if (lu_desert(ilu)) then
               !---------------------------------------------------------------------------!
               !    This land use type existed but became a desert.  Create an empty       !
               ! patch. EJL - adding all C and N to thickest layer for now                 !
               !---------------------------------------------------------------------------!
               ipa = ipa + 1
               csite%dist_type         (ipa) = ilu
               csite%area              (ipa) = area(ilu)
               csite%age               (ipa) = 0.0
               csite%fast_soil_C       (1,ipa) = fsc (ilu) / area(ilu)
               csite%slow_soil_C       (1,ipa) = ssc (ilu) / area(ilu)
               csite%structural_soil_C (1,ipa) = stsc(ilu) / area(ilu)
               csite%structural_soil_L (1,ipa) = stsl(ilu) / area(ilu)
               csite%mineralized_soil_N(1,ipa) = msn (ilu) / area(ilu)
               csite%fast_soil_N       (1,ipa) = fsn (ilu) / area(ilu)
               csite%sum_dgd           (ipa) = 0.0
               csite%sum_chd           (ipa) = 0.0
               csite%cohort_count      (ipa) = 0
               !---------------------------------------------------------------------------!

            elseif (lu_npatch(ilu) /= 0) then
               !---------------------------------------------------------------------------!
               !    This land use type has cohorts.   The first step is to find the total  !
               ! area that each patch/cohort/PFT is entitled.                              !
               !---------------------------------------------------------------------------!
               pft_area(:) = lai(:,ilu) / sum(lai(:,ilu))

               !----- Loop over all PFTs, and add those that have something. --------------!
               addpftloop: do ipft = 1,n_pft
                  if (lai(ipft,ilu) > 0.0) then
                     ipa = ipa + 1
                     !----- Initialise the patch-level properties. ------------------------!
                     ! EJL - Adding C and N to first layer of pool for now
                     csite%dist_type         (ipa) = ilu
                     csite%area              (ipa) = pft_area(ipft) * area(ilu)
                     csite%age               (ipa) = 0.0
                     csite%fast_soil_C       (1,ipa) = fsc (ilu) / area(ilu)
                     csite%slow_soil_C       (1,ipa) = ssc (ilu) / area(ilu)
                     csite%structural_soil_C (1,ipa) = stsc(ilu) / area(ilu)
                     csite%structural_soil_L (1,ipa) = stsl(ilu) / area(ilu)
                     csite%mineralized_soil_N(1,ipa) = msn (ilu) / area(ilu)
                     csite%fast_soil_N       (1,ipa) = fsn (ilu) / area(ilu)
                     csite%sum_dgd           (ipa) = sum_dgd(ilu)
                     csite%sum_chd           (ipa) = sum_chd(ilu)
                     csite%cohort_count      (ipa) = 1

                     !----- Big-leaf model, therefore only one cohort per patch. ----------!
                     cpatch => csite%patch(ipa)
                     call allocate_patchtype(cpatch,1)
                     cpatch%pft          (1) = ipft
                     cpatch%hite         (1) = hgt_max(ipft)
                     cpatch%dbh          (1) = dbh_bigleaf(ipft)
                     cpatch%bleaf        (1) = size2bl(cpatch%dbh(1),cpatch%hite(1),ipft)
                     cpatch%bdead        (1) = dbh2bd (cpatch%dbh(1),ipft)
                     cpatch%broot        (1) = cpatch%bleaf(1) * q(ipft)
                     cpatch%bsapwooda    (1) = agf_bs(ipft) * cpatch%bleaf(1)              &
                                             * qsw(ipft) * cpatch%hite(1)
                     cpatch%bsapwoodb    (1) = (1.0 - agf_bs(ipft)) * cpatch%bleaf(1)      &
                                             * qsw(ipft) * cpatch%hite(1)
                     cpatch%balive       (1) = cpatch%bleaf    (1) + cpatch%broot    (1)   &
                                             + cpatch%bsapwooda(1) + cpatch%bsapwoodb(1)
                     cpatch%sla          (1) = sla     (ipft)
                     cpatch%nplant       (1) = lai (ipft,ilu)                              &
                                             / ( cpatch%sla(1) * cpatch%bleaf(1)           &
                                               * csite%area(ipa) )
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !     Start plants with full phenology and in great carbon balance,   !
                     ! we will take care of phenology after this sub-routine.              !
                     !---------------------------------------------------------------------!
                     cpatch%cb               (1:12,1) = 1.0
                     cpatch%cb_lightmax      (1:12,1) = 1.0
                     cpatch%cb_moistmax      (1:12,1) = 1.0
                     cpatch%cb_mlmax         (1:12,1) = 1.0
                     cpatch%cb                 (13,1) = 0.0
                     cpatch%cb_lightmax        (13,1) = 0.0
                     cpatch%cb_moistmax        (13,1) = 0.0
                     cpatch%cb_mlmax           (13,1) = 0.0
                     cpatch%phenology_status      (1) = 0
                     cpatch%bstorage              (1) = 0.0
                     !---------------------------------------------------------------------!



                     !----- Assign LAI, WAI, and CAI --------------------------------------!
                     call area_indices(cpatch, 1)

                     !----- Above ground biomass, use the allometry. ----------------------!
                     cpatch%agb(1)     = ed_biomass(cpatch, 1)
                     cpatch%basarea(1) = pio4 * cpatch%dbh(1) * cpatch%dbh(1)

                     !----- Growth rates, start with zero. --------------------------------!
                     cpatch%dagb_dt  (1)  = 0.
                     cpatch%dlnagb_dt(1)  = 0.
                     cpatch%dba_dt   (1)  = 0.
                     cpatch%dlnba_dt (1)  = 0.
                     cpatch%ddbh_dt  (1)  = 0.
                     cpatch%dlndbh_dt(1)  = 0.

                     !---------------------------------------------------------------------!
                     !      Initialise other cohort variables.  Some of them won't be      !
                     ! updated unless the lai exceeds lai_min.                             !
                     !---------------------------------------------------------------------!
                     cpatch%fsw  (1) = 1.0
                     cpatch%gpp  (1) = 0.0
                     cpatch%par_l(1) = 0.0
                     !---------------------------------------------------------------------!

                     !----- Update the patch level above-ground biomass. ------------------!
                     csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)             &
                                                 + cpatch%agb(1) * cpatch%nplant(1)
                     !---------------------------------------------------------------------!
                  end if
               end do addpftloop
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do addluloop
         !---------------------------------------------------------------------------------!
      end do firstsiteloop
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Initialise all the other site , patch , and cohort level variables.           !
      !------------------------------------------------------------------------------------!
      initsiteloop: do isi = 1,cpoly%nsites
         area_sum = 0.0
         ncohorts = 0


         !----- Make sure that the total patch area is 1. ---------------------------------!
         csite => cpoly%site(isi)
         area_tot      = sum(csite%area(1:csite%npatches))
         csite%area(:) = csite%area(:)/area_tot

         !----- Find the patch-level LAI, WAI, and CAI. -----------------------------------!
         do ipa=1,csite%npatches
            area_sum        = area_sum + csite%area(ipa)
            patch_lai       = 0.0
            patch_wai       = 0.0
            cpatch => csite%patch(ipa)
            do ico = 1,cpatch%ncohorts
               patch_lai = patch_lai + cpatch%lai(ico)
               patch_wai = patch_wai + cpatch%wai(ico)
               ncohorts  = ncohorts + 1
            end do
         end do

         !----- Initialise the cohort variables, then sort them by size. ------------------!
         do ipa = 1,csite%npatches
            cpatch => csite%patch(ipa)
            do ico = 1,cpatch%ncohorts
               call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi))
            end do

            !----- Make sure that cohorts are organised from tallest to shortest. ---------!
            call sort_cohorts(cpatch)
         end do

         !----- Initialise the patch-level variables. -------------------------------------!
         call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))

         !----- Make sure that patches are organised from oldest to youngest. -------------!
         call sort_patches(csite)
      end do initsiteloop
      !------------------------------------------------------------------------------------!



      !----- Initialise site-level variables. ---------------------------------------------!
      call init_ed_site_vars(cpoly)


      !----- Get a diagnostic on the polygon's vegetation. --------------------------------!
      poly_pop  = 0.0
      poly_lai  = 0.0
      poly_agb  = 0.0
      poly_bsa  = 0.0
      poly_fsc  = 0.0
      poly_ssc  = 0.0
      poly_stsc = 0.0
      ncohorts  = 0

      do isi = 1,cpoly%nsites
         nsitepat = 0
         csite => cpoly%site(isi)

         site_pop  = 0.0
         site_lai  = 0.0
         site_agb  = 0.0
         site_bsa  = 0.0
         site_fsc  = 0.0
         site_ssc  = 0.0
         site_stsc = 0.0
         do ipa = 1,csite%npatches
            patch_lai  = 0.0
            patch_wai  = 0.0
            npatchco        = 0

            cpatch => csite%patch(ipa)
            patch_pop       = 0.0
            patch_agb       = 0.0
            patch_bsa       = 0.0
            do ico = 1,cpatch%ncohorts
               ncohorts   = ncohorts+1
               npatchco   = npatchco+1
               patch_lai  = patch_lai  + cpatch%lai(ico)
               patch_wai  = patch_wai  + cpatch%wai(ico)
               patch_pop  = patch_pop  + cpatch%nplant(ico)
               patch_agb  = patch_agb  + cpatch%nplant(ico) * cpatch%agb    (ico)
               patch_bsa  = patch_bsa  + cpatch%nplant(ico) * cpatch%basarea(ico)
            end do
            site_pop  = site_pop  + patch_pop                    * csite%area(ipa)
            site_lai  = site_lai  + patch_lai                    * csite%area(ipa)
            site_agb  = site_agb  + patch_agb                    * csite%area(ipa)
            site_bsa  = site_bsa  + patch_bsa                    * csite%area(ipa)
            site_fsc  = site_fsc  + csite%fast_soil_C      (nzg,ipa) * csite%area(ipa)
            site_ssc  = site_ssc  + csite%slow_soil_C      (nzg,ipa) * csite%area(ipa)
            site_stsc = site_stsc + csite%structural_soil_C(nzg,ipa) * csite%area(ipa)
            csite%cohort_count(ipa) = npatchco
            nsitepat                = nsitepat + 1
         end do

         poly_pop  = poly_pop  + site_pop  * cpoly%area(isi)
         poly_lai  = poly_lai  + site_lai  * cpoly%area(isi)
         poly_agb  = poly_agb  + site_agb  * cpoly%area(isi)
         poly_bsa  = poly_bsa  + site_bsa  * cpoly%area(isi)
         poly_fsc  = poly_fsc  + site_fsc  * cpoly%area(isi)
         poly_ssc  = poly_ssc  + site_ssc  * cpoly%area(isi)
         poly_stsc = poly_stsc + site_stsc * cpoly%area(isi)

         cpoly%patch_count(isi) = nsitepat


         write (unit=*,fmt='(4(i12,1x),7(f12.3,1x))') ipy,cpoly%nsites,nsitepat,ncohorts   &
                                                     ,poly_pop,poly_lai,poly_agb,poly_bsa  &
                                                     ,poly_fsc,poly_ssc,poly_stsc
      end do

      !----- Initialise the polygon-level variables. --------------------------------------!
      call init_ed_poly_vars(cgrid)
   end do polyloop
   write(unit=*,fmt='(141a)') ('-',k=1,141)
   write(unit=*,fmt='(a)') ' '
   write(unit=*,fmt='(a)') ' '
   !---------------------------------------------------------------------------------------!
   return
end subroutine ed_bigleaf_init
!==========================================================================================!
!==========================================================================================!
