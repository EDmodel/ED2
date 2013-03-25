!==========================================================================================!
!==========================================================================================!
!     This subroutine will drive the update of derived properties.                         !
!------------------------------------------------------------------------------------------!
subroutine update_derived_props(cgrid)
   use ed_state_vars , only : edtype      & ! structure
                            , polygontype & ! structure
                            , sitetype    ! ! structure
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)      , target  :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   type(sitetype)    , pointer :: csite
   integer                     :: ipy
   integer                     :: isi
   integer                     :: ipa
   !---------------------------------------------------------------------------------------!
   
   do ipy = 1,cgrid%npolygons
     cpoly => cgrid%polygon(ipy)
     
     do isi = 1,cpoly%nsites
        csite => cpoly%site(isi)

        do ipa = 1,csite%npatches
           call update_patch_derived_props(csite,cpoly%lsl(isi),cpoly%met(isi)%prss,ipa)
        end do

        call update_site_derived_props(cpoly, 0, isi)
     end do

     call update_polygon_derived_props(cgrid)
   end do

   return
end subroutine update_derived_props
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine will take care of derived patch-level structural quantities.  These !
! depend on the results from reproduction, which in turn depends on structural growth      !
! results from all patches.                                                                !
!------------------------------------------------------------------------------------------!
subroutine update_patch_derived_props(csite,lsl,prss,ipa)
  
   use ed_state_vars       , only : sitetype                   & ! structure
                                  , patchtype                  ! ! structure
   use allometry           , only : ed_biomass                 ! ! function
   use fuse_fiss_utils     , only : patch_pft_size_profile     ! ! subroutine
   use canopy_air_coms     , only : veg_height_min             & ! intent(in)
                                  , minimum_canopy_depth       & ! intent(in)
                                  , ez                         & ! intent(in)
                                  , vh2vr                      & ! intent(in)
                                  , vh2dh                      ! ! intent(in)
   use soil_coms           , only : soil_rough                 & ! intent(in)
                                  , ny07_eq04_a                & ! intent(in)
                                  , ny07_eq04_m                & ! intent(in)
                                  , tiny_sfcwater_mass         ! ! intent(in)
   use consts_coms         , only : wdns                       & ! intent(in)
                                  , fsdns                      & ! intent(in)
                                  , fsdnsi                     ! ! intent(in)
   
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)  , target     :: csite
   integer         , intent(in) :: ipa
   integer         , intent(in) :: lsl
   real            , intent(in) :: prss
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype) , pointer    :: cpatch
   real                         :: weight
   real                         :: weight_sum
   real                         :: total_sfcw_mass
   real                         :: bulk_sfcw_dens
   integer                      :: ico
   integer                      :: k
   integer                      :: ksn
   integer                      :: ipft
   !---------------------------------------------------------------------------------------!


   !----- Find the total snow depth. ------------------------------------------------------!
   ksn = csite%nlev_sfcwater(ipa)
   csite%total_sfcw_depth(ipa) = 0.
   total_sfcw_mass             = 0.
   do k=1,ksn
      csite%total_sfcw_depth(ipa) = csite%total_sfcw_depth(ipa)                            &
                                  + csite%sfcwater_depth(k,ipa)
      total_sfcw_mass             = total_sfcw_mass                                        &
                                  + csite%sfcwater_mass(k,ipa)
   end do
   !---------------------------------------------------------------------------------------!




   !----- Reset properties. ---------------------------------------------------------------!
   csite%veg_height(ipa)       = 0.0
   weight_sum                  = 0.0
   csite%opencan_frac(ipa)     = 1.0
   csite%plant_ag_biomass(ipa) = 0.0
   !---------------------------------------------------------------------------------------!


   cpatch => csite%patch(ipa)

   !----- Loop over cohorts and integrate the patch-level properties. ---------------------!
   do ico = 1,cpatch%ncohorts

      ipft = cpatch%pft(ico)


      !----- Compute the patch-level above-ground biomass
      csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)                            &
                                  + ed_biomass(cpatch%bdead(ico),cpatch%bleaf(ico)         &
                                              ,cpatch%bsapwooda(ico),cpatch%pft(ico))                      &
                                  * cpatch%nplant(ico)           
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Compute average vegetation height, weighting using basal area.  We add the     !
      ! cohorts only until when the canopy is closed, this way we will not bias the        !
      ! vegetation height or the canopy depth towards the cohorts that live in the under-  !
      ! storey.  Also, we must take into account the depth of the temporary surface water  !
      ! or snow, because this will make the plants "shorter".                              !
      !------------------------------------------------------------------------------------!
      if (csite%opencan_frac(ipa) > 0.0) then
         weight                  = cpatch%nplant(ico) * cpatch%basarea(ico)
         weight_sum              = weight_sum + weight
         csite%veg_height(ipa)   = csite%veg_height(ipa) + cpatch%hite(ico) * weight
         csite%opencan_frac(ipa) = csite%opencan_frac(ipa) * (1.0 - cpatch%crown_area(ico))
      end if
      !------------------------------------------------------------------------------------!

   end do
   !---------------------------------------------------------------------------------------!



   !----- Normalise the vegetation height, making sure that it is above the minimum. ------!
   if (weight_sum > tiny(1.0)) then
      csite%veg_height(ipa)  = max(veg_height_min,csite%veg_height(ipa) / weight_sum)
   else
      csite%veg_height(ipa)  = veg_height_min
   end if
   !---------------------------------------------------------------------------------------!



   !----- Find the patch roughness due to vegetation. -------------------------------------!
   csite%veg_rough(ipa) = vh2vr * csite%veg_height(ipa)
   !---------------------------------------------------------------------------------------!



   !----- Find the 0-plane displacement due to vegetation. --------------------------------!
   csite%veg_displace(ipa) = vh2dh * csite%veg_height(ipa)
   !---------------------------------------------------------------------------------------!



   !----- Update the canopy depth, and impose the minimum if needed be. -------------------!
   csite%can_depth(ipa) = max(csite%veg_height(ipa), minimum_canopy_depth)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the fraction of the canopy covered in snow.  I could not find any            !
   ! reference for the original method (commented out), so I implemented the method used   !
   ! in CLM-4, which is based on:                                                          !
   !                                                                                       !
   ! Niu, G.-Y., and Z.-L. Yang (2007), An observation-based formulation of snow cover     !
   !    fraction and its evaluation over large North American river basins,                !
   !    J. Geophys. Res., 112, D21101, doi:10.1029/2007JD008674                            !
   !---------------------------------------------------------------------------------------!
   ! csite%snowfac(ipa) = min(0.99, csite%total_sfcw_depth(ipa)/csite%veg_height(ipa))
   if (total_sfcw_mass > tiny_sfcwater_mass) then
      bulk_sfcw_dens     = max( fsdns, min( wdns                                           &
                              , total_sfcw_mass / csite%total_sfcw_depth(ipa)))
      csite%snowfac(ipa) = max( 0.0, min( 0.99                                             &
                              , tanh( csite%total_sfcw_depth(ipa)                          &
                                    / ( ny07_eq04_a * soil_rough                           &
                                      * (bulk_sfcw_dens * fsdnsi) ** ny07_eq04_m ) ) ) )
   else
      csite%snowfac(ipa) = 0.0
   end if
   !---------------------------------------------------------------------------------------!



   !----- Find the PFT-dependent size distribution of this patch. -------------------------!
   call patch_pft_size_profile(csite,ipa)
   !---------------------------------------------------------------------------------------!


   !----- Update the cohort count (may be redundant as well...) ---------------------------!
   csite%cohort_count(ipa) = cpatch%ncohorts
   !---------------------------------------------------------------------------------------!

   return
end subroutine update_patch_derived_props
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine will take care of some diagnostic thermodynamic properties.         !
!------------------------------------------------------------------------------------------!
subroutine update_patch_thermo_props(csite,ipaa,ipaz,mzg,mzs,ntext_soil)
  
   use ed_state_vars, only : sitetype         ! ! structure
   use therm_lib    , only : idealdenssh      & ! function
                           , press2exner      & ! function
                           , extheta2temp     & ! function
                           , uextcm2tl        & ! function
                           , uint2tl          ! ! function
   use consts_coms  , only : t00              & ! intent(in)
                           , wdns             ! ! intent(in)
   use soil_coms    , only : soil             & ! intent(in)
                           , matric_potential ! ! function
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)                , target     :: csite
   integer                       , intent(in) :: ipaa
   integer                       , intent(in) :: ipaz
   integer                       , intent(in) :: mzg
   integer                       , intent(in) :: mzs
   integer       , dimension(mzg), intent(in) :: ntext_soil
   !----- Local variables. ----------------------------------------------------------------!
   integer                                    :: ipa
   integer                                    :: nsoil
   integer                                    :: ksn
   integer                                    :: k
   real                                       :: soilhcap
   real                                       :: can_exner
   !---------------------------------------------------------------------------------------!


   do ipa=ipaa,ipaz

      !----- Canopy air temperature and density. ------------------------------------------!
      can_exner           = press2exner (csite%can_prss(ipa))
      csite%can_temp(ipa) = extheta2temp(can_exner,csite%can_theta(ipa))
      csite%can_rhos(ipa) = idealdenssh ( csite%can_prss  (ipa)                            &
                                        , csite%can_temp  (ipa)                            &
                                        , csite%can_shv   (ipa)                            )
      !------------------------------------------------------------------------------------!


      !----- Update soil temperature and liquid water fraction. ---------------------------!
      do k = 1, mzg
         nsoil    = ntext_soil(k)
         soilhcap = soil(nsoil)%slcpd
         call uextcm2tl(csite%soil_energy(k,ipa),csite%soil_water(k,ipa)*wdns,soilhcap     &
                       ,csite%soil_tempk(k,ipa),csite%soil_fracliq(k,ipa))
         csite%soil_mstpot(k,ipa) = matric_potential(nsoil,csite%soil_water(k,ipa))
      end do
      !------------------------------------------------------------------------------------!



      !----- Update temporary surface water temperature and liquid water fraction. --------!
      ksn = csite%nlev_sfcwater(ipa)
      csite%total_sfcw_depth(ipa) = 0.
      do k = 1, ksn
         call uint2tl(csite%sfcwater_energy(k,ipa),csite%sfcwater_tempk(k,ipa)             &
                     ,csite%sfcwater_fracliq(k,ipa))
         csite%total_sfcw_depth(ipa) =  csite%total_sfcw_depth(ipa)                        &
                                     +  csite%sfcwater_depth(k,ipa)
      end do
      do k = ksn+1,mzs
         if (k == 1) then
            csite%sfcwater_tempk  (k,ipa) = csite%soil_tempk  (mzg,ipa)
            csite%sfcwater_fracliq(k,ipa) = csite%soil_fracliq(mzg,ipa)
         else
            csite%sfcwater_tempk  (k,ipa) = csite%sfcwater_tempk  (k-1,ipa)
            csite%sfcwater_fracliq(k,ipa) = csite%sfcwater_fracliq(k-1,ipa)
         end if
      end do
      !------------------------------------------------------------------------------------!
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine update_patch_thermo_props
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine will update the fast mean properties, similarly to the routine      !
! above.                                                                                   !
!------------------------------------------------------------------------------------------!
subroutine update_patch_thermo_fmean(csite,ipaa,ipaz,mzg,ntext_soil)
  
   use ed_state_vars, only : sitetype           ! ! structure
   use therm_lib    , only : idealdenssh        & ! function
                           , press2exner        & ! function
                           , extheta2temp       & ! function
                           , uextcm2tl          & ! function
                           , uint2tl            ! ! function
   use consts_coms  , only : t00                & ! intent(in)
                           , wdns               ! ! intent(in)
   use soil_coms    , only : soil               & ! intent(in)
                           , tiny_sfcwater_mass & ! intent(in)
                           , matric_potential   ! ! function
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)                , target     :: csite
   integer                       , intent(in) :: ipaa
   integer                       , intent(in) :: ipaz
   integer                       , intent(in) :: mzg
   integer       , dimension(mzg), intent(in) :: ntext_soil
   !----- Local variables. ----------------------------------------------------------------!
   integer                                    :: ipa
   integer                                    :: nsoil
   integer                                    :: ksn
   integer                                    :: k
   real                                       :: soilhcap
   real                                       :: can_exner
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   do ipa=ipaa,ipaz

      !----- Canopy air temperature and density. ------------------------------------------!
      can_exner                 = press2exner (csite%fmean_can_prss(ipa))
      csite%fmean_can_temp(ipa) = extheta2temp(can_exner,csite%fmean_can_theta(ipa))
      csite%fmean_can_rhos(ipa) = idealdenssh ( csite%fmean_can_prss  (ipa)                &
                                              , csite%fmean_can_temp  (ipa)                &
                                              , csite%fmean_can_shv   (ipa)                )
      !------------------------------------------------------------------------------------!


      !----- Update soil temperature and liquid water fraction. ---------------------------!
      do k = 1, mzg
         nsoil    = ntext_soil(k)
         soilhcap = soil(nsoil)%slcpd
         call uextcm2tl( csite%fmean_soil_energy(k,ipa)                                    &
                       , csite%fmean_soil_water (k,ipa) * wdns                             &
                       , soilhcap                                                          &
                       , csite%fmean_soil_temp  (k,ipa)                                    &
                       , csite%fmean_soil_fliq  (k,ipa) )
         csite%fmean_soil_mstpot(k,ipa) = matric_potential( nsoil                          &
                                                          , csite%fmean_soil_water(k,ipa))
      end do
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !   If the patch had some temporary snow/pounding layer, convert the mean energy to  !
      ! J/kg, then find the mean temperature and liquid fraction.  Otherwise, set them to  !
      ! either zero or default values.                                                     !
      !------------------------------------------------------------------------------------!
      if (csite%fmean_sfcw_mass(ipa) > tiny_sfcwater_mass) then
         csite%fmean_sfcw_energy(ipa) = csite%fmean_sfcw_energy(ipa)                       &
                                      / csite%fmean_sfcw_mass(ipa)
         call uint2tl(csite%fmean_sfcw_energy(ipa),csite%fmean_sfcw_temp(ipa)              &
                     ,csite%fmean_sfcw_fliq(ipa))
      else
         csite%fmean_sfcw_mass  (ipa)  = 0.
         csite%fmean_sfcw_depth (ipa)  = 0.
         csite%fmean_sfcw_energy(ipa)  = 0.
         csite%fmean_sfcw_temp  (ipa)  = csite%fmean_soil_temp(mzg,ipa)
         csite%fmean_sfcw_fliq  (ipa)  = csite%fmean_soil_fliq(mzg,ipa)
      end if
      !-----------------------------------------------------------------------------------!
   end do
   return
end subroutine update_patch_thermo_fmean
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will update the derived properties at the site level.                !
!------------------------------------------------------------------------------------------!
subroutine update_site_derived_props(cpoly,census_flag,isi)
  
   use ed_state_vars , only : polygontype  & ! structure
                            , sitetype     & ! structure
                            , patchtype    ! ! structure
   use consts_coms   , only : pio4         ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(polygontype) , target     :: cpoly
   integer           , intent(in) :: census_flag
   integer           , intent(in) :: isi
   !----- Local variables -----------------------------------------------------------------!
   type(sitetype)    , pointer    :: csite
   type(patchtype)   , pointer    :: cpatch
   integer                        :: bdbh
   integer                        :: ipa
   integer                        :: ico
   integer                        :: ipft
   integer                        :: ilu
   !---------------------------------------------------------------------------------------!
   
   !----- Initialise the variables before looping. ----------------------------------------!
   cpoly%basal_area(:,:,isi) = 0.0
   cpoly%agb       (:,:,isi) = 0.0

   csite => cpoly%site(isi)

   !----- Loop over patches. --------------------------------------------------------------!
   do ipa = 1,csite%npatches
      ilu = csite%dist_type(ipa)
      cpatch => csite%patch(ipa)

      !----- Loop over cohorts. -----------------------------------------------------------!
      do ico = 1,cpatch%ncohorts
         ipft = cpatch%pft(ico)

         !----- Update basal area and above-ground biomass. -------------------------------!
         if(census_flag == 0 .or. cpatch%first_census(ico) == 1)then
            bdbh = max(0,min( int(cpatch%dbh(ico) * 0.1), 10)) + 1

            cpoly%basal_area(ipft,bdbh,isi) = cpoly%basal_area(ipft, bdbh,isi)             &
                                            + cpatch%basarea(ico) * cpatch%nplant(ico)     &
                                            * csite%area(ipa)   
            cpoly%agb(ipft,bdbh,isi)        = cpoly%agb(ipft, bdbh,isi)                    &
                                            + cpatch%agb(ico)     * cpatch%nplant(ico)     &
                                            * csite%area(ipa)
         end if
      end do
   end do
   
   return
end subroutine update_site_derived_props
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     The following subroutine finds the polygon averages from site-, patch-, and cohort-  !
! -level properties whose time step is longer than DTLSM (days, months, years).  Fluxes,   !
! meteorological input, thermodynamic properties, and radiation are aggregated in  sub-    !
! -routine aggregate_polygon_fmean.                                                        !
!------------------------------------------------------------------------------------------!
subroutine update_polygon_derived_props(cgrid)
   use ed_state_vars         , only : edtype             & ! structure
                                    , polygontype        & ! structure
                                    , sitetype           & ! structure
                                    , patchtype          ! ! structure
   use soil_coms             , only : soil               & ! intent(in)
                                    , dslz               ! ! intent(in)
   use grid_coms             , only : nzg                & ! intent(in)
                                    , nzs                ! ! intent(in)
   use ed_max_dims           , only : n_dbh              ! ! intent(in)
   use ed_misc_coms          , only : ddbhi              ! ! intent(in)
   use pft_coms              , only : c2n_leaf           & ! intent(in)
                                    , c2n_stem           & ! intent(in)
                                    , c2n_storage        & ! intent(in)
                                    , c2n_recruit        & ! intent(in)
                                    , c2n_slow           & ! intent(in)
                                    , c2n_structural     ! ! intent(in)
   use decomp_coms           , only : cwd_frac           ! ! intent(in)

   implicit none
   !----- Arguments.      -----------------------------------------------------------------!
   type(edtype)         , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype)    , pointer :: cpoly
   type(sitetype)       , pointer :: csite
   type(patchtype)      , pointer :: cpatch
   integer                        :: ipy
   integer                        :: isi
   integer                        :: ipa
   integer                        :: ico
   integer                        :: p
   integer                        :: d
   integer                        :: k
   real                           :: poly_area_i
   real                           :: site_area_i
   real                           :: site_wgt
   real                           :: patch_wgt
   real                           :: rdepth
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!   !
   !---------------------------------------------------------------------------------------!
   !     Please, don't initialise polygon-level (cgrid) variables outside polyloop.  This  !
   ! works in off-line runs, but it causes memory leaks (and crashes) in the coupled runs  !
   ! over the ocean, where cgrid%npolygons can be 0 if one of the sub-domains falls        !
   ! entirely over the ocean.  Thanks!                                                     !
   !---------------------------------------------------------------------------------------!
   ! cgrid%blah = 0. !<<--- This is a bad way of doing, look inside the loop for the
   !                 !      safe way of initialising the variable.
   !---------------------------------------------------------------------------------------!
   polyloop: do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)

      !------------------------------------------------------------------------------------!
      !                                                                                    !
      !     Initialise properties.                                                         !
      !                                                                                    !
      !     This is the right and safe place to initialise polygon-level (cgrid) vari-     !
      ! ables, so in case npolygons is zero this will not cause memory leaks.  I know,     !
      ! this never happens in off-line runs, but it is quite common in coupled runs...     !
      ! Whenever one of the nodes receives a sub-domain where all the points are over the  !
      ! ocean, ED will not assign any polygon in that sub-domain, which means that that    !
      ! node will have 0 polygons, and the variables cannot be allocated.  If you try to   !
      ! access the polygon level variable outside the loop, then the model crashes due to  !
      ! segmentation violation (a bad thing), whereas by putting the variables here both   !
      ! the off-line model and the coupled runs will work, because this loop will be       !
      ! skipped when there is no polygon.                                                  !
      !------------------------------------------------------------------------------------!
      ! cgrid%blah(ipy) = 0.0 ! <<- This way works for all cases. 
      !------------------------------------------------------------------------------------!
      cgrid%nplant              (:,:,ipy) = 0.0
      cgrid%agb                 (:,:,ipy) = 0.0
      cgrid%lai                 (:,:,ipy) = 0.0
      cgrid%wai                 (:,:,ipy) = 0.0
      cgrid%basal_area          (:,:,ipy) = 0.0
      cgrid%bdead               (:,:,ipy) = 0.0
      cgrid%balive              (:,:,ipy) = 0.0
      cgrid%bleaf               (:,:,ipy) = 0.0
      cgrid%broot               (:,:,ipy) = 0.0
      cgrid%bsapwooda           (:,:,ipy) = 0.0
      cgrid%bsapwoodb           (:,:,ipy) = 0.0
      cgrid%bseeds              (:,:,ipy) = 0.0
      cgrid%bstorage            (:,:,ipy) = 0.0
      cgrid%bdead_n             (:,:,ipy) = 0.0
      cgrid%balive_n            (:,:,ipy) = 0.0
      cgrid%bleaf_n             (:,:,ipy) = 0.0
      cgrid%broot_n             (:,:,ipy) = 0.0
      cgrid%bsapwooda_n         (:,:,ipy) = 0.0
      cgrid%bsapwoodb_n         (:,:,ipy) = 0.0
      cgrid%bseeds_n            (:,:,ipy) = 0.0
      cgrid%bstorage_n          (:,:,ipy) = 0.0
      cgrid%leaf_maintenance    (:,:,ipy) = 0.0
      cgrid%root_maintenance    (:,:,ipy) = 0.0
      cgrid%leaf_drop           (:,:,ipy) = 0.0
      cgrid%fast_soil_c             (ipy) = 0.0
      cgrid%slow_soil_c             (ipy) = 0.0
      cgrid%struct_soil_c           (ipy) = 0.0
      cgrid%struct_soil_l           (ipy) = 0.0
      cgrid%cwd_c                   (ipy) = 0.0
      cgrid%fast_soil_n             (ipy) = 0.0
      cgrid%mineral_soil_n          (ipy) = 0.0
      cgrid%cwd_n                   (ipy) = 0.0
      !------------------------------------------------------------------------------------!
      !     Some of these variables are redundant with variables above.  Perhaps we should !
      ! find a single way to report them.                                                  !
      !------------------------------------------------------------------------------------!
      cgrid%Nbiomass_uptake         (ipy) = 0.0
      cgrid%Cleaf_litter_flux       (ipy) = 0.0
      cgrid%Croot_litter_flux       (ipy) = 0.0
      cgrid%Nleaf_litter_flux       (ipy) = 0.0
      cgrid%Nroot_litter_flux       (ipy) = 0.0
      cgrid%Ngross_min              (ipy) = 0.0
      cgrid%Nnet_min                (ipy) = 0.0
      !------------------------------------------------------------------------------------!




      !----- Inverse of this polygon area (it should be always 1.) ------------------------!
      poly_area_i = 1./sum(cpoly%area)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop over sites.                                                               !
      !------------------------------------------------------------------------------------!
      siteloop: do isi=1,cpoly%nsites
         csite => cpoly%site(isi)

         !----- Inverse of this site area (it should be always 1.) ------------------------!
         site_area_i=1./sum(csite%area)
         !---------------------------------------------------------------------------------!


         !----- Site weight. --------------------------------------------------------------!
         site_wgt = cpoly%area(isi) * poly_area_i
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Loop over patches.                                                          !
         !---------------------------------------------------------------------------------!
         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)


            !----- Site weight. -----------------------------------------------------------!
            patch_wgt = csite%area(ipa) * site_area_i * site_wgt
            !------------------------------------------------------------------------------!


            !----- Integrate soil properties. ---------------------------------------------!
            cgrid%fast_soil_c   (ipy) = cgrid%fast_soil_c        (ipy)                     &
                                      + csite%fast_soil_c        (ipa)                     &
                                      * patch_wgt
            cgrid%slow_soil_c   (ipy) = cgrid%slow_soil_c        (ipy)                     &
                                      + csite%slow_soil_c        (ipa)                     &
                                      * patch_wgt
            cgrid%struct_soil_c (ipy) = cgrid%struct_soil_c      (ipy)                     &
                                      + csite%structural_soil_c  (ipa)                     &
                                      * patch_wgt
            cgrid%struct_soil_l (ipy) = cgrid%struct_soil_l      (ipy)                     &
                                      + csite%structural_soil_l  (ipa)                     &
                                      * patch_wgt
            cgrid%fast_soil_n   (ipy) = cgrid%fast_soil_n        (ipy)                     &
                                      + csite%fast_soil_n        (ipa)                     &
                                      * patch_wgt
            cgrid%mineral_soil_n(ipy) = cgrid%mineral_soil_n     (ipy)                     &
                                      + csite%mineralized_soil_n (ipa)                     &
                                      * patch_wgt
            !------------------------------------------------------------------------------!
            !     I am definitely not sure about the way I did CWD.  I just used the same  !
            ! fraction used to compute cwd_rh in soil_rh (soil_respiration.f90) to         !
            ! determine the fraction of the pools that are CWD.                            !
            !------------------------------------------------------------------------------!
            cgrid%cwd_c         (ipy) = cgrid%cwd_c              (ipy)                     &
                                      + ( csite%slow_soil_c      (ipa)                     &
                                        + csite%structural_soil_c(ipa) )                   &
                                      * cwd_frac * patch_wgt
            cgrid%cwd_n         (ipy) = cgrid%cwd_n              (ipy)                     &
                                      + ( csite%slow_soil_c      (ipa) / c2n_slow          &
                                        + csite%structural_soil_c(ipa) / c2n_structural )  &
                                      * cwd_frac * patch_wgt
            !------------------------------------------------------------------------------!



            !----- Zero the root fraction (patch-level diagnostic). -----------------------!
            csite%rootdense(:,ipa) = 0.
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Loop over cohorts.                                                       !
            !------------------------------------------------------------------------------!
            cohortloop: do ico=1,cpatch%ncohorts
               !----- Find the PFT and DBH class to which this cohort belongs. ------------!
               p = cpatch%pft(ico)
               d = max(1,min(n_dbh,ceiling(cpatch%dbh(ico)*ddbhi)))
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Rooting fraction: step 1, find root biomass per cubic meter            !
               !    broot*nplant/rooting_depth   [kg/plant]*[plant/m2]/[m]                 !
               !---------------------------------------------------------------------------!
               rdepth = sum(dslz(cpatch%krdepth(ico):nzg))
               do k=cpatch%krdepth(ico),nzg
                  csite%rootdense(k,ipa) = csite%rootdense(k,ipa)                          &
                                         + cpatch%broot(ico)*cpatch%nplant(ico) / rdepth
               end do
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Integrate cohort-based properties.  Make sure that the polygon-        !
               ! -level gets the right units (i.e., no /plant, but /m2).                   !
               !---------------------------------------------------------------------------!
               cgrid%nplant          (p,d,ipy) = cgrid%nplant          (p,d,ipy)           &
                                               + cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%lai             (p,d,ipy) = cgrid%lai             (p,d,ipy)           &
                                               + cpatch%lai                (ico)
               cgrid%wai             (p,d,ipy) = cgrid%wai             (p,d,ipy)           &
                                               + cpatch%wai                (ico)
               cgrid%agb             (p,d,ipy) = cgrid%agb             (p,d,ipy)           &
                                               + cpatch%agb                (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%basal_area      (p,d,ipy) = cgrid%basal_area      (p,d,ipy)           &
                                               + cpatch%basarea            (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%bdead           (p,d,ipy) = cgrid%bdead           (p,d,ipy)           &
                                               + cpatch%bdead              (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%balive          (p,d,ipy) = cgrid%balive          (p,d,ipy)           &
                                               + cpatch%balive             (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%bleaf           (p,d,ipy) = cgrid%bleaf           (p,d,ipy)           &
                                               + cpatch%bleaf              (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%broot           (p,d,ipy) = cgrid%broot           (p,d,ipy)           &
                                               + cpatch%broot              (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%bsapwooda       (p,d,ipy) = cgrid%bsapwooda       (p,d,ipy)           &
                                               + cpatch%bsapwooda          (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%bsapwoodb       (p,d,ipy) = cgrid%bsapwoodb       (p,d,ipy)           &
                                               + cpatch%bsapwoodb          (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%bseeds          (p,d,ipy) = cgrid%bseeds          (p,d,ipy)           &
                                               + cpatch%bseeds             (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%bstorage        (p,d,ipy) = cgrid%bstorage        (p,d,ipy)           &
                                               + cpatch%bstorage           (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%bdead_n         (p,d,ipy) = cgrid%bdead_n         (p,d,ipy)           &
                                               + cpatch%bdead              (ico)           &
                                               / c2n_stem                    (p)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%balive_n        (p,d,ipy) = cgrid%balive_n        (p,d,ipy)           &
                                               + ( ( cpatch%bleaf          (ico)           &
                                                   + cpatch%broot          (ico) )         &
                                                   / c2n_leaf                (p)           &
                                                 + ( cpatch%bsapwooda      (ico)           &
                                                   + cpatch%bsapwoodb      (ico) )         &
                                                   / c2n_stem                (p)   )       &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%bleaf_n         (p,d,ipy) = cgrid%bleaf_n         (p,d,ipy)           &
                                               + cpatch%bleaf              (ico)           &
                                               / c2n_leaf                    (p)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%broot_n         (p,d,ipy) = cgrid%broot_n         (p,d,ipy)           &
                                               + cpatch%broot              (ico)           &
                                               / c2n_leaf                    (p)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%bsapwooda_n     (p,d,ipy) = cgrid%bsapwooda_n     (p,d,ipy)           &
                                               + cpatch%bsapwooda          (ico)           &
                                               / c2n_stem                    (p)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%bsapwoodb_n     (p,d,ipy) = cgrid%bsapwoodb_n     (p,d,ipy)           &
                                               + cpatch%bsapwoodb          (ico)           &
                                               / c2n_stem                    (p)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%bseeds_n        (p,d,ipy) = cgrid%bseeds_n        (p,d,ipy)           &
                                               + cpatch%bseeds             (ico)           &
                                               / c2n_recruit                 (p)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%bstorage_n      (p,d,ipy) = cgrid%bstorage_n      (p,d,ipy)           &
                                               + cpatch%bstorage           (ico)           &
                                               / c2n_storage                               &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%leaf_maintenance(p,d,ipy) = cgrid%leaf_maintenance(p,d,ipy)           &
                                               + cpatch%leaf_maintenance   (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%root_maintenance(p,d,ipy) = cgrid%root_maintenance(p,d,ipy)           &
                                               + cpatch%root_maintenance   (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               cgrid%leaf_drop       (p,d,ipy) = cgrid%leaf_drop       (p,d,ipy)           &
                                               + cpatch%leaf_drop          (ico)           &
                                               * cpatch%nplant             (ico)           &
                                               * patch_wgt
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !   CHECK! CHECK! CHECK! CHECK! CHECK! CHECK! CHECK! CHECK! CHECK! CHECK!   !
               !---------------------------------------------------------------------------!
               !   These variables were originally in integrate_ed_daily_output_flux, I    !
               ! moved them to here just to be consistent (they are not "dmean"            !
               ! variables).  However, they are slightly different than the original,      !
               ! which were missing nplant (so the results were in not in units of         !
               ! kg/m2/yr, but kg/plant/yr).                                               !
               !---------------------------------------------------------------------------!
               cgrid%Cleaf_litter_flux(ipy) = cgrid%Cleaf_litter_flux(ipy)                 &
                                            + cpatch%leaf_maintenance(ico)                 &
                                            * cpatch%nplant          (ico)                 &
                                            * patch_wgt
               cgrid%Croot_litter_flux(ipy) = cgrid%Croot_litter_flux(ipy)                 &
                                            + cpatch%root_maintenance(ico)                 &
                                            * cpatch%nplant          (ico)                 &
                                            * patch_wgt
               cgrid%Nleaf_litter_flux(ipy) = cgrid%Cleaf_litter_flux(ipy)                 &
                                            + cpatch%leaf_maintenance(ico) / c2n_leaf(p)   &
                                            * cpatch%nplant          (ico)                 &
                                            * patch_wgt
               cgrid%Nroot_litter_flux(ipy) = cgrid%Croot_litter_flux(ipy)                 &
                                            + cpatch%root_maintenance(ico) / c2n_leaf(p)   &
                                            * cpatch%nplant          (ico)                 &
                                            * patch_wgt
               !---------------------------------------------------------------------------!
            end do cohortloop
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Update the patch-related N budget variables.                              !
            !------------------------------------------------------------------------------!
            cgrid%Ngross_min     (ipy) = cgrid%Ngross_min(ipy)                             &
                                       + csite%mineralized_N_input   (ipa)   * patch_wgt
            cgrid%Ngross_min     (ipy) = cgrid%Ngross_min(ipy)                             &
                                       + ( csite%mineralized_N_input (ipa)                 &
                                         - csite%mineralized_N_loss  (ipa) ) * patch_wgt
            cgrid%Nbiomass_uptake(ipy) = cgrid%Ngross_min(ipy)                             &
                                       + csite%total_plant_nitrogen_uptake(ipa)* patch_wgt
            !------------------------------------------------------------------------------!

         end do patchloop
         !---------------------------------------------------------------------------------!
      end do siteloop
      !------------------------------------------------------------------------------------!
   end do polyloop
   !---------------------------------------------------------------------------------------!
   return
end subroutine update_polygon_derived_props
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will read the regular soil moisture and temperature dataset.          !
!------------------------------------------------------------------------------------------!
subroutine read_soil_moist_temp(cgrid,igr)

   use ed_state_vars , only : edtype       & ! structure
                            , polygontype  & ! structure
                            , sitetype     & ! structure
                            , patchtype    ! ! structure
   use soil_coms     , only : soilstate_db & ! intent(in)
                            , soil         & ! intent(in)
                            , slz          ! ! intent(in)
   use consts_coms   , only : wdns         & ! intent(in)
                            , t3ple        ! ! intent(in)
   use therm_lib     , only : cmtl2uext    ! ! function
   use grid_coms     , only : nzg          & ! intent(in)
                            , nzs          & ! intent(in)
                            , ngrids       ! ! intent(in)
   use ed_therm_lib  , only : ed_grndvap   ! ! subroutine
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)      , target     :: cgrid         ! Alias for current ED grid
   integer           , intent(in) :: igr           ! The grid number
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer    :: cpoly          ! Alias for current polygon
   type(sitetype)    , pointer    :: csite          ! Alias for current site
   type(patchtype)   , pointer    :: cpatch         ! Alias for current patch
   integer                        :: ntext          !
   integer                        :: ilat           !
   integer                        :: ilon           !
   integer                        :: ilatf          !
   integer                        :: ilonf          !
   integer                        :: nls            !
   integer                        :: nlsw1          !
   integer                        :: k              !
   integer                        :: ipy            !
   integer                        :: isi            !
   integer                        :: ipa            !
   logical                        :: l1             !
   real                           :: glat           !
   real                           :: glon           !
   real                           :: soil_tempaux   !
   real                           :: tmp1           !
   real                           :: tmp2           !
   real                           :: soilw1         !
   real                           :: soilw2         !
   !----- Local constants.  ---------------------------------------------------------------!
   logical           , parameter  :: harvard_override = .false.
   integer           , parameter  :: nlon = 144
   integer           , parameter  :: nlat = 73
   real              , parameter  :: dlon = 2.5
   real              , parameter  :: dlat = 2.5
   !---------------------------------------------------------------------------------------!

   !----- First thing, check whether the dataset exists and crash the run if it doesn´t. --!
   inquire(file=trim(soilstate_db),exist=l1)
   if (.not.l1) then
      write (unit=*,fmt='(a)') ' Your namelist has ISOILSTATEINIT set to read the initial'
      write (unit=*,fmt='(a)') ' soil moisture and temperature from a file.  The file'
      write (unit=*,fmt='(a)') ' specified by SOILSTATE_DB, however, doesn''t exist!'
      call fatal_error('Soil database '//trim(soilstate_db)//' not found!'                 &
                     &,'read_soil_moist_temp','update_derived_props.f90')
   end if

   open (unit=12,file=trim(soilstate_db),form='formatted',status='old',position='rewind')

   !---- Loop over latitude levels, from north pole then southwards. ----------------------!
   latloop: do ilatf = 1,nlat

      !----- Loop over longitude levels, from Greenwich Meridian then eastwards. ----------!
      lonloop: do ilonf = 1,nlon 
 
         !---------------------------------------------------------------------------------!
         !     Read in reanalysis: two temperatures and moistures, corresponding to        !
         ! different depths.                                                               !
         ! + soilw1, soilw2 are relative porosities and thus range from [0-1].             !
         ! + tmp1, tmp2 are temperature in kelvin.                                         !
         !---------------------------------------------------------------------------------!
         read (unit=12,fmt=*) tmp1,tmp2,soilw1,soilw2

         !----- Make sure the numbers make sense... ---------------------------------------!
         if (tmp1 > 0.0 .and. tmp2 > 0.0 .and. soilw1 > 0.0 .and. soilw2 > 0.0) then

            !----- Loop over land points. -------------------------------------------------!
            polyloop: do ipy=1,cgrid%npolygons
               cpoly => cgrid%polygon(ipy)

               !----- Land point lat, lon. ------------------------------------------------!
               glat = cgrid%lat(ipy)
               glon = cgrid%lon(ipy)
               
               if(glon < 0.0) glon = glon + 360.0
               
               !----- Find reanalysis point corresponding to this land point. -------------!
               if(glat >= 0.0)then
                  ilat = nint((90.0 - glat)/dlat) + 1
               else
                  ilat = nlat - nint((90.0 - abs(glat))/dlat)
               end if
               ilon = int(glon/dlon) + 1
               
               !----- If we are at the right point, fill the array. -----------------------!
               if(ilat == ilatf .and. ilon == ilonf)then

                  !------ Loop over sites and patches. ------------------------------------!
                  siteloop: do isi=1,cpoly%nsites
                     csite => cpoly%site(isi)

                     patchloop: do ipa=1,csite%npatches
                        cpatch => csite%patch(ipa)

                        do k=1,nzg
                           ntext = cpoly%ntext_soil(k,isi)

                           if(abs(slz(k)) < 0.1)then
                              csite%soil_tempk(k,ipa) = tmp1
                              csite%soil_water(k,ipa) = max(soil(ntext)%soilcp             &
                                                           ,soilw1 * soil(ntext)%slmsts)
                           else
                              csite%soil_tempk(k,ipa) = tmp2
                              csite%soil_water(k,ipa) = max(soil(ntext)%soilcp             &
                                                           ,soilw2 * soil(ntext)%slmsts)
                           endif
                           if (csite%soil_tempk(k,ipa) > t3ple) then
                              csite%soil_fracliq(k,ipa) = 1.0
                           elseif (csite%soil_tempk(k,ipa) < t3ple) then
                              csite%soil_fracliq(k,ipa) = 0.0
                           else
                              csite%soil_fracliq(k,ipa) = 0.5
                           end if
                           csite%soil_energy(k,ipa) = cmtl2uext( soil(ntext)%slcpd         &
                                                               , csite%soil_water(k,ipa)   &
                                                               * wdns                      &
                                                               , csite%soil_tempk(k,ipa)   &
                                                               , csite%soil_fracliq(k,ipa))
                        end do


                       !----- Initial condition is with no snow/pond. ---------------------!
                       csite%nlev_sfcwater(ipa)    = 0
                       csite%total_sfcw_depth(ipa) = 0.
                        do k=1,nzs
                           csite%sfcwater_energy (k,ipa) = 0.
                           csite%sfcwater_depth  (k,ipa) = 0.
                           csite%sfcwater_mass   (k,ipa) = 0.
                           csite%sfcwater_tempk  (k,ipa) = csite%soil_tempk(nzg,ipa)
                           csite%sfcwater_fracliq(k,ipa) = csite%sfcwater_fracliq(nzg,ipa)
                        end do

                        if(harvard_override)then
                           csite%soil_tempk(1,ipa)     = 277.6
                           csite%soil_tempk(2:4,ipa)   = 276.0
                           csite%soil_energy(1,ipa)    =   1.5293664e8
                           csite%soil_energy(2,ipa)    =   1.4789957e8
                           csite%soil_energy(3:4,ipa)  =   1.4772002e8
                           csite%soil_water(1:4,ipa)   =   0.41595e+0
                           csite%soil_fracliq(1:4,ipa) =   1.0
                        endif
                        
                        !----- Compute the ground specific humidity. ----------------------!
                        ntext = cpoly%ntext_soil(k,isi)
                        nls   = csite%nlev_sfcwater(ipa)
                        nlsw1 = max(1,nls)
                        call ed_grndvap(nls,ntext,csite%soil_water(nzg,ipa)                &
                                       ,csite%soil_tempk(nzg,ipa)                          &
                                       ,csite%soil_fracliq(nzg,ipa)                        &
                                       ,csite%sfcwater_tempk(nlsw1,ipa)                    &
                                       ,csite%sfcwater_fracliq(nlsw1,ipa)                  &
                                       ,csite%snowfac(ipa),csite%can_prss(ipa)             &
                                       ,csite%can_shv(ipa),csite%ground_shv(ipa)           &
                                       ,csite%ground_ssh(ipa),csite%ground_temp(ipa)       &
                                       ,csite%ground_fliq(ipa),csite%ggsoil(ipa))

                     end do patchloop
                  end do siteloop
               end if
            end do polyloop
         end if
      end do lonloop
   end do latloop

   close(unit=12,status='keep')
   return

end subroutine read_soil_moist_temp
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will convert the integrated number of time steps in steps/day, then  !
! it will update the monthly mean workload.                                                !
!------------------------------------------------------------------------------------------!
subroutine update_workload(cgrid)
   use ed_state_vars, only : edtype        ! ! structure
   use ed_misc_coms , only : current_time  & ! intent(in)
                           , simtime       ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype), target     :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(simtime)            :: lastmonth
   integer                  :: lmon
   integer                  :: ipy
   real                     :: ndaysi
   !---------------------------------------------------------------------------------------!

   !----- Find last month information. ----------------------------------------------------!
   call lastmonthdate(current_time,lastmonth,ndaysi)
   lmon = lastmonth%month

   !---------------------------------------------------------------------------------------!
   !     Loop over all polygons, normalise the workload, then copy it to the corresponding !
   ! month.  Then copy the scratch column (13) to the appropriate month, and reset it.     !
   !---------------------------------------------------------------------------------------!
   do ipy=1,cgrid%npolygons
      cgrid%workload(lmon,ipy) = cgrid%workload(13,ipy) * ndaysi
      cgrid%workload(13,ipy)   = 0.
   end do

   return
end subroutine update_workload
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine updates the model time.  It was moved to here to reduce the number  !
! of routines that must be written twice (off-line and coupled model).                     !
!------------------------------------------------------------------------------------------!
subroutine update_model_time_dm(ctime,dtlong)

   use ed_misc_coms, only : simtime ! ! variable type
   use consts_coms , only : day_sec & ! intent(in)
                          , hr_sec  & ! intent(in)
                          , min_sec ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(simtime)         , intent(inout) :: ctime  ! Current time
   real                  , intent(in)    :: dtlong ! Model time step
   !----- Local variables. ----------------------------------------------------------------!
   integer               , dimension(12) :: daymax
   !----- External functions. -------------------------------------------------------------!
   logical               , external      :: isleap ! This function check for leap years
   !---------------------------------------------------------------------------------------!



   !----- Assume that the year is not leap. -----------------------------------------------!
   daymax=(/31,28,31,30,31,30,31,31,30,31,30,31/)
   !---------------------------------------------------------------------------------------!


   !----- Update the time. ----------------------------------------------------------------!
   ctime%time = ctime%time + dtlong
   !---------------------------------------------------------------------------------------!

   if (ctime%time >= day_sec)then

      !----- If time is greater than one day, update the day. -----------------------------!
      ctime%time = ctime%time - day_sec
      ctime%date = ctime%date + 1

      !----- If the year is leap, correct February's number of days. ----------------------!
      if (isleap(ctime%year)) daymax(2) = 29

      !----- If we have reached the end of the month, update the month. -------------------!
      if (ctime%date > daymax(ctime%month)) then
         ctime%date  = 1
         ctime%month = ctime%month + 1

         !------ If we have reached the end of the year, update the year. -----------------!
         if (ctime%month == 13) then
            ctime%month = 1
            ctime%year = ctime%year + 1
         end if
      end if

   elseif (ctime%time < 0.0) then
      !----- Time is negative, go back one day. -------------------------------------------!
      ctime%time = ctime%time + day_sec
      ctime%date = ctime%date - 1

      !----- Day is zero, which means we must go to the previous month. -------------------!
      if (ctime%date == 0) then
         ctime%month = ctime%month - 1

         !----- Month is zero, which means that we must go to the previous year. ----------!
         if (ctime%month == 0) then
            ctime%month = 12
            ctime%year = ctime%year - 1
         end if

         !---------------------------------------------------------------------------------!
         !     Fix the month in case it it a leap year, and make the day the last of the   !
         ! previous month.                                                                 !
         !---------------------------------------------------------------------------------!
         if (isleap(ctime%year)) daymax(2) = 29
         ctime%date = daymax(ctime%month)
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !----- Update the hours, minutes, and seconds. -----------------------------------------!
   ctime%hour = floor(ctime%time / hr_sec)
   ctime%min  = floor((ctime%time - real(ctime%hour)*hr_sec)/min_sec)
   ctime%sec  = floor(ctime%time - real(ctime%hour)*hr_sec - real(ctime%min)*min_sec)
   !---------------------------------------------------------------------------------------!


   return
end subroutine update_model_time_dm
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!       This subroutine scales all the "extensive" cohort properties when cohort area or   !
! population changes.                                                                      !
! IMPORTANT: Only cohort-level variables that have units per area (m2 ground) should be    !
!            rescaled.  Variables whose units are per plant, m2 leaf, or m2 wood           !
!            SHOULD NOT be included here.                                                  !
!------------------------------------------------------------------------------------------!
subroutine update_cohort_extensive_props(cpatch,aco,zco,mult)
   use ed_state_vars, only : patchtype    ! ! structure
   use ed_misc_coms , only : writing_long & ! intent(in)
                           , writing_eorq & ! intent(in)
                           , writing_dcyc ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(patchtype), target     :: cpatch    ! Current patch
   integer        , intent(in) :: aco       ! First cohort to be rescaled
   integer        , intent(in) :: zco       ! Last  cohort to be rescaled
   real           , intent(in) :: mult      ! Scale factor
   !----- Local variables. ----------------------------------------------------------------!
   integer                     :: ico       ! Cohort counter
   real                        :: mult_2    ! Square of the scale factor
   !---------------------------------------------------------------------------------------!


   !----- Set up the scale factor for monthly mean sum of squares. ------------------------!
   mult_2 = mult * mult
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Loop over all cohorts.                                                            !
   !---------------------------------------------------------------------------------------!
   cohortloop: do ico=aco,zco
      cpatch%lai                (ico) = cpatch%lai                (ico) * mult
      cpatch%wai                (ico) = cpatch%wai                (ico) * mult
      cpatch%nplant             (ico) = cpatch%nplant             (ico) * mult
      cpatch%today_gpp          (ico) = cpatch%today_gpp          (ico) * mult
      cpatch%today_nppleaf      (ico) = cpatch%today_nppleaf      (ico) * mult
      cpatch%today_nppfroot     (ico) = cpatch%today_nppfroot     (ico) * mult
      cpatch%today_nppsapwood   (ico) = cpatch%today_nppsapwood   (ico) * mult
      cpatch%today_nppcroot     (ico) = cpatch%today_nppcroot     (ico) * mult
      cpatch%today_nppseeds     (ico) = cpatch%today_nppseeds     (ico) * mult
      cpatch%today_nppwood      (ico) = cpatch%today_nppwood      (ico) * mult
      cpatch%today_nppdaily     (ico) = cpatch%today_nppdaily     (ico) * mult
      cpatch%today_gpp_pot      (ico) = cpatch%today_gpp_pot      (ico) * mult
      cpatch%today_gpp_lightmax (ico) = cpatch%today_gpp_lightmax (ico) * mult
      cpatch%today_gpp_moistmax (ico) = cpatch%today_gpp_moistmax (ico) * mult
      cpatch%today_leaf_resp    (ico) = cpatch%today_leaf_resp    (ico) * mult
      cpatch%today_root_resp    (ico) = cpatch%today_root_resp    (ico) * mult
      cpatch%gpp                (ico) = cpatch%gpp                (ico) * mult
      cpatch%leaf_respiration   (ico) = cpatch%leaf_respiration   (ico) * mult
      cpatch%root_respiration   (ico) = cpatch%root_respiration   (ico) * mult
      cpatch%monthly_dndt       (ico) = cpatch%monthly_dndt       (ico) * mult
      cpatch%leaf_water         (ico) = cpatch%leaf_water         (ico) * mult
      cpatch%leaf_hcap          (ico) = cpatch%leaf_hcap          (ico) * mult
      cpatch%leaf_energy        (ico) = cpatch%leaf_energy        (ico) * mult
      cpatch%wood_water         (ico) = cpatch%wood_water         (ico) * mult
      cpatch%wood_hcap          (ico) = cpatch%wood_hcap          (ico) * mult
      cpatch%wood_energy        (ico) = cpatch%wood_energy        (ico) * mult
      !----- Crown area shall not exceed 1. -----------------------------------------------!
      cpatch%crown_area         (ico) = min(1.,cpatch%crown_area  (ico) * mult)
      !----- Fast-scale means. ------------------------------------------------------------!
      cpatch%fmean_leaf_energy   (ico) = cpatch%fmean_leaf_energy   (ico) * mult
      cpatch%fmean_leaf_water    (ico) = cpatch%fmean_leaf_water    (ico) * mult
      cpatch%fmean_leaf_hcap     (ico) = cpatch%fmean_leaf_hcap     (ico) * mult
      cpatch%fmean_wood_energy   (ico) = cpatch%fmean_wood_energy   (ico) * mult
      cpatch%fmean_wood_water    (ico) = cpatch%fmean_wood_water    (ico) * mult
      cpatch%fmean_wood_hcap     (ico) = cpatch%fmean_wood_hcap     (ico) * mult
      cpatch%fmean_water_supply  (ico) = cpatch%fmean_water_supply  (ico) * mult
      cpatch%fmean_par_l         (ico) = cpatch%fmean_par_l         (ico) * mult
      cpatch%fmean_par_l_beam    (ico) = cpatch%fmean_par_l_beam    (ico) * mult
      cpatch%fmean_par_l_diff    (ico) = cpatch%fmean_par_l_diff    (ico) * mult
      cpatch%fmean_rshort_l      (ico) = cpatch%fmean_rshort_l      (ico) * mult
      cpatch%fmean_rlong_l       (ico) = cpatch%fmean_rlong_l       (ico) * mult
      cpatch%fmean_sensible_lc   (ico) = cpatch%fmean_sensible_lc   (ico) * mult
      cpatch%fmean_vapor_lc      (ico) = cpatch%fmean_vapor_lc      (ico) * mult
      cpatch%fmean_transp        (ico) = cpatch%fmean_transp        (ico) * mult
      cpatch%fmean_intercepted_al(ico) = cpatch%fmean_intercepted_al(ico) * mult
      cpatch%fmean_wshed_lg      (ico) = cpatch%fmean_wshed_lg      (ico) * mult
      cpatch%fmean_rshort_w      (ico) = cpatch%fmean_rshort_w      (ico) * mult
      cpatch%fmean_rlong_w       (ico) = cpatch%fmean_rlong_w       (ico) * mult
      cpatch%fmean_rad_profile (:,ico) = cpatch%fmean_rad_profile (:,ico) * mult
      cpatch%fmean_sensible_wc   (ico) = cpatch%fmean_sensible_wc   (ico) * mult
      cpatch%fmean_vapor_wc      (ico) = cpatch%fmean_vapor_wc      (ico) * mult
      cpatch%fmean_intercepted_aw(ico) = cpatch%fmean_intercepted_aw(ico) * mult
      cpatch%fmean_wshed_wg      (ico) = cpatch%fmean_wshed_wg      (ico) * mult
      !----- Daily means. -----------------------------------------------------------------!
      if (writing_long) then
         cpatch%dmean_leaf_energy   (ico) = cpatch%dmean_leaf_energy   (ico) * mult
         cpatch%dmean_leaf_water    (ico) = cpatch%dmean_leaf_water    (ico) * mult
         cpatch%dmean_leaf_hcap     (ico) = cpatch%dmean_leaf_hcap     (ico) * mult
         cpatch%dmean_wood_energy   (ico) = cpatch%dmean_wood_energy   (ico) * mult
         cpatch%dmean_wood_water    (ico) = cpatch%dmean_wood_water    (ico) * mult
         cpatch%dmean_wood_hcap     (ico) = cpatch%dmean_wood_hcap     (ico) * mult
         cpatch%dmean_water_supply  (ico) = cpatch%dmean_water_supply  (ico) * mult
         cpatch%dmean_par_l         (ico) = cpatch%dmean_par_l         (ico) * mult
         cpatch%dmean_par_l_beam    (ico) = cpatch%dmean_par_l_beam    (ico) * mult
         cpatch%dmean_par_l_diff    (ico) = cpatch%dmean_par_l_diff    (ico) * mult
         cpatch%dmean_rshort_l      (ico) = cpatch%dmean_rshort_l      (ico) * mult
         cpatch%dmean_rlong_l       (ico) = cpatch%dmean_rlong_l       (ico) * mult
         cpatch%dmean_sensible_lc   (ico) = cpatch%dmean_sensible_lc   (ico) * mult
         cpatch%dmean_vapor_lc      (ico) = cpatch%dmean_vapor_lc      (ico) * mult
         cpatch%dmean_transp        (ico) = cpatch%dmean_transp        (ico) * mult
         cpatch%dmean_intercepted_al(ico) = cpatch%dmean_intercepted_al(ico) * mult
         cpatch%dmean_wshed_lg      (ico) = cpatch%dmean_wshed_lg      (ico) * mult
         cpatch%dmean_rshort_w      (ico) = cpatch%dmean_rshort_w      (ico) * mult
         cpatch%dmean_rlong_w       (ico) = cpatch%dmean_rlong_w       (ico) * mult
         cpatch%dmean_rad_profile (:,ico) = cpatch%dmean_rad_profile (:,ico) * mult
         cpatch%dmean_sensible_wc   (ico) = cpatch%dmean_sensible_wc   (ico) * mult
         cpatch%dmean_vapor_wc      (ico) = cpatch%dmean_vapor_wc      (ico) * mult
         cpatch%dmean_intercepted_aw(ico) = cpatch%dmean_intercepted_aw(ico) * mult
         cpatch%dmean_wshed_wg      (ico) = cpatch%dmean_wshed_wg      (ico) * mult
      end if
      !----- Monthly means. ---------------------------------------------------------------!
      if (writing_eorq) then
         cpatch%mmean_lai           (ico) = cpatch%mmean_lai           (ico) * mult
         cpatch%mmean_leaf_energy   (ico) = cpatch%mmean_leaf_energy   (ico) * mult
         cpatch%mmean_leaf_water    (ico) = cpatch%mmean_leaf_water    (ico) * mult
         cpatch%mmean_leaf_hcap     (ico) = cpatch%mmean_leaf_hcap     (ico) * mult
         cpatch%mmean_wood_energy   (ico) = cpatch%mmean_wood_energy   (ico) * mult
         cpatch%mmean_wood_water    (ico) = cpatch%mmean_wood_water    (ico) * mult
         cpatch%mmean_wood_hcap     (ico) = cpatch%mmean_wood_hcap     (ico) * mult
         cpatch%mmean_water_supply  (ico) = cpatch%mmean_water_supply  (ico) * mult
         cpatch%mmean_par_l         (ico) = cpatch%mmean_par_l         (ico) * mult
         cpatch%mmean_par_l_beam    (ico) = cpatch%mmean_par_l_beam    (ico) * mult
         cpatch%mmean_par_l_diff    (ico) = cpatch%mmean_par_l_diff    (ico) * mult
         cpatch%mmean_rshort_l      (ico) = cpatch%mmean_rshort_l      (ico) * mult
         cpatch%mmean_rlong_l       (ico) = cpatch%mmean_rlong_l       (ico) * mult
         cpatch%mmean_sensible_lc   (ico) = cpatch%mmean_sensible_lc   (ico) * mult
         cpatch%mmean_vapor_lc      (ico) = cpatch%mmean_vapor_lc      (ico) * mult
         cpatch%mmean_transp        (ico) = cpatch%mmean_transp        (ico) * mult
         cpatch%mmean_intercepted_al(ico) = cpatch%mmean_intercepted_al(ico) * mult
         cpatch%mmean_wshed_lg      (ico) = cpatch%mmean_wshed_lg      (ico) * mult
         cpatch%mmean_rshort_w      (ico) = cpatch%mmean_rshort_w      (ico) * mult
         cpatch%mmean_rlong_w       (ico) = cpatch%mmean_rlong_w       (ico) * mult
         cpatch%mmean_rad_profile (:,ico) = cpatch%mmean_rad_profile (:,ico) * mult
         cpatch%mmean_sensible_wc   (ico) = cpatch%mmean_sensible_wc   (ico) * mult
         cpatch%mmean_vapor_wc      (ico) = cpatch%mmean_vapor_wc      (ico) * mult
         cpatch%mmean_intercepted_aw(ico) = cpatch%mmean_intercepted_aw(ico) * mult
         cpatch%mmean_wshed_wg      (ico) = cpatch%mmean_wshed_wg      (ico) * mult
         cpatch%mmsqu_gpp           (ico) = cpatch%mmsqu_gpp           (ico) * mult_2
         cpatch%mmsqu_npp           (ico) = cpatch%mmsqu_npp           (ico) * mult_2
         cpatch%mmsqu_plresp        (ico) = cpatch%mmsqu_plresp        (ico) * mult_2
         cpatch%mmsqu_sensible_lc   (ico) = cpatch%mmsqu_sensible_lc   (ico) * mult_2
         cpatch%mmsqu_vapor_lc      (ico) = cpatch%mmsqu_vapor_lc      (ico) * mult_2
         cpatch%mmsqu_transp        (ico) = cpatch%mmsqu_transp        (ico) * mult_2
         cpatch%mmsqu_sensible_wc   (ico) = cpatch%mmsqu_sensible_wc   (ico) * mult_2
         cpatch%mmsqu_vapor_wc      (ico) = cpatch%mmsqu_vapor_wc      (ico) * mult_2
      end if
      !----- Mean diel. -------------------------------------------------------------------!
      if (writing_dcyc) then
         cpatch%qmean_leaf_energy   (:,ico) = cpatch%qmean_leaf_energy   (:,ico) * mult
         cpatch%qmean_leaf_water    (:,ico) = cpatch%qmean_leaf_water    (:,ico) * mult
         cpatch%qmean_leaf_hcap     (:,ico) = cpatch%qmean_leaf_hcap     (:,ico) * mult
         cpatch%qmean_wood_energy   (:,ico) = cpatch%qmean_wood_energy   (:,ico) * mult
         cpatch%qmean_wood_water    (:,ico) = cpatch%qmean_wood_water    (:,ico) * mult
         cpatch%qmean_wood_hcap     (:,ico) = cpatch%qmean_wood_hcap     (:,ico) * mult
         cpatch%qmean_water_supply  (:,ico) = cpatch%qmean_water_supply  (:,ico) * mult
         cpatch%qmean_par_l         (:,ico) = cpatch%qmean_par_l         (:,ico) * mult
         cpatch%qmean_par_l_beam    (:,ico) = cpatch%qmean_par_l_beam    (:,ico) * mult
         cpatch%qmean_par_l_diff    (:,ico) = cpatch%qmean_par_l_diff    (:,ico) * mult
         cpatch%qmean_rshort_l      (:,ico) = cpatch%qmean_rshort_l      (:,ico) * mult
         cpatch%qmean_rlong_l       (:,ico) = cpatch%qmean_rlong_l       (:,ico) * mult
         cpatch%qmean_sensible_lc   (:,ico) = cpatch%qmean_sensible_lc   (:,ico) * mult
         cpatch%qmean_vapor_lc      (:,ico) = cpatch%qmean_vapor_lc      (:,ico) * mult
         cpatch%qmean_transp        (:,ico) = cpatch%qmean_transp        (:,ico) * mult
         cpatch%qmean_intercepted_al(:,ico) = cpatch%qmean_intercepted_al(:,ico) * mult
         cpatch%qmean_wshed_lg      (:,ico) = cpatch%qmean_wshed_lg      (:,ico) * mult
         cpatch%qmean_rshort_w      (:,ico) = cpatch%qmean_rshort_w      (:,ico) * mult
         cpatch%qmean_rlong_w       (:,ico) = cpatch%qmean_rlong_w       (:,ico) * mult
         cpatch%qmean_rad_profile (:,:,ico) = cpatch%qmean_rad_profile (:,:,ico) * mult
         cpatch%qmean_sensible_wc   (:,ico) = cpatch%qmean_sensible_wc   (:,ico) * mult
         cpatch%qmean_vapor_wc      (:,ico) = cpatch%qmean_vapor_wc      (:,ico) * mult
         cpatch%qmean_intercepted_aw(:,ico) = cpatch%qmean_intercepted_aw(:,ico) * mult
         cpatch%qmean_wshed_wg      (:,ico) = cpatch%qmean_wshed_wg      (:,ico) * mult
         cpatch%qmsqu_gpp           (:,ico) = cpatch%qmsqu_gpp           (:,ico) * mult_2
         cpatch%qmsqu_npp           (:,ico) = cpatch%qmsqu_npp           (:,ico) * mult_2
         cpatch%qmsqu_plresp        (:,ico) = cpatch%qmsqu_plresp        (:,ico) * mult_2
         cpatch%qmsqu_sensible_lc   (:,ico) = cpatch%qmsqu_sensible_lc   (:,ico) * mult_2
         cpatch%qmsqu_vapor_lc      (:,ico) = cpatch%qmsqu_vapor_lc      (:,ico) * mult_2
         cpatch%qmsqu_transp        (:,ico) = cpatch%qmsqu_transp        (:,ico) * mult_2
         cpatch%qmsqu_sensible_wc   (:,ico) = cpatch%qmsqu_sensible_wc   (:,ico) * mult_2
         cpatch%qmsqu_vapor_wc      (:,ico) = cpatch%qmsqu_vapor_wc      (:,ico) * mult_2
      end if
      !------------------------------------------------------------------------------------!

   end do cohortloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine update_cohort_extensive_props
!==========================================================================================!
!==========================================================================================!
