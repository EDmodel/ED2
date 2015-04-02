!==========================================================================================!
!==========================================================================================!
!     This module contains sub-routines used to create the averaged output variables.      !
!------------------------------------------------------------------------------------------!
module average_utils

   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !                           |----------------------------------|                        !
   !                           |** FREQUENT AVERAGE SUBROUTINES **|                        !
   !                           |----------------------------------|                        !
   !=======================================================================================!
   !=======================================================================================!
   !     The following subroutine finds the polygon averages from site-, patch-, and       !
   ! cohort-level properties that have fmean variables associated.                         !
   !---------------------------------------------------------------------------------------!
   subroutine aggregate_polygon_fmean(cgrid)
      use ed_state_vars         , only : edtype             & ! structure
                                       , polygontype        & ! structure
                                       , sitetype           & ! structure
                                       , patchtype          ! ! structure
      use grid_coms             , only : ngrids             & ! intent(in)
                                       , nzg                & ! intent(in)
                                       , nzs                ! ! intent(in)
      use consts_coms           , only : wdns               & ! intent(in)
                                       , t00                ! ! intent(in)
      use ed_misc_coms          , only : frqsum             ! ! intent(in)
      use therm_lib             , only : uextcm2tl          & ! subroutine
                                       , uint2tl            & ! subroutine
                                       , idealdenssh        & ! function
                                       , press2exner        & ! function
                                       , extheta2temp       ! ! function
      use soil_coms             , only : tiny_sfcwater_mass & ! intent(in)
                                       , isoilbc            & ! intent(in)
                                       , soil               & ! intent(in)
                                       , dslz               ! ! intent(in)
      use ed_max_dims           , only : n_pft              ! ! intent(in)


      implicit none
      !----- Arguments.      --------------------------------------------------------------!
      type(edtype)         , target  :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype)    , pointer :: cpoly
      type(sitetype)       , pointer :: csite
      type(patchtype)      , pointer :: cpatch
      real, dimension(nzg)           :: cgrid_fmean_soil_hcap
      integer                        :: ipy
      integer                        :: isi
      integer                        :: ipa
      integer                        :: ico
      integer                        :: k
      integer                        :: nsoil
      real                           :: site_area_i
      real                           :: poly_area_i
      real                           :: poly_lai
      real                           :: poly_wai
      real                           :: poly_nplant
      real                           :: site_wgt
      real                           :: patch_wgt
      real                           :: skin_energy
      real                           :: skin_water
      real                           :: skin_hcap
      real                           :: skin_fliq
      real                           :: dslzsum_i
      real                           :: can_exner
      real                           :: atm_exner
      !------------------------------------------------------------------------------------!





      !------------------------------------------------------------------------------------!
      !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !
      !------------------------------------------------------------------------------------!
      !     Please, don't initialise polygon-level (cgrid) variables outside polyloop.     !
      ! This works in off-line runs, but it causes memory leaks (and crashes) in the       !
      ! coupled runs over the ocean, where cgrid%npolygons can be 0 if one of the sub-     !
      ! -domains falls entirely over the ocean.  Thanks!                                   !
      !------------------------------------------------------------------------------------!
      ! cgrid%blah = 0. !<<--- This is a bad way of doing, look inside the loop for the
      !                 !      safe way of initialising the variable.
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         !---------------------------------------------------------------------------------!
         !     This is the right and safe place to initialise polygon-level (cgrid) vari-  !
         ! ables, so in case npolygons is zero this will not cause memory leaks.  I know,  !
         ! this never happens in off-line runs, but it is quite common in coupled runs...  !
         ! Whenever one of the nodes receives a sub-domain where all the points are over   !
         ! the ocean, ED will not assign any polygon in that sub-domain, which means that  !
         ! that node will have 0 polygons, and the variables cannot be allocated.  If you  !
         ! try to access the polygon level variable outside the loop, then the model       !
         ! crashes due to segmentation violation (a bad thing), whereas by putting the     !
         ! variables here both the off-line model and the coupled runs will work, because  !
         ! this loop will be skipped when there is no polygon.                             !
         !---------------------------------------------------------------------------------!
         ! cgrid%blah(ipy) = 0. ! <<- This way works for all cases. 
         !---------------------------------------------------------------------------------!


         !----- Inverse of this polygon area (it should be always 1.) ---------------------!
         poly_area_i = 1./sum(cpoly%area)
         !---------------------------------------------------------------------------------!

         !----- Re-set some support variables. --------------------------------------------!
         poly_lai                 = 0.0
         poly_wai                 = 0.0
         poly_nplant              = 0.0
         cgrid_fmean_soil_hcap(:) = 0.0
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Loop over sites.                                                            !
         !---------------------------------------------------------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            !----- Inverse of this site area (it should be always 1.) ---------------------!
            site_area_i=1./sum(csite%area)
            !------------------------------------------------------------------------------!


            !----- Site weight. -----------------------------------------------------------!
            site_wgt = cpoly%area(isi) * poly_area_i
            !------------------------------------------------------------------------------!

            !----- Inverse of the soil depth. ---------------------------------------------!
            dslzsum_i = 1./ sum(dslz(cpoly%lsl(isi):nzg))
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     Loop over patches.                                                       !
            !------------------------------------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)


               !----- Site weight. --------------------------------------------------------!
               patch_wgt = csite%area(ipa) * site_area_i * site_wgt
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !     Loop over cohorts.                                                    !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts
                  !------------------------------------------------------------------------!
                  !    Aggregate AIs and density, they may be used to normalise averages.  !
                  !------------------------------------------------------------------------!
                  poly_nplant = poly_nplant + cpatch%nplant(ico) * patch_wgt
                  poly_lai    = poly_lai    + cpatch%lai   (ico) * patch_wgt
                  poly_wai    = poly_wai    + cpatch%wai   (ico) * patch_wgt
                  !------------------------------------------------------------------------!


                  !----- Aggregate the other properties. ----------------------------------!
                  cgrid%fmean_gpp           (ipy) = cgrid%fmean_gpp            (ipy)       &
                                                  + cpatch%fmean_gpp           (ico)       &
                                                  * cpatch%nplant              (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_npp           (ipy) = cgrid%fmean_npp            (ipy)       &
                                                  + cpatch%fmean_npp           (ico)       &
                                                  * cpatch%nplant              (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_leaf_resp     (ipy) = cgrid%fmean_leaf_resp      (ipy)       &
                                                  + cpatch%fmean_leaf_resp     (ico)       &
                                                  * cpatch%nplant              (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_root_resp     (ipy) = cgrid%fmean_root_resp      (ipy)       &
                                                  + cpatch%fmean_root_resp     (ico)       &
                                                  * cpatch%nplant              (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_growth_resp   (ipy) = cgrid%fmean_growth_resp    (ipy)       &
                                                  + cpatch%fmean_growth_resp   (ico)       &
                                                  * cpatch%nplant              (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_storage_resp  (ipy) = cgrid%fmean_storage_resp   (ipy)       &
                                                  + cpatch%fmean_storage_resp  (ico)       &
                                                  * cpatch%nplant              (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_vleaf_resp    (ipy) = cgrid%fmean_vleaf_resp     (ipy)       &
                                                  + cpatch%fmean_vleaf_resp    (ico)       &
                                                  * cpatch%nplant              (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_plresp        (ipy) = cgrid%fmean_plresp         (ipy)       &
                                                  + cpatch%fmean_plresp        (ico)       &
                                                  * cpatch%nplant              (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_leaf_energy   (ipy) = cgrid%fmean_leaf_energy    (ipy)       &
                                                  + cpatch%fmean_leaf_energy   (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_leaf_water    (ipy) = cgrid%fmean_leaf_water     (ipy)       &
                                                  + cpatch%fmean_leaf_water    (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_leaf_hcap     (ipy) = cgrid%fmean_leaf_hcap      (ipy)       &
                                                  + cpatch%fmean_leaf_hcap     (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_leaf_vpdef    (ipy) = cgrid%fmean_leaf_vpdef     (ipy)       &
                                                  + cpatch%fmean_leaf_vpdef    (ico)       &
                                                  * cpatch%lai                 (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_leaf_gsw      (ipy) = cgrid%fmean_leaf_gsw       (ipy)       &
                                                  + cpatch%fmean_leaf_gsw      (ico)       &
                                                  * cpatch%lai                 (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_leaf_gbw      (ipy) = cgrid%fmean_leaf_gbw       (ipy)       &
                                                  + cpatch%fmean_leaf_gbw      (ico)       &
                                                  * cpatch%lai                 (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_wood_energy   (ipy) = cgrid%fmean_wood_energy    (ipy)       &
                                                  + cpatch%fmean_wood_energy   (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_wood_water    (ipy) = cgrid%fmean_wood_water     (ipy)       &
                                                  + cpatch%fmean_wood_water    (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_wood_hcap     (ipy) = cgrid%fmean_wood_hcap      (ipy)       &
                                                  + cpatch%fmean_wood_hcap     (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_wood_gbw      (ipy) = cgrid%fmean_wood_gbw       (ipy)       &
                                                  + cpatch%fmean_wood_gbw      (ico)       &
                                                  * cpatch%wai                 (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_fs_open       (ipy) = cgrid%fmean_fs_open        (ipy)       &
                                                  + cpatch%fmean_fs_open       (ico)       &
                                                  * cpatch%lai                 (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_fsw           (ipy) = cgrid%fmean_fsw            (ipy)       &
                                                  + cpatch%fmean_fsw           (ico)       &
                                                  * cpatch%lai                 (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_fsn           (ipy) = cgrid%fmean_fsn            (ipy)       &
                                                  + cpatch%fmean_fsn           (ico)       &
                                                  * cpatch%lai                 (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_a_light       (ipy) = cgrid%fmean_a_light        (ipy)       &
                                                  + cpatch%fmean_a_light       (ico)       &
                                                  * cpatch%lai                 (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_a_rubp        (ipy) = cgrid%fmean_a_rubp         (ipy)       &
                                                  + cpatch%fmean_a_rubp        (ico)       &
                                                  * cpatch%lai                 (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_a_co2         (ipy) = cgrid%fmean_a_co2          (ipy)       &
                                                  + cpatch%fmean_a_co2         (ico)       &
                                                  * cpatch%lai                 (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_psi_open      (ipy) = cgrid%fmean_psi_open       (ipy)       &
                                                  + cpatch%fmean_psi_open      (ico)       &
                                                  * cpatch%lai                 (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_psi_closed    (ipy) = cgrid%fmean_psi_closed     (ipy)       &
                                                  + cpatch%fmean_psi_closed    (ico)       &
                                                  * cpatch%lai                 (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_water_supply  (ipy) = cgrid%fmean_water_supply   (ipy)       &
                                                  + cpatch%fmean_water_supply  (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_par_l         (ipy) = cgrid%fmean_par_l          (ipy)       &
                                                  + cpatch%fmean_par_l         (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_par_l_beam    (ipy) = cgrid%fmean_par_l_beam     (ipy)       &
                                                  + cpatch%fmean_par_l_beam    (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_par_l_diff    (ipy) = cgrid%fmean_par_l_diff     (ipy)       &
                                                  + cpatch%fmean_par_l_diff    (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_rshort_l      (ipy) = cgrid%fmean_rshort_l       (ipy)       &
                                                  + cpatch%fmean_rshort_l      (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_rlong_l       (ipy) = cgrid%fmean_rlong_l        (ipy)       &
                                                  + cpatch%fmean_rlong_l       (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_sensible_lc   (ipy) = cgrid%fmean_sensible_lc    (ipy)       &
                                                  + cpatch%fmean_sensible_lc   (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_vapor_lc      (ipy) = cgrid%fmean_vapor_lc       (ipy)       &
                                                  + cpatch%fmean_vapor_lc      (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_transp        (ipy) = cgrid%fmean_transp         (ipy)       &
                                                  + cpatch%fmean_transp        (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_intercepted_al(ipy) = cgrid%fmean_intercepted_al (ipy)       &
                                                  + cpatch%fmean_intercepted_al(ico)       &
                                                  * patch_wgt
                  cgrid%fmean_wshed_lg      (ipy) = cgrid%fmean_wshed_lg       (ipy)       &
                                                  + cpatch%fmean_wshed_lg      (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_rshort_w      (ipy) = cgrid%fmean_rshort_w       (ipy)       &
                                                  + cpatch%fmean_rshort_w      (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_rlong_w       (ipy) = cgrid%fmean_rlong_w        (ipy)       &
                                                  + cpatch%fmean_rlong_w       (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_sensible_wc   (ipy) = cgrid%fmean_sensible_wc    (ipy)       &
                                                  + cpatch%fmean_sensible_wc   (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_vapor_wc      (ipy) = cgrid%fmean_vapor_wc       (ipy)       &
                                                  + cpatch%fmean_vapor_wc      (ico)       &
                                                  * patch_wgt
                  cgrid%fmean_intercepted_aw(ipy) = cgrid%fmean_intercepted_aw (ipy)       &
                                                  + cpatch%fmean_intercepted_aw(ico)       &
                                                  * patch_wgt
                  cgrid%fmean_wshed_wg      (ipy) = cgrid%fmean_wshed_wg       (ipy)       &
                                                  + cpatch%fmean_wshed_wg      (ico)       &
                                                  * patch_wgt
               end do cohortloop
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !     Aggregate the patch-level variables.                                  !
               !---------------------------------------------------------------------------!
               cgrid%fmean_rh             (ipy) = cgrid%fmean_rh             (ipy)         &
                                                + csite%fmean_rh             (ipa)         &
                                                * patch_wgt
               cgrid%fmean_cwd_rh         (ipy) = cgrid%fmean_cwd_rh         (ipy)         &
                                                + csite%fmean_cwd_rh         (ipa)         &
                                                * patch_wgt
               cgrid%fmean_nep            (ipy) = cgrid%fmean_nep            (ipy)         &
                                                + csite%fmean_nep            (ipa)         &
                                                * patch_wgt
               cgrid%fmean_rk4step        (ipy) = cgrid%fmean_rk4step        (ipy)         &
                                                + csite%fmean_rk4step        (ipa)         &
                                                * patch_wgt
               cgrid%fmean_available_water(ipy) = cgrid%fmean_available_water(ipy)         &
                                                + csite%fmean_available_water(ipa)         &
                                                * patch_wgt
               cgrid%fmean_can_theiv      (ipy) = cgrid%fmean_can_theiv      (ipy)         &
                                                + csite%fmean_can_theiv      (ipa)         &
                                                * patch_wgt
               cgrid%fmean_can_theta      (ipy) = cgrid%fmean_can_theta      (ipy)         &
                                                + csite%fmean_can_theta      (ipa)         &
                                                * patch_wgt
               cgrid%fmean_can_vpdef      (ipy) = cgrid%fmean_can_vpdef      (ipy)         &
                                                + csite%fmean_can_vpdef      (ipa)         &
                                                * patch_wgt
               cgrid%fmean_can_shv        (ipy) = cgrid%fmean_can_shv        (ipy)         &
                                                + csite%fmean_can_shv        (ipa)         &
                                                * patch_wgt
               cgrid%fmean_can_co2        (ipy) = cgrid%fmean_can_co2        (ipy)         &
                                                + csite%fmean_can_co2        (ipa)         &
                                                * patch_wgt
               cgrid%fmean_can_prss       (ipy) = cgrid%fmean_can_prss       (ipy)         &
                                                + csite%fmean_can_prss       (ipa)         &
                                                * patch_wgt
               cgrid%fmean_gnd_temp       (ipy) = cgrid%fmean_gnd_temp       (ipy)         &
                                                + csite%fmean_gnd_temp       (ipa)         &
                                                * patch_wgt
               cgrid%fmean_gnd_shv        (ipy) = cgrid%fmean_gnd_shv        (ipy)         &
                                                + csite%fmean_gnd_shv        (ipa)         &
                                                * patch_wgt
               cgrid%fmean_can_ggnd       (ipy) = cgrid%fmean_can_ggnd       (ipy)         &
                                                + csite%fmean_can_ggnd       (ipa)         &
                                                * patch_wgt
               cgrid%fmean_sfcw_depth     (ipy) = cgrid%fmean_sfcw_depth     (ipy)         &
                                                + csite%fmean_sfcw_depth     (ipa)         &
                                                * patch_wgt
               !----- Temporarily convert pounding internal energy to J/m2. ---------------!
               cgrid%fmean_sfcw_energy    (ipy) = cgrid%fmean_sfcw_energy    (ipy)         &
                                                + csite%fmean_sfcw_energy    (ipa)         &
                                                * csite%fmean_sfcw_mass      (ipa)         &
                                                * patch_wgt
               !---------------------------------------------------------------------------!
               cgrid%fmean_sfcw_mass      (ipy) = cgrid%fmean_sfcw_mass      (ipy)         &
                                                + csite%fmean_sfcw_mass      (ipa)         &
                                                * patch_wgt
               cgrid%fmean_rshort_gnd     (ipy) = cgrid%fmean_rshort_gnd     (ipy)         &
                                                + csite%fmean_rshort_gnd     (ipa)         &
                                                * patch_wgt
               cgrid%fmean_par_gnd        (ipy) = cgrid%fmean_par_gnd        (ipy)         &
                                                + csite%fmean_par_gnd        (ipa)         &
                                                * patch_wgt
               cgrid%fmean_rlong_gnd      (ipy) = cgrid%fmean_rlong_gnd      (ipy)         &
                                                + csite%fmean_rlong_gnd      (ipa)         &
                                                * patch_wgt
               cgrid%fmean_rlongup        (ipy) = cgrid%fmean_rlongup        (ipy)         &
                                                + csite%fmean_rlongup        (ipa)         &
                                                * patch_wgt
               cgrid%fmean_parup          (ipy) = cgrid%fmean_parup          (ipy)         &
                                                + csite%fmean_parup          (ipa)         &
                                                * patch_wgt
               cgrid%fmean_nirup          (ipy) = cgrid%fmean_nirup          (ipy)         &
                                                + csite%fmean_nirup          (ipa)         &
                                                * patch_wgt
               cgrid%fmean_rshortup       (ipy) = cgrid%fmean_rshortup       (ipy)         &
                                                + csite%fmean_rshortup       (ipa)         &
                                                * patch_wgt
               cgrid%fmean_rnet           (ipy) = cgrid%fmean_rnet           (ipy)         &
                                                + csite%fmean_rnet           (ipa)         &
                                                * patch_wgt
               cgrid%fmean_albedo         (ipy) = cgrid%fmean_albedo         (ipy)         &
                                                + csite%fmean_albedo         (ipa)         &
                                                * patch_wgt
               cgrid%fmean_albedo_par     (ipy) = cgrid%fmean_albedo_par     (ipy)         &
                                                + csite%fmean_albedo_par     (ipa)         &
                                                * patch_wgt
               cgrid%fmean_albedo_nir     (ipy) = cgrid%fmean_albedo_nir     (ipy)         &
                                                + csite%fmean_albedo_nir     (ipa)         &
                                                * patch_wgt
               cgrid%fmean_rlong_albedo   (ipy) = cgrid%fmean_rlong_albedo   (ipy)         &
                                                + csite%fmean_rlong_albedo   (ipa)         &
                                                * patch_wgt
               cgrid%fmean_ustar          (ipy) = cgrid%fmean_ustar          (ipy)         &
                                                + csite%fmean_ustar          (ipa)         &
                                                * patch_wgt
               cgrid%fmean_tstar          (ipy) = cgrid%fmean_tstar          (ipy)         &
                                                + csite%fmean_tstar          (ipa)         &
                                                * patch_wgt
               cgrid%fmean_qstar          (ipy) = cgrid%fmean_qstar          (ipy)         &
                                                + csite%fmean_qstar          (ipa)         &
                                                * patch_wgt
               cgrid%fmean_cstar          (ipy) = cgrid%fmean_cstar          (ipy)         &
                                                + csite%fmean_cstar          (ipa)         &
                                                * patch_wgt
               cgrid%fmean_carbon_ac      (ipy) = cgrid%fmean_carbon_ac      (ipy)         &
                                                + csite%fmean_carbon_ac      (ipa)         &
                                                * patch_wgt
               cgrid%fmean_carbon_st      (ipy) = cgrid%fmean_carbon_st      (ipy)         &
                                                + csite%fmean_carbon_st      (ipa)         &
                                                * patch_wgt
               cgrid%fmean_vapor_gc       (ipy) = cgrid%fmean_vapor_gc       (ipy)         &
                                                + csite%fmean_vapor_gc       (ipa)         &
                                                * patch_wgt
               cgrid%fmean_vapor_ac       (ipy) = cgrid%fmean_vapor_ac       (ipy)         &
                                                + csite%fmean_vapor_ac       (ipa)         &
                                                * patch_wgt
               cgrid%fmean_throughfall    (ipy) = cgrid%fmean_throughfall    (ipy)         &
                                                + csite%fmean_throughfall    (ipa)         &
                                                * patch_wgt
               cgrid%fmean_runoff         (ipy) = cgrid%fmean_runoff         (ipy)         &
                                                + csite%fmean_runoff         (ipa)         &
                                                * patch_wgt
               cgrid%fmean_drainage       (ipy) = cgrid%fmean_drainage       (ipy)         &
                                                + csite%fmean_drainage       (ipa)         &
                                                * patch_wgt
               cgrid%fmean_sensible_gc    (ipy) = cgrid%fmean_sensible_gc    (ipy)         &
                                                + csite%fmean_sensible_gc    (ipa)         &
                                                * patch_wgt
               cgrid%fmean_sensible_ac    (ipy) = cgrid%fmean_sensible_ac    (ipy)         &
                                                + csite%fmean_sensible_ac    (ipa)         &
                                                * patch_wgt
               cgrid%fmean_qthroughfall   (ipy) = cgrid%fmean_qthroughfall   (ipy)         &
                                                + csite%fmean_qthroughfall   (ipa)         &
                                                * patch_wgt
               cgrid%fmean_qrunoff        (ipy) = cgrid%fmean_qrunoff        (ipy)         &
                                                + csite%fmean_qrunoff        (ipa)         &
                                                * patch_wgt
               cgrid%fmean_qdrainage      (ipy) = cgrid%fmean_qdrainage      (ipy)         &
                                                + csite%fmean_qdrainage      (ipa)         &
                                                * patch_wgt

               !----- Soil (extensive) properties. ----------------------------------------!
               do k=1,nzg
                  nsoil = cpoly%ntext_soil(k,isi)
                  cgrid%fmean_soil_energy (k,ipy) = cgrid%fmean_soil_energy (k,ipy)        &
                                                  + csite%fmean_soil_energy (k,ipa)        &
                                                  * patch_wgt
                  cgrid%fmean_soil_mstpot (k,ipy) = cgrid%fmean_soil_mstpot (k,ipy)        &
                                                  + csite%fmean_soil_mstpot (k,ipa)        &
                                                  * patch_wgt
                  cgrid%fmean_soil_water  (k,ipy) = cgrid%fmean_soil_water  (k,ipy)        &
                                                  + csite%fmean_soil_water  (k,ipa)        &
                                                  * patch_wgt
                  cgrid%fmean_smoist_gg   (k,ipy) = cgrid%fmean_smoist_gg   (k,ipy)        &
                                                  + csite%fmean_smoist_gg   (k,ipa)        &
                                                  * patch_wgt
                  cgrid%fmean_transloss   (k,ipy) = cgrid%fmean_transloss   (k,ipy)        &
                                                  + csite%fmean_transloss   (k,ipa)        &
                                                  * patch_wgt
                  cgrid%fmean_sensible_gg (k,ipy) = cgrid%fmean_sensible_gg (k,ipy)        &
                                                  + csite%fmean_sensible_gg (k,ipa)        &
                                                  * patch_wgt
                  cgrid_fmean_soil_hcap   (k)     = cgrid_fmean_soil_hcap   (k)            &
                                                  + soil(nsoil)%slcpd                      &
                                                  * patch_wgt

                  !----- Find the mean soil wetness. --------------------------------------!
                  cgrid%fmean_soil_wetness  (ipy) = cgrid%fmean_soil_wetness    (ipy)      &
                                                  + ( ( csite%fmean_soil_water(k,ipa)      &
                                                      - soil(nsoil)%soilwp) )              &
                                                  / ( soil(nsoil)%slmsts                   &
                                                    - soil(nsoil)%soilwp)                  &
                                                  * dslz(k) * dslzsum_i * patch_wgt
                  !------------------------------------------------------------------------!
               end do
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Aggregate the patch-level variables.                                     !
            !------------------------------------------------------------------------------!
            cgrid%fmean_atm_theiv      (ipy) = cgrid%fmean_atm_theiv      (ipy)            &
                                             + cpoly%fmean_atm_theiv      (isi)            &
                                             * site_wgt
            cgrid%fmean_atm_theta      (ipy) = cgrid%fmean_atm_theta      (ipy)            &
                                             + cpoly%fmean_atm_theta      (isi)            &
                                             * site_wgt
            cgrid%fmean_atm_temp       (ipy) = cgrid%fmean_atm_temp       (ipy)            &
                                             + cpoly%fmean_atm_temp       (isi)            &
                                             * site_wgt
            cgrid%fmean_atm_vpdef      (ipy) = cgrid%fmean_atm_vpdef      (ipy)            &
                                             + cpoly%fmean_atm_vpdef      (isi)            &
                                             * site_wgt
            cgrid%fmean_atm_shv        (ipy) = cgrid%fmean_atm_shv        (ipy)            &
                                             + cpoly%fmean_atm_shv        (isi)            &
                                             * site_wgt
            cgrid%fmean_atm_rshort     (ipy) = cgrid%fmean_atm_rshort     (ipy)            &
                                             + cpoly%fmean_atm_rshort     (isi)            &
                                             * site_wgt
            cgrid%fmean_atm_rshort_diff(ipy) = cgrid%fmean_atm_rshort_diff(ipy)            &
                                             + cpoly%fmean_atm_rshort_diff(isi)            &
                                             * site_wgt
            cgrid%fmean_atm_par        (ipy) = cgrid%fmean_atm_par        (ipy)            &
                                             + cpoly%fmean_atm_par        (isi)            &
                                             * site_wgt
            cgrid%fmean_atm_par_diff   (ipy) = cgrid%fmean_atm_par_diff   (ipy)            &
                                             + cpoly%fmean_atm_par_diff   (isi)            &
                                             * site_wgt
            cgrid%fmean_atm_rlong      (ipy) = cgrid%fmean_atm_rlong      (ipy)            &
                                             + cpoly%fmean_atm_rlong      (isi)            &
                                             * site_wgt
            cgrid%fmean_atm_vels       (ipy) = cgrid%fmean_atm_vels       (ipy)            &
                                             + cpoly%fmean_atm_vels       (isi)            &
                                             * site_wgt
            cgrid%fmean_atm_rhos       (ipy) = cgrid%fmean_atm_rhos       (ipy)            &
                                             + cpoly%fmean_atm_rhos       (isi)            &
                                             * site_wgt
            cgrid%fmean_atm_prss       (ipy) = cgrid%fmean_atm_prss       (ipy)            &
                                             + cpoly%fmean_atm_prss       (isi)            &
                                             * site_wgt
            cgrid%fmean_atm_co2        (ipy) = cgrid%fmean_atm_co2        (ipy)            &
                                             + cpoly%fmean_atm_co2        (isi)            &
                                             * site_wgt
            cgrid%fmean_pcpg           (ipy) = cgrid%fmean_pcpg           (ipy)            &
                                             + cpoly%fmean_pcpg           (isi)            &
                                             * site_wgt
            cgrid%fmean_qpcpg          (ipy) = cgrid%fmean_qpcpg          (ipy)            &
                                             + cpoly%fmean_qpcpg          (isi)            &
                                             * site_wgt
            cgrid%fmean_dpcpg          (ipy) = cgrid%fmean_dpcpg          (ipy)            &
                                             + cpoly%fmean_dpcpg          (isi)            &
                                             * site_wgt
         end do siteloop
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !      Find the derived properties for the air above canopy.                      !
         !---------------------------------------------------------------------------------!
         atm_exner                 = press2exner (cgrid%fmean_atm_prss(ipy))
         cgrid%fmean_atm_temp(ipy) = extheta2temp(atm_exner,cgrid%fmean_atm_theta(ipy))
         cgrid%fmean_atm_rhos(ipy) = idealdenssh ( cgrid%fmean_atm_prss  (ipy)             &
                                                 , cgrid%fmean_atm_temp  (ipy)             &
                                                 , cgrid%fmean_atm_shv   (ipy) )
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !      Find the derived properties for the canopy air space.                      !
         !---------------------------------------------------------------------------------!
         can_exner                 = press2exner (cgrid%fmean_can_prss(ipy))
         cgrid%fmean_can_temp(ipy) = extheta2temp(can_exner,cgrid%fmean_can_theta(ipy))
         cgrid%fmean_can_rhos(ipy) = idealdenssh ( cgrid%fmean_can_prss  (ipy)             &
                                                 , cgrid%fmean_can_temp  (ipy)             &
                                                 , cgrid%fmean_can_shv   (ipy) )
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !   If the patch had some temporary snow/pounding layer, convert the mean energy  !
         ! to J/kg, then find the mean temperature and liquid fraction.  Otherwise, set    !
         ! them to either zero or default values.                                          !
         !---------------------------------------------------------------------------------!
         if (cgrid%fmean_sfcw_mass(ipy) > tiny_sfcwater_mass) then
            cgrid%fmean_sfcw_energy(ipy) = cgrid%fmean_sfcw_energy(ipy)                    &
                                         / cgrid%fmean_sfcw_mass(ipy)
            call uint2tl(cgrid%fmean_sfcw_energy(ipy),cgrid%fmean_sfcw_temp(ipy)           &
                        ,cgrid%fmean_sfcw_fliq(ipy))
         else
            cgrid%fmean_sfcw_mass  (ipy)  = 0.
            cgrid%fmean_sfcw_depth (ipy)  = 0.
            cgrid%fmean_sfcw_energy(ipy)  = 0.
            cgrid%fmean_sfcw_temp  (ipy)  = cgrid%fmean_soil_temp(nzg,ipy)
            cgrid%fmean_sfcw_fliq  (ipy)  = cgrid%fmean_soil_fliq(nzg,ipy)
         end if
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Find the temperature and the fraction of liquid water.                      !
         !---------------------------------------------------------------------------------!
         do k=1,nzg
            call uextcm2tl( cgrid%fmean_soil_energy(k,ipy)                                 &
                          , cgrid%fmean_soil_water (k,ipy) * wdns                          &
                          , cgrid_fmean_soil_hcap  (k)                                     & 
                          , cgrid%fmean_soil_temp  (k,ipy)                                 &
                          , cgrid%fmean_soil_fliq  (k,ipy) )
         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the vegetation temperature and liquid fraction.                        !
         !---------------------------------------------------------------------------------!
         !----- Leaf. ---------------------------------------------------------------------!
         if (cgrid%fmean_leaf_hcap(ipy) > 0.) then
            call uextcm2tl( cgrid%fmean_leaf_energy(ipy), cgrid%fmean_leaf_water (ipy)     &
                          , cgrid%fmean_leaf_hcap  (ipy), cgrid%fmean_leaf_temp  (ipy)     &
                          , cgrid%fmean_leaf_fliq  (ipy) )
         else
            cgrid%fmean_leaf_temp (ipy) = cgrid%fmean_can_temp (ipy)
            if (cgrid%fmean_can_temp(ipy) > t00) then
               cgrid%fmean_leaf_fliq(ipy) = 1.0
            elseif (cgrid%fmean_can_temp(ipy) == t00) then
               cgrid%fmean_leaf_fliq(ipy) = 0.5
            else
               cgrid%fmean_leaf_fliq(ipy) = 0.0
            end if
         end if
         !----- Wood. ---------------------------------------------------------------------!
         if (cgrid%fmean_wood_hcap(ipy) > 0.) then
            call uextcm2tl( cgrid%fmean_wood_energy(ipy)                                   &
                          , cgrid%fmean_wood_water (ipy)                                   &
                          , cgrid%fmean_wood_hcap  (ipy)                                   &
                          , cgrid%fmean_wood_temp  (ipy)                                   &
                          , cgrid%fmean_wood_fliq  (ipy) )
         else
            cgrid%fmean_wood_temp(ipy) = cgrid%fmean_can_temp(ipy)
            if (cgrid%fmean_can_temp(ipy) > t00) then
               cgrid%fmean_wood_fliq(ipy) = 1.0
            elseif (cgrid%fmean_can_temp(ipy) == t00) then
               cgrid%fmean_wood_fliq(ipy) = 0.5
            else
               cgrid%fmean_wood_fliq(ipy) = 0.0
            end if
         end if
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Normalise the "intensive" properties.  The weight was either the LAI, WAI,   !
         ! or plant density.  In case none of the cohorts qualified to contribute, then we !
         ! assign either the canopy air space property, or a default number.               !
         !---------------------------------------------------------------------------------!
         if (poly_lai > 0.) then
            cgrid%fmean_leaf_vpdef(ipy) = cgrid%fmean_leaf_vpdef (ipy) / poly_lai
            cgrid%fmean_leaf_gsw  (ipy) = cgrid%fmean_leaf_gsw   (ipy) / poly_lai
            cgrid%fmean_leaf_gbw  (ipy) = cgrid%fmean_leaf_gbw   (ipy) / poly_lai
            cgrid%fmean_fs_open   (ipy) = cgrid%fmean_fs_open    (ipy) / poly_lai
            cgrid%fmean_fsw       (ipy) = cgrid%fmean_fsw        (ipy) / poly_lai
            cgrid%fmean_fsn       (ipy) = cgrid%fmean_fsn        (ipy) / poly_lai
         else
            cgrid%fmean_leaf_vpdef(ipy) = cgrid%fmean_can_vpdef  (ipy)
            cgrid%fmean_leaf_gsw  (ipy) = 0.0
            cgrid%fmean_leaf_gbw  (ipy) = 0.0
            cgrid%fmean_fs_open   (ipy) = 0.5
            cgrid%fmean_fsw       (ipy) = 0.5
            cgrid%fmean_fsn       (ipy) = 0.5
         end if
         if (poly_wai > 0.) then
            cgrid%fmean_wood_gbw  (ipy) = cgrid%fmean_wood_gbw   (ipy) / poly_wai
         else
            cgrid%fmean_wood_gbw  (ipy) = 0.0
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Compute the average amongst all surfaces (soil, temporary surface water, and !
         ! vegetation, the last two only if they really exist).  All energy terms are      !
         ! converted to J/m2, all water terms to kg/m2, and the heat capacities of every-  !
         ! thing that is not water is in J/m2/K.                                           !
         !---------------------------------------------------------------------------------!
         skin_energy = cgrid%fmean_leaf_energy(ipy) + cgrid%fmean_wood_energy(ipy)         &
                     + cgrid%fmean_sfcw_energy(ipy) * cgrid%fmean_sfcw_mass  (ipy)         &
                     + cgrid%fmean_soil_energy(nzg,ipy) * dslz(nzg)
         skin_water  = cgrid%fmean_leaf_water(ipy)  + cgrid%fmean_wood_water(ipy)          &
                     + cgrid%fmean_sfcw_mass (ipy)                                         &
                     + cgrid%fmean_soil_water(nzg,ipy) * dslz(nzg) * wdns
         skin_hcap   = cgrid%fmean_leaf_hcap(ipy) + cgrid%fmean_wood_hcap(ipy)             &
                     + cgrid_fmean_soil_hcap(nzg) * dslz(nzg)
         call uextcm2tl(skin_energy,skin_water,skin_hcap                                   &
                       ,cgrid%fmean_skin_temp(ipy),skin_fliq)
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!
      return
   end subroutine aggregate_polygon_fmean
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine increments the time averaged site met-forcing variables.  The     !
   ! polygon-level averages are found after the site-level are normalised.                 !
   !---------------------------------------------------------------------------------------!
   subroutine integrate_ed_fmean_met_vars(cgrid)
      use ed_state_vars  , only : edtype          & ! structure
                                , polygontype     ! ! structure
      use met_driver_coms, only : met_driv_state  ! ! structure
      use ed_misc_coms   , only : dtlsm           & ! intent(in)
                                , frqsum          ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)      , target  :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)   , pointer :: cpoly
      type(met_driv_state), pointer :: cmet
      integer                       :: ipy
      integer                       :: isi
      !----- Locally saved variables. -----------------------------------------------------!
      real              , save      :: dtlsm_o_frqsum = 1.e34
      logical           , save      :: first_time     = .true.
      !------------------------------------------------------------------------------------!


      !----- Assign the constant scaling factor. ------------------------------------------!
      if (first_time) then
         first_time     = .false.
         dtlsm_o_frqsum = dtlsm / frqsum
      end if
      !------------------------------------------------------------------------------------!

      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            cmet  => cpoly%met(isi)

            !----- The site-level averages. -----------------------------------------------!
            cpoly%fmean_atm_theiv      (isi) = cpoly%fmean_atm_theiv      (isi)            &
                                             + cmet%atm_theiv                              &
                                             * dtlsm_o_frqsum
            cpoly%fmean_atm_theta      (isi) = cpoly%fmean_atm_theta      (isi)            &
                                             + cmet%atm_theta                              &
                                             * dtlsm_o_frqsum
            cpoly%fmean_atm_vpdef      (isi) = cpoly%fmean_atm_vpdef      (isi)            &
                                             + cmet%atm_vpdef                              &
                                             * dtlsm_o_frqsum
            cpoly%fmean_atm_shv        (isi) = cpoly%fmean_atm_shv        (isi)            &
                                             + cmet%atm_shv                                &
                                             * dtlsm_o_frqsum
            cpoly%fmean_atm_rshort     (isi) = cpoly%fmean_atm_rshort     (isi)            &
                                             + cmet%rshort                                 &
                                             * dtlsm_o_frqsum
            cpoly%fmean_atm_rshort_diff(isi) = cpoly%fmean_atm_rshort_diff(isi)            &
                                             + cmet%rshort_diffuse                         &
                                             * dtlsm_o_frqsum
            cpoly%fmean_atm_par        (isi) = cpoly%fmean_atm_par        (isi)            &
                                             + ( cmet%par_beam + cmet%par_diffuse )        &
                                             * dtlsm_o_frqsum
            cpoly%fmean_atm_par_diff   (isi) = cpoly%fmean_atm_par_diff   (isi)            &
                                             + cmet%par_diffuse                            &
                                             * dtlsm_o_frqsum
            cpoly%fmean_atm_rlong      (isi) = cpoly%fmean_atm_rlong      (isi)            &
                                             + cmet%rlong                                  &
                                             * dtlsm_o_frqsum
            cpoly%fmean_atm_vels       (isi) = cpoly%fmean_atm_vels       (isi)            &
                                             + cmet%vels                                   &
                                             * dtlsm_o_frqsum
            cpoly%fmean_atm_prss       (isi) = cpoly%fmean_atm_prss       (isi)            &
                                             + cmet%prss                                   &
                                             * dtlsm_o_frqsum
            cpoly%fmean_atm_co2        (isi) = cpoly%fmean_atm_co2        (isi)            &
                                             + cmet%atm_co2                                &
                                             * dtlsm_o_frqsum
            cpoly%fmean_pcpg           (isi) = cpoly%fmean_pcpg           (isi)            &
                                             + cmet%pcpg                                   &
                                             * dtlsm_o_frqsum
            cpoly%fmean_qpcpg          (isi) = cpoly%fmean_qpcpg          (isi)            &
                                             + cmet%qpcpg                                  &
                                             * dtlsm_o_frqsum
            cpoly%fmean_dpcpg          (isi) = cpoly%fmean_dpcpg          (isi)            &
                                             + cmet%dpcpg                                  &
                                             * dtlsm_o_frqsum
            !------------------------------------------------------------------------------!
         end do siteloop
      end do polyloop
      return
   end subroutine integrate_ed_fmean_met_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !       The following sub-routine scales several variables that are integrated during   !
   ! one output step (frqsum) to actual rates, and find derived properties.                !
   !---------------------------------------------------------------------------------------!
   subroutine normalize_ed_fmean_vars(cgrid)
      use grid_coms    , only : nzg                ! ! intent(in)
      use ed_misc_coms , only : dtlsm              & ! intent(in)
                              , frqsum             & ! intent(in)
                              , radfrq             & ! intent(in)
                              , current_time       ! ! intent(in)
      use ed_state_vars, only : edtype             & ! structure
                              , polygontype        & ! structure
                              , sitetype           & ! structure
                              , patchtype          ! ! structure
      use therm_lib    , only : uextcm2tl          & ! subroutine
                              , uint2tl            & ! subroutine
                              , idealdenssh        & ! function
                              , press2exner        & ! function
                              , extheta2temp       ! ! function
      use consts_coms  , only : t00                & ! intent(in)
                              , wdns               ! ! intent(in)
      use soil_coms    , only : tiny_sfcwater_mass & ! intent(in)
                              , isoilbc            & ! intent(in)
                              , soil               & ! intent(in)
                              , dslz               & ! intent(in)
                              , matric_potential   ! ! function
      implicit none
      !----- Arguments.  ------------------------------------------------------------------!
      type(edtype)      , target     :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype) , pointer    :: cpoly
      type(sitetype)    , pointer    :: csite
      type(patchtype)   , pointer    :: cpatch
      integer                        :: ipy
      integer                        :: isi
      integer                        :: ipa
      integer                        :: ico
      integer                        :: nsoil
      real                           :: dtlsm_o_frqsum
      real                           :: radfrq_o_frqsum
      real                           :: frqsumi
      real                           :: pss_npp
      real                           :: pss_lai
      real                           :: atm_exner
      real                           :: can_exner
      integer                        :: k
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Find some useful conversion factors.                                           !
      ! 1. FRQSUMI         -- inverse of the elapsed time between two analyses (or one     !
      !                       day).  This should be used by variables that are fluxes and  !
      !                       are solved by RK4, they are holding the integral over the    !
      !                       past frqsum seconds.                                         !
      ! 2. DTLSM_O_FRQSUM  -- inverse of the number of the main time steps (DTLSM) since   !
      !                       previous analysis.  Only photosynthesis- and decomposition-  !
      !                       related variables, or STATE VARIABLES should use this        !
      !                       factor.  Do not use this for energy and water fluxes, CO2    !
      !                       eddy flux, and CO2 storage.                                  !
      ! 3. RADFRQ_O_FRQSUM -- inverse of the number of radiation time steps since the      !
      !                       previous analysis.  Only radiation-related variables should  !
      !                       use this factor.                                             !
      !------------------------------------------------------------------------------------!
      frqsumi         =    1.0 / frqsum
      dtlsm_o_frqsum  =  dtlsm * frqsumi
      radfrq_o_frqsum = radfrq * frqsumi
      !------------------------------------------------------------------------------------!




      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !------------------------------------------------------------------------------!
            !      Now we find the derived properties for the air above canopy.            !
            !------------------------------------------------------------------------------!
            atm_exner                 = press2exner (cpoly%fmean_atm_prss(isi))
            cpoly%fmean_atm_temp(isi) = extheta2temp(atm_exner,cpoly%fmean_atm_theta(isi))
            cpoly%fmean_atm_rhos(isi) = idealdenssh ( cpoly%fmean_atm_prss  (isi)          &
                                                    , cpoly%fmean_atm_temp  (isi)          &
                                                    , cpoly%fmean_atm_shv   (isi) )
            !------------------------------------------------------------------------------!



            patchloop: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !---------------------------------------------------------------------------!
               !      Reset the patch-level GPP and plant respiration, used to find NEP.   !
               !---------------------------------------------------------------------------!
               pss_npp    = 0.0
               pss_lai    = 0.0
               !---------------------------------------------------------------------------!





               !---------------------------------------------------------------------------!
               !     The following variables are fluxes that cam from the RK4 and there-   !
               ! fore hold the integral over time, divide them by the total time to obtain !
               ! the mean fluxes.                                                          !
               !---------------------------------------------------------------------------!
               csite%fmean_rk4step       (ipa) = csite%fmean_rk4step       (ipa) * frqsumi
               csite%fmean_ustar         (ipa) = csite%fmean_ustar         (ipa) * frqsumi
               csite%fmean_tstar         (ipa) = csite%fmean_tstar         (ipa) * frqsumi
               csite%fmean_qstar         (ipa) = csite%fmean_qstar         (ipa) * frqsumi
               csite%fmean_cstar         (ipa) = csite%fmean_cstar         (ipa) * frqsumi
               csite%fmean_carbon_ac     (ipa) = csite%fmean_carbon_ac     (ipa) * frqsumi
               csite%fmean_carbon_st     (ipa) = csite%fmean_carbon_st     (ipa) * frqsumi
               csite%fmean_vapor_gc      (ipa) = csite%fmean_vapor_gc      (ipa) * frqsumi
               csite%fmean_vapor_ac      (ipa) = csite%fmean_vapor_ac      (ipa) * frqsumi
               csite%fmean_throughfall   (ipa) = csite%fmean_throughfall   (ipa) * frqsumi
               csite%fmean_runoff        (ipa) = csite%fmean_runoff        (ipa) * frqsumi
               csite%fmean_drainage      (ipa) = csite%fmean_drainage      (ipa) * frqsumi
               csite%fmean_sensible_gc   (ipa) = csite%fmean_sensible_gc   (ipa) * frqsumi
               csite%fmean_sensible_ac   (ipa) = csite%fmean_sensible_ac   (ipa) * frqsumi
               csite%fmean_qthroughfall  (ipa) = csite%fmean_qthroughfall  (ipa) * frqsumi
               csite%fmean_qrunoff       (ipa) = csite%fmean_qrunoff       (ipa) * frqsumi
               csite%fmean_qdrainage     (ipa) = csite%fmean_qdrainage     (ipa) * frqsumi
               !------ Soil flux. ---------------------------------------------------------!
               do k=cpoly%lsl(isi),nzg
                  csite%fmean_sensible_gg(k,ipa) = csite%fmean_sensible_gg   (k,ipa)       &
                                                 * frqsumi
                  csite%fmean_smoist_gg  (k,ipa) = csite%fmean_smoist_gg     (k,ipa)       &
                                                 * frqsumi
                  csite%fmean_transloss  (k,ipa) = csite%fmean_transloss     (k,ipa)       &
                                                 * frqsumi
               end do
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !---------------------------------------------------------------------------!
               !     Most state variables are already normalised.  All that we need to do  !
               ! is to find the derived properties.                                        !
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !     Soil matric potential, temperature, and liquid water.                 !
               !---------------------------------------------------------------------------!
               do k=1,nzg
                  nsoil = cpoly%ntext_soil(k,isi)
                  call uextcm2tl( csite%fmean_soil_energy(k,ipa)                           &
                                , csite%fmean_soil_water (k,ipa) * wdns                    &
                                , soil(nsoil)%slcpd                                        &
                                , csite%fmean_soil_temp  (k,ipa)                           &
                                , csite%fmean_soil_fliq  (k,ipa))

                  csite%fmean_soil_mstpot (k,ipa) =                                        &
                                     matric_potential(nsoil,csite%fmean_soil_water (k,ipa))
                  
               end do
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !      Now we find the derived properties for the canopy air space.         !
               !---------------------------------------------------------------------------!
               can_exner                 = press2exner ( csite%fmean_can_prss  (ipa) )
               csite%fmean_can_temp(ipa) = extheta2temp( can_exner                         &
                                                       , csite%fmean_can_theta (ipa) )
               csite%fmean_can_rhos(ipa) = idealdenssh ( csite%fmean_can_prss  (ipa)       &
                                                       , csite%fmean_can_temp  (ipa)       &
                                                       , csite%fmean_can_shv   (ipa)       )
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !   If the patch had some temporary snow/pounding layer, convert the mean   !
               ! energy to J/kg, then find the mean temperature and liquid fraction.       !
               ! Otherwise, set them to either zero or default values.                     !
               !---------------------------------------------------------------------------!
               if (csite%fmean_sfcw_mass(ipa) > tiny_sfcwater_mass) then
                  csite%fmean_sfcw_energy(ipa) = csite%fmean_sfcw_energy(ipa)              &
                                               / csite%fmean_sfcw_mass(ipa)
                  call uint2tl(csite%fmean_sfcw_energy(ipa),csite%fmean_sfcw_temp(ipa)     &
                              ,csite%fmean_sfcw_fliq(ipa))
               else
                  csite%fmean_sfcw_mass  (ipa)  = 0.
                  csite%fmean_sfcw_depth (ipa)  = 0.
                  csite%fmean_sfcw_energy(ipa)  = 0.
                  csite%fmean_sfcw_temp  (ipa)  = csite%fmean_soil_temp(nzg,ipa)
                  csite%fmean_sfcw_fliq  (ipa)  = csite%fmean_soil_fliq(nzg,ipa)
               end if
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !      Loop over the cohorts and find the mean for derived properties.      !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts



                  !------------------------------------------------------------------------!
                  !    Energy and water fluxes were integrated over the past frqsum        !
                  ! interval.   Use frqsumi to normalise them.  Energy fluxes will become  !
                  ! W/m², and water fluxes will become kg/m²/s.                            !
                  !------------------------------------------------------------------------!
                  cpatch%fmean_sensible_lc   (ico) = cpatch%fmean_sensible_lc   (ico)      &
                                                   * frqsumi
                  cpatch%fmean_vapor_lc      (ico) = cpatch%fmean_vapor_lc      (ico)      &
                                                   * frqsumi
                  cpatch%fmean_transp        (ico) = cpatch%fmean_transp        (ico)      &
                                                   * frqsumi
                  cpatch%fmean_intercepted_al(ico) = cpatch%fmean_intercepted_al(ico)      &
                                                   * frqsumi
                  cpatch%fmean_wshed_lg      (ico) = cpatch%fmean_wshed_lg      (ico)      &
                                                   * frqsumi
                  cpatch%fmean_sensible_wc   (ico) = cpatch%fmean_sensible_wc   (ico)      &
                                                   * frqsumi
                  cpatch%fmean_vapor_wc      (ico) = cpatch%fmean_vapor_wc      (ico)      &
                                                   * frqsumi
                  cpatch%fmean_intercepted_aw(ico) = cpatch%fmean_intercepted_aw(ico)      &
                                                   * frqsumi
                  cpatch%fmean_wshed_wg      (ico) = cpatch%fmean_wshed_wg      (ico)      &
                                                   * frqsumi
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Find the vegetation temperature and liquid fraction.               !
                  !------------------------------------------------------------------------!
                  !----- Leaf. ------------------------------------------------------------!
                  if (cpatch%fmean_leaf_hcap(ico) > 0.) then
                     call uextcm2tl( cpatch%fmean_leaf_energy(ico)                         &
                                   , cpatch%fmean_leaf_water (ico)                         &
                                   , cpatch%fmean_leaf_hcap  (ico)                         &
                                   , cpatch%fmean_leaf_temp  (ico)                         &
                                   , cpatch%fmean_leaf_fliq  (ico) )
                  else
                     cpatch%fmean_leaf_vpdef(ico) = csite%fmean_can_vpdef(ipa)
                     cpatch%fmean_leaf_temp (ico) = csite%fmean_can_temp (ipa)
                     if (csite%fmean_can_temp(ipa) > t00) then
                        cpatch%fmean_leaf_fliq(ico) = 1.0
                     elseif (csite%fmean_can_temp(ipa) == t00) then
                        cpatch%fmean_leaf_fliq(ico) = 0.5
                     else
                        cpatch%fmean_leaf_fliq(ico) = 0.0
                     end if
                  end if
                  !----- Wood. ------------------------------------------------------------!
                  if (cpatch%fmean_wood_hcap(ico) > 0.) then
                     call uextcm2tl( cpatch%fmean_wood_energy(ico)                         &
                                   , cpatch%fmean_wood_water (ico)                         &
                                   , cpatch%fmean_wood_hcap  (ico)                         &
                                   , cpatch%fmean_wood_temp  (ico)                         &
                                   , cpatch%fmean_wood_fliq  (ico) )
                  else
                     cpatch%fmean_wood_temp(ico) = csite%fmean_can_temp(ipa)
                     if (csite%fmean_can_temp(ipa) > t00) then
                        cpatch%fmean_wood_fliq(ico) = 1.0
                     elseif (csite%fmean_can_temp(ipa) == t00) then
                        cpatch%fmean_wood_fliq(ico) = 0.5
                     else
                        cpatch%fmean_wood_fliq(ico) = 0.0
                     end if
                  end if
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !      Integrate the total plant respiration and net primary             !
                  ! productivity.                                                          !
                  !------------------------------------------------------------------------!
                  cpatch%fmean_plresp(ico) = cpatch%fmean_leaf_resp   (ico)                &
                                           + cpatch%fmean_root_resp   (ico)                &
                                           + cpatch%fmean_storage_resp(ico)                &
                                           + cpatch%fmean_growth_resp (ico)                &
                                           + cpatch%fmean_vleaf_resp  (ico)
                  cpatch%fmean_npp   (ico) = cpatch%fmean_gpp         (ico)                &
                                           - cpatch%fmean_plresp      (ico)
                  !------------------------------------------------------------------------!



                  !----- Add LAI and extensive NPP to compute NEP. ------------------------!
                  pss_lai    = pss_lai    + cpatch%lai       (ico)
                  pss_npp    = pss_npp    + cpatch%fmean_npp (ico) * cpatch%nplant(ico)
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Net Ecosystem Productivity found by combining the cohort-level terms  !
               ! (gross primary productivity and total plant respiration), with hetero-    !
               ! trophic respiration, a patch-level variable.                              !
               !---------------------------------------------------------------------------!
               csite%fmean_nep(ipa) = pss_npp - csite%fmean_rh(ipa)
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !      Budget variables.  They contain integral values, so we must divide   !
               ! by the elapsed time to get them in flux units.                            !
               !---------------------------------------------------------------------------!
               csite%co2budget_gpp        (ipa) = csite%co2budget_gpp        (ipa) * frqsumi
               csite%co2budget_plresp     (ipa) = csite%co2budget_plresp     (ipa) * frqsumi
               csite%co2budget_rh         (ipa) = csite%co2budget_rh         (ipa) * frqsumi
               csite%co2budget_loss2atm   (ipa) = csite%co2budget_loss2atm   (ipa) * frqsumi
               csite%co2budget_denseffect (ipa) = csite%co2budget_denseffect (ipa) * frqsumi
               csite%co2budget_residual   (ipa) = csite%co2budget_residual   (ipa) * frqsumi
               csite%ebudget_precipgain   (ipa) = csite%ebudget_precipgain   (ipa) * frqsumi
               csite%ebudget_netrad       (ipa) = csite%ebudget_netrad       (ipa) * frqsumi
               csite%ebudget_denseffect   (ipa) = csite%ebudget_denseffect   (ipa) * frqsumi
               csite%ebudget_prsseffect   (ipa) = csite%ebudget_prsseffect   (ipa) * frqsumi
               csite%ebudget_loss2atm     (ipa) = csite%ebudget_loss2atm     (ipa) * frqsumi
               csite%ebudget_loss2drainage(ipa) = csite%ebudget_loss2drainage(ipa) * frqsumi
               csite%ebudget_loss2runoff  (ipa) = csite%ebudget_loss2runoff  (ipa) * frqsumi
               csite%ebudget_residual     (ipa) = csite%ebudget_residual     (ipa) * frqsumi
               csite%wbudget_precipgain   (ipa) = csite%wbudget_precipgain   (ipa) * frqsumi
               csite%wbudget_loss2atm     (ipa) = csite%wbudget_loss2atm     (ipa) * frqsumi
               csite%wbudget_loss2drainage(ipa) = csite%wbudget_loss2drainage(ipa) * frqsumi
               csite%wbudget_loss2runoff  (ipa) = csite%wbudget_loss2runoff  (ipa) * frqsumi
               csite%wbudget_denseffect   (ipa) = csite%wbudget_denseffect   (ipa) * frqsumi
               csite%wbudget_residual     (ipa) = csite%wbudget_residual     (ipa) * frqsumi
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!
      return
   end subroutine normalize_ed_fmean_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine zero_ed_fmean_vars(cgrid)

      use ed_state_vars, only : edtype      & ! structure
                              , polygontype & ! structure
                              , sitetype    & ! structure
                              , patchtype   ! ! structure
      
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target  :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer :: cpoly
      type(sitetype)   , pointer :: csite
      type(patchtype)  , pointer :: cpatch
      integer                    :: ipy
      integer                    :: isi
      integer                    :: ipa
      integer                    :: ico
      !------------------------------------------------------------------------------------!

      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         cgrid%fmean_gpp             (  ipy) = 0.0
         cgrid%fmean_npp             (  ipy) = 0.0
         cgrid%fmean_leaf_resp       (  ipy) = 0.0
         cgrid%fmean_root_resp       (  ipy) = 0.0
         cgrid%fmean_growth_resp     (  ipy) = 0.0
         cgrid%fmean_storage_resp    (  ipy) = 0.0
         cgrid%fmean_vleaf_resp      (  ipy) = 0.0
         cgrid%fmean_plresp          (  ipy) = 0.0
         cgrid%fmean_leaf_energy     (  ipy) = 0.0
         cgrid%fmean_leaf_water      (  ipy) = 0.0
         cgrid%fmean_leaf_hcap       (  ipy) = 0.0
         cgrid%fmean_leaf_vpdef      (  ipy) = 0.0
         cgrid%fmean_leaf_temp       (  ipy) = 0.0
         cgrid%fmean_leaf_fliq       (  ipy) = 0.0
         cgrid%fmean_leaf_gsw        (  ipy) = 0.0
         cgrid%fmean_leaf_gbw        (  ipy) = 0.0
         cgrid%fmean_wood_energy     (  ipy) = 0.0
         cgrid%fmean_wood_water      (  ipy) = 0.0
         cgrid%fmean_wood_hcap       (  ipy) = 0.0
         cgrid%fmean_wood_temp       (  ipy) = 0.0
         cgrid%fmean_wood_fliq       (  ipy) = 0.0
         cgrid%fmean_wood_gbw        (  ipy) = 0.0
         cgrid%fmean_fs_open         (  ipy) = 0.0
         cgrid%fmean_fsw             (  ipy) = 0.0
         cgrid%fmean_fsn             (  ipy) = 0.0
         cgrid%fmean_a_light         (  ipy) = 0.0
         cgrid%fmean_a_rubp          (  ipy) = 0.0
         cgrid%fmean_a_co2           (  ipy) = 0.0
         cgrid%fmean_psi_open        (  ipy) = 0.0
         cgrid%fmean_psi_closed      (  ipy) = 0.0
         cgrid%fmean_water_supply    (  ipy) = 0.0
         cgrid%fmean_par_l           (  ipy) = 0.0
         cgrid%fmean_par_l_beam      (  ipy) = 0.0
         cgrid%fmean_par_l_diff      (  ipy) = 0.0
         cgrid%fmean_rshort_l        (  ipy) = 0.0
         cgrid%fmean_rlong_l         (  ipy) = 0.0
         cgrid%fmean_sensible_lc     (  ipy) = 0.0
         cgrid%fmean_vapor_lc        (  ipy) = 0.0
         cgrid%fmean_transp          (  ipy) = 0.0
         cgrid%fmean_intercepted_al  (  ipy) = 0.0
         cgrid%fmean_wshed_lg        (  ipy) = 0.0
         cgrid%fmean_rshort_w        (  ipy) = 0.0
         cgrid%fmean_rlong_w         (  ipy) = 0.0
         cgrid%fmean_sensible_wc     (  ipy) = 0.0
         cgrid%fmean_vapor_wc        (  ipy) = 0.0
         cgrid%fmean_intercepted_aw  (  ipy) = 0.0
         cgrid%fmean_wshed_wg        (  ipy) = 0.0
         cgrid%fmean_rh              (  ipy) = 0.0
         cgrid%fmean_cwd_rh          (  ipy) = 0.0
         cgrid%fmean_nep             (  ipy) = 0.0
         cgrid%fmean_rk4step         (  ipy) = 0.0
         cgrid%fmean_available_water (  ipy) = 0.0
         cgrid%fmean_can_theiv       (  ipy) = 0.0
         cgrid%fmean_can_theta       (  ipy) = 0.0
         cgrid%fmean_can_vpdef       (  ipy) = 0.0
         cgrid%fmean_can_temp        (  ipy) = 0.0
         cgrid%fmean_can_shv         (  ipy) = 0.0
         cgrid%fmean_can_co2         (  ipy) = 0.0
         cgrid%fmean_can_rhos        (  ipy) = 0.0
         cgrid%fmean_can_prss        (  ipy) = 0.0
         cgrid%fmean_gnd_temp        (  ipy) = 0.0
         cgrid%fmean_gnd_shv         (  ipy) = 0.0
         cgrid%fmean_can_ggnd        (  ipy) = 0.0
         cgrid%fmean_sfcw_depth      (  ipy) = 0.0
         cgrid%fmean_sfcw_energy     (  ipy) = 0.0
         cgrid%fmean_sfcw_mass       (  ipy) = 0.0
         cgrid%fmean_sfcw_temp       (  ipy) = 0.0
         cgrid%fmean_sfcw_fliq       (  ipy) = 0.0
         cgrid%fmean_soil_energy     (:,ipy) = 0.0
         cgrid%fmean_soil_mstpot     (:,ipy) = 0.0
         cgrid%fmean_soil_water      (:,ipy) = 0.0
         cgrid%fmean_soil_temp       (:,ipy) = 0.0
         cgrid%fmean_soil_fliq       (:,ipy) = 0.0
         cgrid%fmean_rshort_gnd      (  ipy) = 0.0
         cgrid%fmean_par_gnd         (  ipy) = 0.0
         cgrid%fmean_rlong_gnd       (  ipy) = 0.0
         cgrid%fmean_rlongup         (  ipy) = 0.0
         cgrid%fmean_parup           (  ipy) = 0.0
         cgrid%fmean_nirup           (  ipy) = 0.0
         cgrid%fmean_rshortup        (  ipy) = 0.0
         cgrid%fmean_rnet            (  ipy) = 0.0
         cgrid%fmean_albedo          (  ipy) = 0.0
         cgrid%fmean_albedo_par      (  ipy) = 0.0
         cgrid%fmean_albedo_nir      (  ipy) = 0.0
         cgrid%fmean_rlong_albedo    (  ipy) = 0.0
         cgrid%fmean_ustar           (  ipy) = 0.0
         cgrid%fmean_tstar           (  ipy) = 0.0
         cgrid%fmean_qstar           (  ipy) = 0.0
         cgrid%fmean_cstar           (  ipy) = 0.0
         cgrid%fmean_carbon_ac       (  ipy) = 0.0
         cgrid%fmean_carbon_st       (  ipy) = 0.0
         cgrid%fmean_vapor_gc        (  ipy) = 0.0
         cgrid%fmean_vapor_ac        (  ipy) = 0.0
         cgrid%fmean_smoist_gg       (:,ipy) = 0.0
         cgrid%fmean_throughfall     (  ipy) = 0.0
         cgrid%fmean_transloss       (:,ipy) = 0.0
         cgrid%fmean_runoff          (  ipy) = 0.0
         cgrid%fmean_drainage        (  ipy) = 0.0
         cgrid%fmean_sensible_gc     (  ipy) = 0.0
         cgrid%fmean_sensible_ac     (  ipy) = 0.0
         cgrid%fmean_sensible_gg     (:,ipy) = 0.0
         cgrid%fmean_qthroughfall    (  ipy) = 0.0
         cgrid%fmean_qrunoff         (  ipy) = 0.0
         cgrid%fmean_qdrainage       (  ipy) = 0.0
         cgrid%fmean_atm_theiv       (  ipy) = 0.0
         cgrid%fmean_atm_theta       (  ipy) = 0.0
         cgrid%fmean_atm_temp        (  ipy) = 0.0
         cgrid%fmean_atm_vpdef       (  ipy) = 0.0
         cgrid%fmean_atm_shv         (  ipy) = 0.0
         cgrid%fmean_atm_rshort      (  ipy) = 0.0
         cgrid%fmean_atm_rshort_diff (  ipy) = 0.0
         cgrid%fmean_atm_par         (  ipy) = 0.0
         cgrid%fmean_atm_par_diff    (  ipy) = 0.0
         cgrid%fmean_atm_rlong       (  ipy) = 0.0
         cgrid%fmean_atm_vels        (  ipy) = 0.0
         cgrid%fmean_atm_rhos        (  ipy) = 0.0
         cgrid%fmean_atm_prss        (  ipy) = 0.0
         cgrid%fmean_atm_co2         (  ipy) = 0.0
         cgrid%fmean_pcpg            (  ipy) = 0.0
         cgrid%fmean_qpcpg           (  ipy) = 0.0
         cgrid%fmean_dpcpg           (  ipy) = 0.0
         cgrid%fmean_soil_wetness    (  ipy) = 0.0
         cgrid%fmean_skin_temp       (  ipy) = 0.0

         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            cpoly%fmean_atm_theiv      (isi) = 0.0
            cpoly%fmean_atm_theta      (isi) = 0.0
            cpoly%fmean_atm_temp       (isi) = 0.0
            cpoly%fmean_atm_vpdef      (isi) = 0.0
            cpoly%fmean_atm_shv        (isi) = 0.0
            cpoly%fmean_atm_rshort     (isi) = 0.0
            cpoly%fmean_atm_rshort_diff(isi) = 0.0
            cpoly%fmean_atm_par        (isi) = 0.0
            cpoly%fmean_atm_par_diff   (isi) = 0.0
            cpoly%fmean_atm_rlong      (isi) = 0.0
            cpoly%fmean_atm_vels       (isi) = 0.0
            cpoly%fmean_atm_rhos       (isi) = 0.0
            cpoly%fmean_atm_prss       (isi) = 0.0
            cpoly%fmean_atm_co2        (isi) = 0.0
            cpoly%fmean_pcpg           (isi) = 0.0
            cpoly%fmean_qpcpg          (isi) = 0.0
            cpoly%fmean_dpcpg          (isi) = 0.0

            patchloop: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Budget variables. ---------------------------------------------------!
               csite%co2budget_gpp                 (ipa) = 0.0
               csite%co2budget_rh                  (ipa) = 0.0
               csite%co2budget_plresp              (ipa) = 0.0
               csite%co2budget_residual            (ipa) = 0.0
               csite%co2budget_loss2atm            (ipa) = 0.0
               csite%co2budget_denseffect          (ipa) = 0.0
               csite%wbudget_precipgain            (ipa) = 0.0
               csite%wbudget_loss2atm              (ipa) = 0.0
               csite%wbudget_loss2runoff           (ipa) = 0.0
               csite%wbudget_loss2drainage         (ipa) = 0.0
               csite%wbudget_denseffect            (ipa) = 0.0
               csite%wbudget_residual              (ipa) = 0.0
               csite%ebudget_precipgain            (ipa) = 0.0
               csite%ebudget_netrad                (ipa) = 0.0
               csite%ebudget_loss2atm              (ipa) = 0.0
               csite%ebudget_loss2runoff           (ipa) = 0.0
               csite%ebudget_loss2drainage         (ipa) = 0.0
               csite%ebudget_denseffect            (ipa) = 0.0
               csite%ebudget_prsseffect            (ipa) = 0.0
               csite%ebudget_residual              (ipa) = 0.0
               !---------------------------------------------------------------------------!




               !----- Fast average variables. ---------------------------------------------!
               csite%fmean_rh             (  ipa) = 0.0
               csite%fmean_cwd_rh         (  ipa) = 0.0
               csite%fmean_nep            (  ipa) = 0.0
               csite%fmean_rk4step        (  ipa) = 0.0
               csite%fmean_available_water(  ipa) = 0.0
               csite%fmean_can_theiv      (  ipa) = 0.0
               csite%fmean_can_theta      (  ipa) = 0.0
               csite%fmean_can_vpdef      (  ipa) = 0.0
               csite%fmean_can_temp       (  ipa) = 0.0
               csite%fmean_can_shv        (  ipa) = 0.0
               csite%fmean_can_co2        (  ipa) = 0.0
               csite%fmean_can_rhos       (  ipa) = 0.0
               csite%fmean_can_prss       (  ipa) = 0.0
               csite%fmean_gnd_temp       (  ipa) = 0.0
               csite%fmean_gnd_shv        (  ipa) = 0.0
               csite%fmean_can_ggnd       (  ipa) = 0.0
               csite%fmean_sfcw_depth     (  ipa) = 0.0
               csite%fmean_sfcw_energy    (  ipa) = 0.0
               csite%fmean_sfcw_mass      (  ipa) = 0.0
               csite%fmean_sfcw_temp      (  ipa) = 0.0
               csite%fmean_sfcw_fliq      (  ipa) = 0.0
               csite%fmean_soil_energy    (:,ipa) = 0.0
               csite%fmean_soil_mstpot    (:,ipa) = 0.0
               csite%fmean_soil_water     (:,ipa) = 0.0
               csite%fmean_soil_temp      (:,ipa) = 0.0
               csite%fmean_soil_fliq      (:,ipa) = 0.0
               csite%fmean_rshort_gnd     (  ipa) = 0.0
               csite%fmean_par_gnd        (  ipa) = 0.0
               csite%fmean_rlong_gnd      (  ipa) = 0.0
               csite%fmean_rlongup        (  ipa) = 0.0
               csite%fmean_parup          (  ipa) = 0.0
               csite%fmean_nirup          (  ipa) = 0.0
               csite%fmean_rshortup       (  ipa) = 0.0
               csite%fmean_rnet           (  ipa) = 0.0
               csite%fmean_albedo         (  ipa) = 0.0
               csite%fmean_albedo_par     (  ipa) = 0.0
               csite%fmean_albedo_nir     (  ipa) = 0.0
               csite%fmean_rlong_albedo   (  ipa) = 0.0
               csite%fmean_ustar          (  ipa) = 0.0
               csite%fmean_tstar          (  ipa) = 0.0
               csite%fmean_qstar          (  ipa) = 0.0
               csite%fmean_cstar          (  ipa) = 0.0
               csite%fmean_carbon_ac      (  ipa) = 0.0
               csite%fmean_carbon_st      (  ipa) = 0.0
               csite%fmean_vapor_gc       (  ipa) = 0.0
               csite%fmean_vapor_ac       (  ipa) = 0.0
               csite%fmean_smoist_gg      (:,ipa) = 0.0
               csite%fmean_throughfall    (  ipa) = 0.0
               csite%fmean_transloss      (:,ipa) = 0.0
               csite%fmean_runoff         (  ipa) = 0.0
               csite%fmean_drainage       (  ipa) = 0.0
               csite%fmean_sensible_gc    (  ipa) = 0.0
               csite%fmean_sensible_ac    (  ipa) = 0.0
               csite%fmean_sensible_gg    (:,ipa) = 0.0
               csite%fmean_qthroughfall   (  ipa) = 0.0
               csite%fmean_qrunoff        (  ipa) = 0.0
               csite%fmean_qdrainage      (  ipa) = 0.0
               !---------------------------------------------------------------------------!




               !----- Cohort level variables. ---------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts
                  cpatch%fmean_gpp               (ico) = 0.0
                  cpatch%fmean_npp               (ico) = 0.0
                  cpatch%fmean_leaf_resp         (ico) = 0.0
                  cpatch%fmean_root_resp         (ico) = 0.0
                  cpatch%fmean_growth_resp       (ico) = 0.0
                  cpatch%fmean_storage_resp      (ico) = 0.0
                  cpatch%fmean_vleaf_resp        (ico) = 0.0
                  cpatch%fmean_plresp            (ico) = 0.0
                  cpatch%fmean_leaf_energy       (ico) = 0.0
                  cpatch%fmean_leaf_water        (ico) = 0.0
                  cpatch%fmean_leaf_hcap         (ico) = 0.0
                  cpatch%fmean_leaf_vpdef        (ico) = 0.0
                  cpatch%fmean_leaf_temp         (ico) = 0.0
                  cpatch%fmean_leaf_fliq         (ico) = 0.0
                  cpatch%fmean_leaf_gsw          (ico) = 0.0
                  cpatch%fmean_leaf_gbw          (ico) = 0.0
                  cpatch%fmean_wood_energy       (ico) = 0.0
                  cpatch%fmean_wood_water        (ico) = 0.0
                  cpatch%fmean_wood_hcap         (ico) = 0.0
                  cpatch%fmean_wood_temp         (ico) = 0.0
                  cpatch%fmean_wood_fliq         (ico) = 0.0
                  cpatch%fmean_wood_gbw          (ico) = 0.0
                  cpatch%fmean_fs_open           (ico) = 0.0
                  cpatch%fmean_fsw               (ico) = 0.0
                  cpatch%fmean_fsn               (ico) = 0.0
                  cpatch%fmean_a_light           (ico) = 0.0
                  cpatch%fmean_a_rubp            (ico) = 0.0
                  cpatch%fmean_a_co2             (ico) = 0.0
                  cpatch%fmean_psi_open          (ico) = 0.0
                  cpatch%fmean_psi_closed        (ico) = 0.0
                  cpatch%fmean_water_supply      (ico) = 0.0
                  cpatch%fmean_light_level       (ico) = 0.0
                  cpatch%fmean_light_level_beam  (ico) = 0.0
                  cpatch%fmean_light_level_diff  (ico) = 0.0

                  cpatch%fmean_par_level_beam    (ico) = 0.0                  
                  cpatch%fmean_par_level_diffd   (ico) = 0.0                  
                  cpatch%fmean_par_level_diffu   (ico) = 0.0

                  cpatch%fmean_par_l             (ico) = 0.0
                  cpatch%fmean_par_l_beam        (ico) = 0.0
                  cpatch%fmean_par_l_diff        (ico) = 0.0
                  cpatch%fmean_rshort_l          (ico) = 0.0
                  cpatch%fmean_rlong_l           (ico) = 0.0
                  cpatch%fmean_sensible_lc       (ico) = 0.0
                  cpatch%fmean_vapor_lc          (ico) = 0.0
                  cpatch%fmean_transp            (ico) = 0.0
                  cpatch%fmean_intercepted_al    (ico) = 0.0
                  cpatch%fmean_wshed_lg          (ico) = 0.0
                  cpatch%fmean_rshort_w          (ico) = 0.0
                  cpatch%fmean_rlong_w           (ico) = 0.0
                  cpatch%fmean_rad_profile     (:,ico) = 0.0
                  cpatch%fmean_sensible_wc       (ico) = 0.0
                  cpatch%fmean_vapor_wc          (ico) = 0.0
                  cpatch%fmean_intercepted_aw    (ico) = 0.0
                  cpatch%fmean_wshed_wg          (ico) = 0.0
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!


      return
   end subroutine zero_ed_fmean_vars
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!










   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !                             |-------------------------------|                         !
   !                             |** DAILY AVERAGE SUBROUTINES **|                         !
   !                             |-------------------------------|                         !
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine integrates most of the daily averages.  This is called after the   !
   ! "fmean" variables are normalised, so we take advantage of these "fmean" variables.    !
   !                                                                                       !
   !    A few variables are _NOT_ integrated here:                                         !
   ! 1. Variables that should be integrated only during daylight hours                     !
   ! 2. The NPP breakdown variables are integrated in a separate routine                   !
   ! 3. Variables such as temperature and liquid fraction, which are found after the daily !
   !    means are normalised.                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine integrate_ed_dmean_vars(cgrid)
      use ed_state_vars        , only : edtype              & ! structure
                                      , polygontype         & ! structure
                                      , sitetype            & ! structure
                                      , patchtype           ! ! structure
      use ed_misc_coms         , only : frqsum              ! ! intent(in)
      use consts_coms          , only : day_sec             ! ! intent(in)
      implicit none

      !----- Argument ---------------------------------------------------------------------!
      type(edtype)      , target  :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype) , pointer :: cpoly
      type(sitetype)    , pointer :: csite
      type(patchtype)   , pointer :: cpatch
      integer                     :: ipy
      integer                     :: isi
      integer                     :: ipa
      integer                     :: ico
      !----- Locally saved variables. -----------------------------------------------------!
      logical                            , save       :: find_factors    = .true.
      real                               , save       :: frqsum_o_daysec = 1.e34
      !------------------------------------------------------------------------------------!


      !----- Compute the normalisation factors. This is done only once. -------------------!
      if (find_factors) then
         frqsum_o_daysec = frqsum / day_sec
         find_factors    = .false.
      end if
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Use the variables that have been already aggregated.                          !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         !---------------------------------------------------------------------------------!
         !    Integrate polygon-level variables.                                           !
         !---------------------------------------------------------------------------------!
         cgrid%dmean_gpp            (ipy) = cgrid%dmean_gpp            (ipy)               &
                                          + cgrid%fmean_gpp            (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_npp            (ipy) = cgrid%dmean_npp            (ipy)               &
                                          + cgrid%fmean_npp            (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_leaf_resp      (ipy) = cgrid%dmean_leaf_resp      (ipy)               &
                                          + cgrid%fmean_leaf_resp      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_root_resp      (ipy) = cgrid%dmean_root_resp      (ipy)               &
                                          + cgrid%fmean_root_resp      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_growth_resp    (ipy) = cgrid%dmean_growth_resp    (ipy)               &
                                          + cgrid%fmean_growth_resp    (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_storage_resp   (ipy) = cgrid%dmean_storage_resp   (ipy)               &
                                          + cgrid%fmean_storage_resp   (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_vleaf_resp     (ipy) = cgrid%dmean_vleaf_resp     (ipy)               &
                                          + cgrid%fmean_vleaf_resp     (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_plresp         (ipy) = cgrid%dmean_plresp         (ipy)               &
                                          + cgrid%fmean_plresp         (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_a_light        (ipy) = cgrid%dmean_a_light        (ipy)               &
                                          + cgrid%fmean_a_light        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_a_rubp         (ipy) = cgrid%dmean_a_rubp         (ipy)               &
                                          + cgrid%fmean_a_rubp         (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_a_co2          (ipy) = cgrid%dmean_a_co2          (ipy)               &
                                          + cgrid%fmean_a_co2          (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_leaf_energy    (ipy) = cgrid%dmean_leaf_energy    (ipy)               &
                                          + cgrid%fmean_leaf_energy    (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_leaf_water     (ipy) = cgrid%dmean_leaf_water     (ipy)               &
                                          + cgrid%fmean_leaf_water     (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_leaf_hcap      (ipy) = cgrid%dmean_leaf_hcap      (ipy)               &
                                          + cgrid%fmean_leaf_hcap      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_leaf_vpdef     (ipy) = cgrid%dmean_leaf_vpdef     (ipy)               &
                                          + cgrid%fmean_leaf_vpdef     (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_leaf_gsw       (ipy) = cgrid%dmean_leaf_gsw       (ipy)               &
                                          + cgrid%fmean_leaf_gsw       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_leaf_gbw       (ipy) = cgrid%dmean_leaf_gbw       (ipy)               &
                                          + cgrid%fmean_leaf_gbw       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_wood_energy    (ipy) = cgrid%dmean_wood_energy    (ipy)               &
                                          + cgrid%fmean_wood_energy    (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_wood_water     (ipy) = cgrid%dmean_wood_water     (ipy)               &
                                          + cgrid%fmean_wood_water     (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_wood_hcap      (ipy) = cgrid%dmean_wood_hcap      (ipy)               &
                                          + cgrid%fmean_wood_hcap      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_wood_gbw       (ipy) = cgrid%dmean_wood_gbw       (ipy)               &
                                          + cgrid%fmean_wood_gbw       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_psi_open       (ipy) = cgrid%dmean_psi_open       (ipy)               &
                                          + cgrid%fmean_psi_open       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_psi_closed     (ipy) = cgrid%dmean_psi_closed     (ipy)               &
                                          + cgrid%fmean_psi_closed     (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_water_supply   (ipy) = cgrid%dmean_water_supply   (ipy)               &
                                          + cgrid%fmean_water_supply   (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_par_l          (ipy) = cgrid%dmean_par_l          (ipy)               &
                                          + cgrid%fmean_par_l          (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_par_l_beam     (ipy) = cgrid%dmean_par_l_beam     (ipy)               &
                                          + cgrid%fmean_par_l_beam     (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_par_l_diff     (ipy) = cgrid%dmean_par_l_diff     (ipy)               &
                                          + cgrid%fmean_par_l_diff     (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_rshort_l       (ipy) = cgrid%dmean_rshort_l       (ipy)               &
                                          + cgrid%fmean_rshort_l       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_rlong_l        (ipy) = cgrid%dmean_rlong_l        (ipy)               &
                                          + cgrid%fmean_rlong_l        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_sensible_lc    (ipy) = cgrid%dmean_sensible_lc    (ipy)               &
                                          + cgrid%fmean_sensible_lc    (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_vapor_lc       (ipy) = cgrid%dmean_vapor_lc       (ipy)               &
                                          + cgrid%fmean_vapor_lc       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_transp         (ipy) = cgrid%dmean_transp         (ipy)               &
                                          + cgrid%fmean_transp         (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_intercepted_al (ipy) = cgrid%dmean_intercepted_al (ipy)               &
                                          + cgrid%fmean_intercepted_al (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_wshed_lg       (ipy) = cgrid%dmean_wshed_lg       (ipy)               &
                                          + cgrid%fmean_wshed_lg       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_rshort_w       (ipy) = cgrid%dmean_rshort_w       (ipy)               &
                                          + cgrid%fmean_rshort_w       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_rlong_w        (ipy) = cgrid%dmean_rlong_w        (ipy)               &
                                          + cgrid%fmean_rlong_w        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_sensible_wc    (ipy) = cgrid%dmean_sensible_wc    (ipy)               &
                                          + cgrid%fmean_sensible_wc    (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_vapor_wc       (ipy) = cgrid%dmean_vapor_wc       (ipy)               &
                                          + cgrid%fmean_vapor_wc       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_intercepted_aw (ipy) = cgrid%dmean_intercepted_aw (ipy)               &
                                          + cgrid%fmean_intercepted_aw (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_wshed_wg       (ipy) = cgrid%dmean_wshed_wg       (ipy)               &
                                          + cgrid%fmean_wshed_wg       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_rh             (ipy) = cgrid%dmean_rh             (ipy)               &
                                          + cgrid%fmean_rh             (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_cwd_rh         (ipy) = cgrid%dmean_cwd_rh         (ipy)               &
                                          + cgrid%fmean_cwd_rh         (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_nep            (ipy) = cgrid%dmean_nep            (ipy)               &
                                          + cgrid%fmean_nep            (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_rk4step        (ipy) = cgrid%dmean_rk4step        (ipy)               &
                                          + cgrid%fmean_rk4step        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_available_water(ipy) = cgrid%dmean_available_water(ipy)               &
                                          + cgrid%fmean_available_water(ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_can_theiv      (ipy) = cgrid%dmean_can_theiv      (ipy)               &
                                          + cgrid%fmean_can_theiv      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_can_theta      (ipy) = cgrid%dmean_can_theta      (ipy)               &
                                          + cgrid%fmean_can_theta      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_can_vpdef      (ipy) = cgrid%dmean_can_vpdef      (ipy)               &
                                          + cgrid%fmean_can_vpdef      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_can_temp       (ipy) = cgrid%dmean_can_temp       (ipy)               &
                                          + cgrid%fmean_can_temp       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_can_shv        (ipy) = cgrid%dmean_can_shv        (ipy)               &
                                          + cgrid%fmean_can_shv        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_can_co2        (ipy) = cgrid%dmean_can_co2        (ipy)               &
                                          + cgrid%fmean_can_co2        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_can_rhos       (ipy) = cgrid%dmean_can_rhos       (ipy)               &
                                          + cgrid%fmean_can_rhos       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_can_prss       (ipy) = cgrid%dmean_can_prss       (ipy)               &
                                          + cgrid%fmean_can_prss       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_gnd_temp       (ipy) = cgrid%dmean_gnd_temp       (ipy)               &
                                          + cgrid%fmean_gnd_temp       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_gnd_shv        (ipy) = cgrid%dmean_gnd_shv        (ipy)               &
                                          + cgrid%fmean_gnd_shv        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_can_ggnd       (ipy) = cgrid%dmean_can_ggnd       (ipy)               &
                                          + cgrid%fmean_can_ggnd       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_sfcw_depth     (ipy) = cgrid%dmean_sfcw_depth     (ipy)               &
                                          + cgrid%fmean_sfcw_depth     (ipy)               &
                                          * frqsum_o_daysec
         !------ Integrate the extensive version of temporary surface water energy. -------!
         cgrid%dmean_sfcw_energy    (ipy) = cgrid%dmean_sfcw_energy    (ipy)               &
                                          + cgrid%fmean_sfcw_energy    (ipy)               &
                                          * cgrid%fmean_sfcw_mass      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_sfcw_mass      (ipy) = cgrid%dmean_sfcw_mass      (ipy)               &
                                          + cgrid%fmean_sfcw_mass      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_rshort_gnd     (ipy) = cgrid%dmean_rshort_gnd     (ipy)               &
                                          + cgrid%fmean_rshort_gnd     (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_par_gnd        (ipy) = cgrid%dmean_par_gnd        (ipy)               &
                                          + cgrid%fmean_par_gnd        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_rlong_gnd      (ipy) = cgrid%dmean_rlong_gnd      (ipy)               &
                                          + cgrid%fmean_rlong_gnd      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_rlongup        (ipy) = cgrid%dmean_rlongup        (ipy)               &
                                          + cgrid%fmean_rlongup        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_parup          (ipy) = cgrid%dmean_parup          (ipy)               &
                                          + cgrid%fmean_parup          (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_nirup          (ipy) = cgrid%dmean_nirup          (ipy)               &
                                          + cgrid%fmean_nirup          (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_rshortup       (ipy) = cgrid%dmean_rshortup       (ipy)               &
                                          + cgrid%fmean_rshortup       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_rnet           (ipy) = cgrid%dmean_rnet           (ipy)               &
                                          + cgrid%fmean_rnet           (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_rlong_albedo   (ipy) = cgrid%dmean_rlong_albedo   (ipy)               &
                                          + cgrid%fmean_rlong_albedo   (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_ustar          (ipy) = cgrid%dmean_ustar          (ipy)               &
                                          + cgrid%fmean_ustar          (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_tstar          (ipy) = cgrid%dmean_tstar          (ipy)               &
                                          + cgrid%fmean_tstar          (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_qstar          (ipy) = cgrid%dmean_qstar          (ipy)               &
                                          + cgrid%fmean_qstar          (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_cstar          (ipy) = cgrid%dmean_cstar          (ipy)               &
                                          + cgrid%fmean_cstar          (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_carbon_ac      (ipy) = cgrid%dmean_carbon_ac      (ipy)               &
                                          + cgrid%fmean_carbon_ac      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_carbon_st      (ipy) = cgrid%dmean_carbon_st      (ipy)               &
                                          + cgrid%fmean_carbon_st      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_vapor_gc       (ipy) = cgrid%dmean_vapor_gc       (ipy)               &
                                          + cgrid%fmean_vapor_gc       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_vapor_ac       (ipy) = cgrid%dmean_vapor_ac       (ipy)               &
                                          + cgrid%fmean_vapor_ac       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_throughfall    (ipy) = cgrid%dmean_throughfall    (ipy)               &
                                          + cgrid%fmean_throughfall    (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_runoff         (ipy) = cgrid%dmean_runoff         (ipy)               &
                                          + cgrid%fmean_runoff         (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_drainage       (ipy) = cgrid%dmean_drainage       (ipy)               &
                                          + cgrid%fmean_drainage       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_sensible_gc    (ipy) = cgrid%dmean_sensible_gc    (ipy)               &
                                          + cgrid%fmean_sensible_gc    (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_sensible_ac    (ipy) = cgrid%dmean_sensible_ac    (ipy)               &
                                          + cgrid%fmean_sensible_ac    (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_qthroughfall   (ipy) = cgrid%dmean_qthroughfall   (ipy)               &
                                          + cgrid%fmean_qthroughfall   (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_qrunoff        (ipy) = cgrid%dmean_qrunoff        (ipy)               &
                                          + cgrid%fmean_qrunoff        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_qdrainage      (ipy) = cgrid%dmean_qdrainage      (ipy)               &
                                          + cgrid%fmean_qdrainage      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_atm_theiv      (ipy) = cgrid%dmean_atm_theiv      (ipy)               &
                                          + cgrid%fmean_atm_theiv      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_atm_theta      (ipy) = cgrid%dmean_atm_theta      (ipy)               &
                                          + cgrid%fmean_atm_theta      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_atm_vpdef      (ipy) = cgrid%dmean_atm_vpdef      (ipy)               &
                                          + cgrid%fmean_atm_vpdef      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_atm_shv        (ipy) = cgrid%dmean_atm_shv        (ipy)               &
                                          + cgrid%fmean_atm_shv        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_atm_rshort     (ipy) = cgrid%dmean_atm_rshort     (ipy)               &
                                          + cgrid%fmean_atm_rshort     (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_atm_rshort_diff(ipy) = cgrid%dmean_atm_rshort_diff(ipy)               &
                                          + cgrid%fmean_atm_rshort_diff(ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_atm_par        (ipy) = cgrid%dmean_atm_par        (ipy)               &
                                          + cgrid%fmean_atm_par        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_atm_par_diff   (ipy) = cgrid%dmean_atm_par_diff   (ipy)               &
                                          + cgrid%fmean_atm_par_diff   (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_atm_rlong      (ipy) = cgrid%dmean_atm_rlong      (ipy)               &
                                          + cgrid%fmean_atm_rlong      (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_atm_vels       (ipy) = cgrid%dmean_atm_vels       (ipy)               &
                                          + cgrid%fmean_atm_vels       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_atm_prss       (ipy) = cgrid%dmean_atm_prss       (ipy)               &
                                          + cgrid%fmean_atm_prss       (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_atm_co2        (ipy) = cgrid%dmean_atm_co2        (ipy)               &
                                          + cgrid%fmean_atm_co2        (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_pcpg           (ipy) = cgrid%dmean_pcpg           (ipy)               &
                                          + cgrid%fmean_pcpg           (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_qpcpg          (ipy) = cgrid%dmean_qpcpg          (ipy)               &
                                          + cgrid%fmean_qpcpg          (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_dpcpg          (ipy) = cgrid%dmean_dpcpg          (ipy)               &
                                          + cgrid%fmean_dpcpg          (ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_soil_energy  (:,ipy) = cgrid%dmean_soil_energy  (:,ipy)               &
                                          + cgrid%fmean_soil_energy  (:,ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_soil_mstpot  (:,ipy) = cgrid%dmean_soil_mstpot  (:,ipy)               &
                                          + cgrid%fmean_soil_mstpot  (:,ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_soil_water   (:,ipy) = cgrid%dmean_soil_water   (:,ipy)               &
                                          + cgrid%fmean_soil_water   (:,ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_smoist_gg    (:,ipy) = cgrid%dmean_smoist_gg    (:,ipy)               &
                                          + cgrid%fmean_smoist_gg    (:,ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_transloss    (:,ipy) = cgrid%dmean_transloss    (:,ipy)               &
                                          + cgrid%fmean_transloss    (:,ipy)               &
                                          * frqsum_o_daysec
         cgrid%dmean_sensible_gg  (:,ipy) = cgrid%dmean_sensible_gg  (:,ipy)               &
                                          + cgrid%fmean_sensible_gg  (:,ipy)               &
                                          * frqsum_o_daysec
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !      Site loop.                                                                 !
         !---------------------------------------------------------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)


            !------------------------------------------------------------------------------!
            !    Integrate site-level variables.                                           !
            !------------------------------------------------------------------------------!
            cpoly%dmean_atm_theiv      (isi) = cpoly%dmean_atm_theiv      (isi)            &
                                             + cpoly%fmean_atm_theiv      (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_theta      (isi) = cpoly%dmean_atm_theta      (isi)            &
                                             + cpoly%fmean_atm_theta      (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_temp       (isi) = cpoly%dmean_atm_temp       (isi)            &
                                             + cpoly%fmean_atm_temp       (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_vpdef      (isi) = cpoly%dmean_atm_vpdef      (isi)            &
                                             + cpoly%fmean_atm_vpdef      (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_shv        (isi) = cpoly%dmean_atm_shv        (isi)            &
                                             + cpoly%fmean_atm_shv        (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_rshort     (isi) = cpoly%dmean_atm_rshort     (isi)            &
                                             + cpoly%fmean_atm_rshort     (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_rshort_diff(isi) = cpoly%dmean_atm_rshort_diff(isi)            &
                                             + cpoly%fmean_atm_rshort_diff(isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_par        (isi) = cpoly%dmean_atm_par        (isi)            &
                                             + cpoly%fmean_atm_par        (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_par_diff   (isi) = cpoly%dmean_atm_par_diff   (isi)            &
                                             + cpoly%fmean_atm_par_diff   (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_rlong      (isi) = cpoly%dmean_atm_rlong      (isi)            &
                                             + cpoly%fmean_atm_rlong      (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_vels       (isi) = cpoly%dmean_atm_vels       (isi)            &
                                             + cpoly%fmean_atm_vels       (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_rhos       (isi) = cpoly%dmean_atm_rhos       (isi)            &
                                             + cpoly%fmean_atm_rhos       (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_prss       (isi) = cpoly%dmean_atm_prss       (isi)            &
                                             + cpoly%fmean_atm_prss       (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_atm_co2        (isi) = cpoly%dmean_atm_co2        (isi)            &
                                             + cpoly%fmean_atm_co2        (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_pcpg           (isi) = cpoly%dmean_pcpg           (isi)            &
                                             + cpoly%fmean_pcpg           (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_qpcpg          (isi) = cpoly%dmean_qpcpg          (isi)            &
                                             + cpoly%fmean_qpcpg          (isi)            &
                                             * frqsum_o_daysec
            cpoly%dmean_dpcpg          (isi) = cpoly%dmean_dpcpg          (isi)            &
                                             + cpoly%fmean_dpcpg          (isi)            &
                                             * frqsum_o_daysec
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Patch loop.                                                             !
            !------------------------------------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa) 


               !---------------------------------------------------------------------------!
               !      Integrate patch-level variables.                                     !
               !---------------------------------------------------------------------------!
               csite%dmean_co2_residual      (ipa) = csite%dmean_co2_residual      (ipa)   &
                                                   + csite%co2budget_residual      (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_energy_residual   (ipa) = csite%dmean_energy_residual   (ipa)   &
                                                   + csite%ebudget_residual        (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_water_residual    (ipa) = csite%dmean_water_residual    (ipa)   &
                                                   + csite%wbudget_residual        (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_rh                (ipa) = csite%dmean_rh                (ipa)   &
                                                   + csite%fmean_rh                (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_cwd_rh            (ipa) = csite%dmean_cwd_rh            (ipa)   &
                                                   + csite%fmean_cwd_rh            (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_nep               (ipa) = csite%dmean_nep               (ipa)   &
                                                   + csite%fmean_nep               (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_rk4step           (ipa) = csite%dmean_rk4step           (ipa)   &
                                                   + csite%fmean_rk4step           (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_available_water   (ipa) = csite%dmean_available_water   (ipa)   &
                                                   + csite%fmean_available_water   (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_can_theiv         (ipa) = csite%dmean_can_theiv         (ipa)   &
                                                   + csite%fmean_can_theiv         (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_can_theta         (ipa) = csite%dmean_can_theta         (ipa)   &
                                                   + csite%fmean_can_theta         (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_can_vpdef         (ipa) = csite%dmean_can_vpdef         (ipa)   &
                                                   + csite%fmean_can_vpdef         (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_can_shv           (ipa) = csite%dmean_can_shv           (ipa)   &
                                                   + csite%fmean_can_shv           (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_can_co2           (ipa) = csite%dmean_can_co2           (ipa)   &
                                                   + csite%fmean_can_co2           (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_can_prss          (ipa) = csite%dmean_can_prss          (ipa)   &
                                                   + csite%fmean_can_prss          (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_gnd_temp          (ipa) = csite%dmean_gnd_temp          (ipa)   &
                                                   + csite%fmean_gnd_temp          (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_gnd_shv           (ipa) = csite%dmean_gnd_shv           (ipa)   &
                                                   + csite%fmean_gnd_shv           (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_can_ggnd          (ipa) = csite%dmean_can_ggnd          (ipa)   &
                                                   + csite%fmean_can_ggnd          (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_sfcw_depth        (ipa) = csite%dmean_sfcw_depth        (ipa)   &
                                                   + csite%fmean_sfcw_depth        (ipa)   &
                                                   * frqsum_o_daysec
               !----- Temporarily make the pounding internal energy extensive [J/m2]. -----!
               csite%dmean_sfcw_energy       (ipa) = csite%dmean_sfcw_energy       (ipa)   &
                                                   + csite%fmean_sfcw_energy       (ipa)   &
                                                   * csite%fmean_sfcw_mass         (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_sfcw_mass         (ipa) = csite%dmean_sfcw_mass         (ipa)   &
                                                   + csite%fmean_sfcw_mass         (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_rshort_gnd        (ipa) = csite%dmean_rshort_gnd        (ipa)   &
                                                   + csite%fmean_rshort_gnd        (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_par_gnd           (ipa) = csite%dmean_par_gnd           (ipa)   &
                                                   + csite%fmean_par_gnd           (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_rlong_gnd         (ipa) = csite%dmean_rlong_gnd         (ipa)   &
                                                   + csite%fmean_rlong_gnd         (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_rlongup           (ipa) = csite%dmean_rlongup           (ipa)   &
                                                   + csite%fmean_rlongup           (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_parup             (ipa) = csite%dmean_parup             (ipa)   &
                                                   + csite%fmean_parup             (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_nirup             (ipa) = csite%dmean_nirup             (ipa)   &
                                                   + csite%fmean_nirup             (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_rshortup          (ipa) = csite%dmean_rshortup          (ipa)   &
                                                   + csite%fmean_rshortup          (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_rnet              (ipa) = csite%dmean_rnet              (ipa)   &
                                                   + csite%fmean_rnet              (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_rlong_albedo      (ipa) = csite%dmean_rlong_albedo      (ipa)   &
                                                   + csite%fmean_rlong_albedo      (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_ustar             (ipa) = csite%dmean_ustar             (ipa)   &
                                                   + csite%fmean_ustar             (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_tstar             (ipa) = csite%dmean_tstar             (ipa)   &
                                                   + csite%fmean_tstar             (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_qstar             (ipa) = csite%dmean_qstar             (ipa)   &
                                                   + csite%fmean_qstar             (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_cstar             (ipa) = csite%dmean_cstar             (ipa)   &
                                                   + csite%fmean_cstar             (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_carbon_ac         (ipa) = csite%dmean_carbon_ac         (ipa)   &
                                                   + csite%fmean_carbon_ac         (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_carbon_st         (ipa) = csite%dmean_carbon_st         (ipa)   &
                                                   + csite%fmean_carbon_st         (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_vapor_gc          (ipa) = csite%dmean_vapor_gc          (ipa)   &
                                                   + csite%fmean_vapor_gc          (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_vapor_ac          (ipa) = csite%dmean_vapor_ac          (ipa)   &
                                                   + csite%fmean_vapor_ac          (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_throughfall       (ipa) = csite%dmean_throughfall       (ipa)   &
                                                   + csite%fmean_throughfall       (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_runoff            (ipa) = csite%dmean_runoff            (ipa)   &
                                                   + csite%fmean_runoff            (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_drainage          (ipa) = csite%dmean_drainage          (ipa)   &
                                                   + csite%fmean_drainage          (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_sensible_gc       (ipa) = csite%dmean_sensible_gc       (ipa)   &
                                                   + csite%fmean_sensible_gc       (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_sensible_ac       (ipa) = csite%dmean_sensible_ac       (ipa)   &
                                                   + csite%fmean_sensible_ac       (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_sensible_gg     (:,ipa) = csite%dmean_sensible_gg     (:,ipa)   &
                                                   + csite%fmean_sensible_gg     (:,ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_qthroughfall      (ipa) = csite%dmean_qthroughfall      (ipa)   &
                                                   + csite%fmean_qthroughfall      (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_qrunoff           (ipa) = csite%dmean_qrunoff           (ipa)   &
                                                   + csite%fmean_qrunoff           (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_qdrainage         (ipa) = csite%dmean_qdrainage         (ipa)   &
                                                   + csite%fmean_qdrainage         (ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_soil_energy     (:,ipa) = csite%dmean_soil_energy     (:,ipa)   &
                                                   + csite%fmean_soil_energy     (:,ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_soil_mstpot     (:,ipa) = csite%dmean_soil_mstpot     (:,ipa)   &
                                                   + csite%fmean_soil_mstpot     (:,ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_soil_water      (:,ipa) = csite%dmean_soil_water      (:,ipa)   &
                                                   + csite%fmean_soil_water      (:,ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_smoist_gg       (:,ipa) = csite%dmean_smoist_gg       (:,ipa)   &
                                                   + csite%fmean_smoist_gg       (:,ipa)   &
                                                   * frqsum_o_daysec
               csite%dmean_transloss       (:,ipa) = csite%dmean_transloss       (:,ipa)   &
                                                   + csite%fmean_transloss       (:,ipa)   &
                                                   * frqsum_o_daysec
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Patch loop.                                                          !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts
                  !------------------------------------------------------------------------!
                  !      Integrate cohort-level variables.                                 !
                  !------------------------------------------------------------------------!
                  cpatch%dmean_gpp           (ico) = cpatch%dmean_gpp           (ico)      &
                                                   + cpatch%fmean_gpp           (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_npp           (ico) = cpatch%dmean_npp           (ico)      &
                                                   + cpatch%fmean_npp           (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_leaf_resp     (ico) = cpatch%dmean_leaf_resp     (ico)      &
                                                   + cpatch%fmean_leaf_resp     (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_root_resp     (ico) = cpatch%dmean_root_resp     (ico)      &
                                                   + cpatch%fmean_root_resp     (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_growth_resp   (ico) = cpatch%dmean_growth_resp   (ico)      &
                                                   + cpatch%fmean_growth_resp   (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_storage_resp  (ico) = cpatch%dmean_storage_resp  (ico)      &
                                                   + cpatch%fmean_storage_resp  (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_vleaf_resp    (ico) = cpatch%dmean_vleaf_resp    (ico)      &
                                                   + cpatch%fmean_vleaf_resp    (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_plresp        (ico) = cpatch%dmean_plresp        (ico)      &
                                                   + cpatch%fmean_plresp        (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_a_light       (ico) = cpatch%dmean_a_light       (ico)      &
                                                   + cpatch%fmean_a_light       (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_a_rubp        (ico) = cpatch%dmean_a_rubp        (ico)      &
                                                   + cpatch%fmean_a_rubp        (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_a_co2         (ico) = cpatch%dmean_a_co2         (ico)      &
                                                   + cpatch%fmean_a_co2         (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_leaf_energy   (ico) = cpatch%dmean_leaf_energy   (ico)      &
                                                   + cpatch%fmean_leaf_energy   (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_leaf_water    (ico) = cpatch%dmean_leaf_water    (ico)      &
                                                   + cpatch%fmean_leaf_water    (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_leaf_hcap     (ico) = cpatch%dmean_leaf_hcap     (ico)      &
                                                   + cpatch%fmean_leaf_hcap     (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_leaf_vpdef    (ico) = cpatch%dmean_leaf_vpdef    (ico)      &
                                                   + cpatch%fmean_leaf_vpdef    (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_leaf_gsw      (ico) = cpatch%dmean_leaf_gsw      (ico)      &
                                                   + cpatch%fmean_leaf_gsw      (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_leaf_gbw      (ico) = cpatch%dmean_leaf_gbw      (ico)      &
                                                   + cpatch%fmean_leaf_gbw      (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_wood_energy   (ico) = cpatch%dmean_wood_energy   (ico)      &
                                                   + cpatch%fmean_wood_energy   (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_wood_water    (ico) = cpatch%dmean_wood_water    (ico)      &
                                                   + cpatch%fmean_wood_water    (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_wood_hcap     (ico) = cpatch%dmean_wood_hcap     (ico)      &
                                                   + cpatch%fmean_wood_hcap     (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_wood_gbw      (ico) = cpatch%dmean_wood_gbw      (ico)      &
                                                   + cpatch%fmean_wood_gbw      (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_psi_open      (ico) = cpatch%dmean_psi_open      (ico)      &
                                                   + cpatch%fmean_psi_open      (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_psi_closed    (ico) = cpatch%dmean_psi_closed    (ico)      &
                                                   + cpatch%fmean_psi_closed    (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_water_supply  (ico) = cpatch%dmean_water_supply  (ico)      &
                                                   + cpatch%fmean_water_supply  (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_par_l         (ico) = cpatch%dmean_par_l         (ico)      &
                                                   + cpatch%fmean_par_l         (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_par_l_beam    (ico) = cpatch%dmean_par_l_beam    (ico)      &
                                                   + cpatch%fmean_par_l_beam    (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_par_l_diff    (ico) = cpatch%dmean_par_l_diff    (ico)      &
                                                   + cpatch%fmean_par_l_diff    (ico)      &
                                                   * frqsum_o_daysec

                  cpatch%dmean_par_level_beam(ico) = cpatch%dmean_par_level_beam(ico)      &
                                                   + cpatch%fmean_par_level_beam(ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_par_level_diffd(ico)= cpatch%dmean_par_level_diffd(ico)     &
                                                   + cpatch%fmean_par_level_diffd(ico)     &
                                                   * frqsum_o_daysec
                  cpatch%dmean_par_level_diffu(ico) = cpatch%dmean_par_level_diffu(ico)    &
                                                   + cpatch%fmean_par_level_diffu(ico)     &
                                                   * frqsum_o_daysec

                  cpatch%dmean_rshort_l      (ico) = cpatch%dmean_rshort_l      (ico)      &
                                                   + cpatch%fmean_rshort_l      (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_rlong_l       (ico) = cpatch%dmean_rlong_l       (ico)      &
                                                   + cpatch%fmean_rlong_l       (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_sensible_lc   (ico) = cpatch%dmean_sensible_lc   (ico)      &
                                                   + cpatch%fmean_sensible_lc   (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_vapor_lc      (ico) = cpatch%dmean_vapor_lc      (ico)      &
                                                   + cpatch%fmean_vapor_lc      (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_transp        (ico) = cpatch%dmean_transp        (ico)      &
                                                   + cpatch%fmean_transp        (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_intercepted_al(ico) = cpatch%dmean_intercepted_al(ico)      &
                                                   + cpatch%fmean_intercepted_al(ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_wshed_lg      (ico) = cpatch%dmean_wshed_lg      (ico)      &
                                                   + cpatch%fmean_wshed_lg      (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_rshort_w      (ico) = cpatch%dmean_rshort_w      (ico)      &
                                                   + cpatch%fmean_rshort_w      (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_rlong_w       (ico) = cpatch%dmean_rlong_w       (ico)      &
                                                   + cpatch%fmean_rlong_w       (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_rad_profile (:,ico) = cpatch%dmean_rad_profile (:,ico)      &
                                                   + cpatch%fmean_rad_profile (:,ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_sensible_wc   (ico) = cpatch%dmean_sensible_wc   (ico)      &
                                                   + cpatch%fmean_sensible_wc   (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_vapor_wc      (ico) = cpatch%dmean_vapor_wc      (ico)      &
                                                   + cpatch%fmean_vapor_wc      (ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_intercepted_aw(ico) = cpatch%dmean_intercepted_aw(ico)      &
                                                   + cpatch%fmean_intercepted_aw(ico)      &
                                                   * frqsum_o_daysec
                  cpatch%dmean_wshed_wg      (ico) = cpatch%dmean_wshed_wg      (ico)      &
                                                   + cpatch%fmean_wshed_wg      (ico)      &
                                                   * frqsum_o_daysec
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!
      return
   end subroutine integrate_ed_dmean_vars
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will scale the daily averages of GPP and some respiration vari-   !
   ! ables to normal units.  These variables are not for output, so they are done          !
   ! separatedly.  There are also some output variables here, because these depend on the  !
   ! average of the gpp, and leaf and root respiration and would need to be calculated     !
   ! again otherwise.                                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine normalize_ed_today_vars(cgrid)
      use ed_state_vars , only : edtype        & ! structure
                               , polygontype   & ! structure
                               , sitetype      & ! structure
                               , patchtype     ! ! structure
      use ed_max_dims   , only : n_pft         & ! intent(in)
                               , n_age         & ! intent(in)
                               , n_dbh         ! ! intent(in)
      use ed_misc_coms  , only : writing_long  & ! intent(in)
                               , writing_eorq  & ! intent(in)
                               , writing_dcyc  & ! intent(in)
                               , dtlsm         ! ! intent(in)
      use consts_coms   , only : umols_2_kgCyr & ! intent(in)
                               , day_sec       ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target     :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer    :: cpoly
      type(sitetype)   , pointer    :: csite
      type(patchtype)  , pointer    :: cpatch
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      real                          :: dtlsm_o_daysec
      !------------------------------------------------------------------------------------!


      !----- Find the time scale factor. --------------------------------------------------!
      dtlsm_o_daysec = dtlsm / day_sec
      !------------------------------------------------------------------------------------!


      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         !----- This part is done only if arrays are sought. ------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            patchloop: do ipa=1,csite%npatches

               csite%today_A_decomp (ipa) = csite%today_A_decomp(ipa)  * dtlsm_o_daysec
               csite%today_Af_decomp(ipa) = csite%today_Af_decomp(ipa) * dtlsm_o_daysec

               !----- Copy the decomposition terms to the daily mean if they are sought. --!
               if (writing_long) then
                  csite%dmean_A_decomp(ipa)  = csite%today_A_decomp(ipa)
                  csite%dmean_Af_decomp(ipa) = csite%today_Af_decomp(ipa)
               end if

               cpatch => csite%patch(ipa)
               
               !----- Included a loop so it won't crash with empty cohorts... -------------!
               cohortloop: do ico=1,cpatch%ncohorts
                  !------------------------------------------------------------------------!
                  !     Normalise the variables used to compute carbon balance.            !
                  !------------------------------------------------------------------------!
                  cpatch%today_gpp          (ico) = cpatch%today_gpp          (ico)        &
                                                  * dtlsm_o_daysec
                  cpatch%today_gpp_pot      (ico) = cpatch%today_gpp_pot      (ico)        &
                                                  * dtlsm_o_daysec
                  cpatch%today_gpp_lightmax (ico) = cpatch%today_gpp_lightmax (ico)        &
                                                  * dtlsm_o_daysec
                  cpatch%today_gpp_moistmax (ico) = cpatch%today_gpp_moistmax (ico)        &
                                                  * dtlsm_o_daysec
                  cpatch%today_leaf_resp    (ico) = cpatch%today_leaf_resp    (ico)        &
                                                  * dtlsm_o_daysec
                  cpatch%today_root_resp    (ico) = cpatch%today_root_resp    (ico)        &
                                                  * dtlsm_o_daysec
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine normalize_ed_today_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will scale the daily NPP allocation terms                         !
   !---------------------------------------------------------------------------------------!
   subroutine normalize_ed_todayNPP_vars(cgrid)
      use ed_state_vars , only : edtype        & ! structure
                               , polygontype   & ! structure
                               , sitetype      & ! structure
                               , patchtype     ! ! structure
      use ed_misc_coms  , only : writing_long  & ! intent(in)
                               , writing_eorq  & ! intent(in)
                               , writing_dcyc  ! ! intent(in)
      use consts_coms   , only : yr_day ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)                                  , target     :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype)                             , pointer    :: cpoly
      type(sitetype)                                , pointer    :: csite
      type(patchtype)                               , pointer    :: cpatch
      integer                                                    :: ipy
      integer                                                    :: isi
      integer                                                    :: ipa
      integer                                                    :: ico

      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            patchloop: do ipa=1,csite%npatches

               cpatch => csite%patch(ipa)
               
               !----- Included a loop so it won't crash with empty cohorts... -------------!
               cohortloop: do ico=1,cpatch%ncohorts

                  !------------------------------------------------------------------------!
                  !    We now update the daily means of NPP allocation terms               !
                  ! and we convert them to kgC/plant/yr                                    !
                  !------------------------------------------------------------------------!
                  if (writing_long) then
                     cpatch%dmean_nppleaf   (ico) = cpatch%today_nppleaf   (ico)           &
                                                  * yr_day / cpatch%nplant (ico)
                     cpatch%dmean_nppfroot  (ico) = cpatch%today_nppfroot  (ico)           &
                                                  * yr_day / cpatch%nplant (ico)
                     cpatch%dmean_nppsapwood(ico) = cpatch%today_nppsapwood(ico)           &
                                                  * yr_day / cpatch%nplant (ico)
                     cpatch%dmean_nppcroot  (ico) = cpatch%today_nppcroot  (ico)           &
                                                  * yr_day / cpatch%nplant (ico)
                     cpatch%dmean_nppseeds  (ico) = cpatch%today_nppseeds  (ico)           &
                                                  * yr_day / cpatch%nplant (ico)
                     cpatch%dmean_nppwood   (ico) = cpatch%today_nppwood   (ico)           &
                                                  * yr_day / cpatch%nplant (ico)
                     cpatch%dmean_nppdaily  (ico) = cpatch%today_nppdaily  (ico)           &
                                                  * yr_day / cpatch%nplant (ico)
                  end if
               end do cohortloop
            end do patchloop
         end do siteloop
      end do polyloop

      return
   end subroutine normalize_ed_todayNPP_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine normalises the daily mean variables of those variables that could  !
   ! not be integrated directly.  This includes temperatures, polygon-level budget vari-   !
   ! ables, and variables that are defined during daytime only.                            !
   !---------------------------------------------------------------------------------------!
   subroutine normalize_ed_dmean_vars(cgrid)
      use ed_state_vars        , only : edtype             & ! structure
                                      , polygontype        & ! structure
                                      , sitetype           & ! structure
                                      , patchtype          ! ! structure
      use grid_coms            , only : nzg                ! ! intent(in)
      use ed_misc_coms         , only : dtlsm              ! ! intent(in)
      use therm_lib            , only : press2exner        & ! function
                                      , extheta2temp       & ! function
                                      , uextcm2tl          & ! subroutine
                                      , uint2tl            & ! subroutine
                                      , idealdenssh        ! ! function
      use soil_coms            , only : tiny_sfcwater_mass & ! intent(in)
                                      , soil               ! ! intent(in)
      use consts_coms          , only : t00                & ! intent(in)
                                      , wdns               ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)                       , target     :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype)                  , pointer    :: cpoly
      type(sitetype)                     , pointer    :: csite
      type(patchtype)                    , pointer    :: cpatch
      real             , dimension(nzg)               :: cgrid_dmean_soil_hcap
      integer                                         :: ipy
      integer                                         :: isi
      integer                                         :: ipa
      integer                                         :: ico
      integer                                         :: k
      integer                                         :: nsoil
      real                                            :: poly_lai
      real                                            :: poly_wai
      real                                            :: poly_nplant
      real                                            :: can_exner
      real                                            :: atm_exner
      real                                            :: daylight_i
      real                                            :: site_area_i
      real                                            :: poly_area_i
      real                                            :: site_wgt
      real                                            :: patch_wgt
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop over polygons.                                                            !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)


         !----- Inverse of this polygon area (it should be always 1.) ---------------------!
         poly_area_i = 1./sum(cpoly%area)
         !---------------------------------------------------------------------------------!


         !----- Re-set some support variables. --------------------------------------------!
         poly_lai                 = 0.0
         poly_wai                 = 0.0
         poly_nplant              = 0.0
         cgrid_dmean_soil_hcap(:) = 0.0
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Loop over sites.                                                            !
         !---------------------------------------------------------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            !----- Inverse of this site area (it should be always 1.) ---------------------!
            site_area_i = 1./sum(csite%area)
            !------------------------------------------------------------------------------!


            !----- Site weight. -----------------------------------------------------------!
            site_wgt = cpoly%area(isi) * poly_area_i
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Find the average day length.  In case of polar night, we set the factor  !
            ! to zero so it won't become a singularity.                                    !
            !------------------------------------------------------------------------------!
            if (cpoly%daylight(isi) >= dtlsm) then
               daylight_i = 1. / cpoly%daylight(isi)
            else
               daylight_i = 0.
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Find the derived properties for the air above canopy.                   !
            !------------------------------------------------------------------------------!
            atm_exner                 = press2exner (cpoly%dmean_atm_prss(isi))
            cpoly%dmean_atm_temp(isi) = extheta2temp(atm_exner,cpoly%dmean_atm_theta(isi))
            cpoly%dmean_atm_rhos(isi) = idealdenssh ( cpoly%dmean_atm_prss  (isi)          &
                                                    , cpoly%dmean_atm_temp  (isi)          &
                                                    , cpoly%dmean_atm_shv   (isi) )
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Loop over patches.                                                       !
            !------------------------------------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)


               !----- Site weight. --------------------------------------------------------!
               patch_wgt = csite%area(ipa) * site_area_i * site_wgt
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !      Normalise daytime only variables.                                    !
               !---------------------------------------------------------------------------!
               csite%dmean_albedo     (ipa)  = csite%dmean_albedo     (ipa) * daylight_i
               csite%dmean_albedo_par (ipa)  = csite%dmean_albedo_par (ipa) * daylight_i
               csite%dmean_albedo_nir (ipa)  = csite%dmean_albedo_nir (ipa) * daylight_i
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !      Now we find the derived properties for the canopy air space.         !
               !---------------------------------------------------------------------------!
               can_exner                 = press2exner ( csite%dmean_can_prss  (ipa) )
               csite%dmean_can_temp(ipa) = extheta2temp( can_exner                         &
                                                       , csite%dmean_can_theta (ipa) )
               csite%dmean_can_rhos(ipa) = idealdenssh ( csite%dmean_can_prss  (ipa)       &
                                                       , csite%dmean_can_temp  (ipa)       &
                                                       , csite%dmean_can_shv   (ipa)       )
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !      Aggregate the polygon-level variables of those variables that were   !
               ! not integrated during the day.                                            !
               !---------------------------------------------------------------------------!
               cgrid%dmean_co2_residual   (ipy) = cgrid%dmean_co2_residual   (ipy)         &
                                                + csite%dmean_co2_residual   (ipa)         &
                                                * patch_wgt
               cgrid%dmean_energy_residual(ipy) = cgrid%dmean_energy_residual(ipy)         &
                                                + csite%dmean_energy_residual(ipa)         &
                                                * patch_wgt
               cgrid%dmean_water_residual (ipy) = cgrid%dmean_water_residual (ipy)         &
                                                + csite%dmean_water_residual (ipa)         &
                                                * patch_wgt
               cgrid%dmean_albedo         (ipy) = cgrid%dmean_albedo         (ipy)         &
                                                + csite%dmean_albedo         (ipa)         &
                                                * patch_wgt
               cgrid%dmean_albedo_par     (ipy) = cgrid%dmean_albedo_par     (ipy)         &
                                                + csite%dmean_albedo_par     (ipa)         &
                                                * patch_wgt
               cgrid%dmean_albedo_nir     (ipy) = cgrid%dmean_albedo_nir     (ipy)         &
                                                + csite%dmean_albedo_nir     (ipa)         &
                                                * patch_wgt
               cgrid%dmean_A_decomp       (ipy) = cgrid%dmean_A_decomp       (ipy)         &
                                                + csite%dmean_A_decomp       (ipa)         &
                                                * patch_wgt
               cgrid%dmean_Af_decomp      (ipy) = cgrid%dmean_Af_decomp      (ipy)         &
                                                + csite%dmean_Af_decomp      (ipa)         &
                                                * patch_wgt
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !     Soil matric potential, temperature, and liquid water.                 !
               !---------------------------------------------------------------------------!
               do k=1,nzg
                  nsoil = cpoly%ntext_soil(k,isi)
                  call uextcm2tl( csite%dmean_soil_energy(k,ipa)                           &
                                , csite%dmean_soil_water (k,ipa) * wdns                    &
                                , soil(nsoil)%slcpd                                        &
                                , csite%dmean_soil_temp  (k,ipa)                           &
                                , csite%dmean_soil_fliq  (k,ipa))
                  cgrid_dmean_soil_hcap   (k)     = cgrid_dmean_soil_hcap(k)               &
                                                  + soil(nsoil)%slcpd * patch_wgt
               end do
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !   If the patch had some temporary snow/pounding layer, convert the mean   !
               ! energy to J/kg, then find the mean temperature and liquid fraction.       !
               ! Otherwise, set them to either zero or default values.                     !
               !---------------------------------------------------------------------------!
               if (csite%dmean_sfcw_mass(ipa) > tiny_sfcwater_mass) then
                  csite%dmean_sfcw_energy(ipa) = csite%dmean_sfcw_energy(ipa)              &
                                               / csite%dmean_sfcw_mass  (ipa)
                  call uint2tl( csite%dmean_sfcw_energy(ipa), csite%dmean_sfcw_temp(ipa)   &
                              , csite%dmean_sfcw_fliq  (ipa))
               else
                  csite%dmean_sfcw_mass  (ipa)  = 0.
                  csite%dmean_sfcw_depth (ipa)  = 0.
                  csite%dmean_sfcw_energy(ipa)  = 0.
                  csite%dmean_sfcw_temp  (ipa)  = csite%dmean_soil_temp(nzg,ipa)
                  csite%dmean_sfcw_fliq  (ipa)  = csite%dmean_soil_fliq(nzg,ipa)
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Loop over the cohorts.                                               !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts
                  !------------------------------------------------------------------------!
                  !     Aggregate AIs and density, they may be used to normalise averages. !
                  !------------------------------------------------------------------------!
                  poly_nplant = poly_nplant + cpatch%nplant(ico) * patch_wgt
                  poly_lai    = poly_lai    + cpatch%lai   (ico) * patch_wgt
                  poly_wai    = poly_wai    + cpatch%wai   (ico) * patch_wgt
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Normalise daytime only variables.                                 !
                  !------------------------------------------------------------------------!
                  cpatch%dmean_light_level     (ico) = cpatch%dmean_light_level     (ico)  &
                                                     * daylight_i
                  cpatch%dmean_light_level_beam(ico) = cpatch%dmean_light_level_beam(ico)  &
                                                     * daylight_i
                  cpatch%dmean_light_level_diff(ico) = cpatch%dmean_light_level_diff(ico)  &
                                                     * daylight_i

                  cpatch%dmean_par_level_beam(ico)   = cpatch%dmean_par_level_beam  (ico)  &
                                                     * daylight_i
                  cpatch%dmean_par_level_diffu(ico)   = cpatch%dmean_par_level_diffu  (ico)  &
                                                     * daylight_i
                  cpatch%dmean_par_level_diffd(ico)   = cpatch%dmean_par_level_diffd  (ico)  &
                                                     * daylight_i

                  cpatch%dmean_fs_open         (ico) = cpatch%dmean_fs_open         (ico)  &
                                                     * daylight_i
                  cpatch%dmean_fsw             (ico) = cpatch%dmean_fsw             (ico)  &
                                                     * daylight_i
                  cpatch%dmean_fsn             (ico) = cpatch%dmean_fsn             (ico)  &
                                                     * daylight_i
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Aggregate the polygon-level variables of the variables that have  !
                  ! not been aggregated during the previous day.                           !
                  !------------------------------------------------------------------------!
                  cgrid%dmean_fs_open          (ipy) = cgrid%dmean_fs_open          (ipy)  &
                                                     + cpatch%dmean_fs_open         (ico)  &
                                                     * cpatch%lai                   (ico)  &
                                                     * patch_wgt
                  cgrid%dmean_fsw              (ipy) = cgrid%dmean_fsw              (ipy)  &
                                                     + cpatch%dmean_fsw             (ico)  &
                                                     * cpatch%lai                   (ico)  &
                                                     * patch_wgt
                  cgrid%dmean_fsn              (ipy) = cgrid%dmean_fsn              (ipy)  &
                                                     + cpatch%dmean_fsn             (ico)  &
                                                     * cpatch%lai                   (ico)  &
                                                     * patch_wgt
                  cgrid%dmean_nppleaf          (ipy) = cgrid%dmean_nppleaf          (ipy)  &
                                                     + cpatch%dmean_nppleaf         (ico)  &
                                                     * cpatch%nplant                (ico)  &
                                                     * patch_wgt
                  cgrid%dmean_nppfroot         (ipy) = cgrid%dmean_nppfroot         (ipy)  &
                                                     + cpatch%dmean_nppfroot        (ico)  &
                                                     * cpatch%nplant                (ico)  &
                                                     * patch_wgt
                  cgrid%dmean_nppsapwood       (ipy) = cgrid%dmean_nppsapwood       (ipy)  &
                                                     + cpatch%dmean_nppsapwood      (ico)  &
                                                     * cpatch%nplant                (ico)  &
                                                     * patch_wgt
                  cgrid%dmean_nppcroot         (ipy) = cgrid%dmean_nppcroot         (ipy)  &
                                                     + cpatch%dmean_nppcroot        (ico)  &
                                                     * cpatch%nplant                (ico)  &
                                                     * patch_wgt
                  cgrid%dmean_nppseeds         (ipy) = cgrid%dmean_nppseeds         (ipy)  &
                                                     + cpatch%dmean_nppseeds        (ico)  &
                                                     * cpatch%nplant                (ico)  &
                                                     * patch_wgt
                  cgrid%dmean_nppwood          (ipy) = cgrid%dmean_nppwood          (ipy)  &
                                                     + cpatch%dmean_nppwood         (ico)  &
                                                     * cpatch%nplant                (ico)  &
                                                     * patch_wgt
                  cgrid%dmean_nppdaily         (ipy) = cgrid%dmean_nppdaily         (ipy)  &
                                                     + cpatch%dmean_nppdaily        (ico)  &
                                                     * cpatch%nplant                (ico)  &
                                                     * patch_wgt
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Find the vegetation temperature and liquid fraction.               !
                  !------------------------------------------------------------------------!
                  !----- Leaf. ------------------------------------------------------------!
                  if (cpatch%dmean_leaf_hcap(ico) > 0.) then
                     call uextcm2tl( cpatch%dmean_leaf_energy(ico)                         &
                                   , cpatch%dmean_leaf_water (ico)                         &
                                   , cpatch%dmean_leaf_hcap  (ico)                         &
                                   , cpatch%dmean_leaf_temp  (ico)                         &
                                   , cpatch%dmean_leaf_fliq  (ico) )
                  else
                     cpatch%dmean_leaf_vpdef(ico) = csite%dmean_can_vpdef(ipa)
                     cpatch%dmean_leaf_temp (ico) = csite%dmean_can_temp (ipa)
                     if (csite%dmean_can_temp(ipa) > t00) then
                        cpatch%dmean_leaf_fliq(ico) = 1.0
                     elseif (csite%dmean_can_temp(ipa) == t00) then
                        cpatch%dmean_leaf_fliq(ico) = 0.5
                     else
                        cpatch%dmean_leaf_fliq(ico) = 0.0
                     end if
                  end if
                  !----- Wood. ------------------------------------------------------------!
                  if (cpatch%dmean_wood_hcap(ico) > 0.) then
                     call uextcm2tl( cpatch%dmean_wood_energy(ico)                         &
                                   , cpatch%dmean_wood_water (ico)                         &
                                   , cpatch%dmean_wood_hcap  (ico)                         &
                                   , cpatch%dmean_wood_temp  (ico)                         &
                                   , cpatch%dmean_wood_fliq  (ico) )
                  else
                     cpatch%dmean_wood_temp(ico) = csite%dmean_can_temp(ipa)
                     if (csite%dmean_can_temp(ipa) > t00) then
                        cpatch%dmean_wood_fliq(ico) = 1.0
                     elseif (csite%dmean_can_temp(ipa) == t00) then
                        cpatch%dmean_wood_fliq(ico) = 0.5
                     else
                        cpatch%dmean_wood_fliq(ico) = 0.0
                     end if
                  end if
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!






         !---------------------------------------------------------------------------------!
         !      Find the derived properties for the air above canopy.                      !
         !---------------------------------------------------------------------------------!
         atm_exner                 = press2exner (cgrid%dmean_atm_prss(ipy))
         cgrid%dmean_atm_temp(ipy) = extheta2temp(atm_exner,cgrid%dmean_atm_theta(ipy))
         cgrid%dmean_atm_rhos(ipy) = idealdenssh ( cgrid%dmean_atm_prss  (ipy)             &
                                                 , cgrid%dmean_atm_temp  (ipy)             &
                                                 , cgrid%dmean_atm_shv   (ipy) )
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !      Find the derived properties for the canopy air space.                      !
         !---------------------------------------------------------------------------------!
         can_exner                 = press2exner (cgrid%dmean_can_prss(ipy))
         cgrid%dmean_can_temp(ipy) = extheta2temp(can_exner,cgrid%dmean_can_theta(ipy))
         cgrid%dmean_can_rhos(ipy) = idealdenssh ( cgrid%dmean_can_prss  (ipy)             &
                                                 , cgrid%dmean_can_temp  (ipy)             &
                                                 , cgrid%dmean_can_shv   (ipy) )
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !   If the patch had some temporary snow/pounding layer, convert the mean energy  !
         ! to J/kg, then find the mean temperature and liquid fraction.  Otherwise, set    !
         ! them to either zero or default values.                                          !
         !---------------------------------------------------------------------------------!
         if (cgrid%dmean_sfcw_mass(ipy) > tiny_sfcwater_mass) then
            cgrid%dmean_sfcw_energy(ipy) = cgrid%dmean_sfcw_energy(ipy)                    &
                                         / cgrid%dmean_sfcw_mass(ipy)
            call uint2tl(cgrid%dmean_sfcw_energy(ipy),cgrid%dmean_sfcw_temp(ipy)           &
                        ,cgrid%dmean_sfcw_fliq(ipy))
         else
            cgrid%dmean_sfcw_mass  (ipy)  = 0.
            cgrid%dmean_sfcw_depth (ipy)  = 0.
            cgrid%dmean_sfcw_energy(ipy)  = 0.
            cgrid%dmean_sfcw_temp  (ipy)  = cgrid%dmean_soil_temp(nzg,ipy)
            cgrid%dmean_sfcw_fliq  (ipy)  = cgrid%dmean_soil_fliq(nzg,ipy)
         end if
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Find the temperature and the fraction of liquid water.                      !
         !---------------------------------------------------------------------------------!
         do k=1,nzg
            call uextcm2tl( cgrid%dmean_soil_energy(k,ipy)                                 &
                          , cgrid%dmean_soil_water (k,ipy) * wdns                          &
                          , cgrid_dmean_soil_hcap  (k)                                     &
                          , cgrid%dmean_soil_temp  (k,ipy)                                 &
                          , cgrid%dmean_soil_fliq  (k,ipy))
         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the vegetation temperature and liquid fraction.                        !
         !---------------------------------------------------------------------------------!
         !----- Leaf. ---------------------------------------------------------------------!
         if (cgrid%dmean_leaf_hcap(ipy) > 0.) then
            call uextcm2tl( cgrid%dmean_leaf_energy(ipy), cgrid%dmean_leaf_water (ipy)     &
                          , cgrid%dmean_leaf_hcap  (ipy), cgrid%dmean_leaf_temp  (ipy)     &
                          , cgrid%dmean_leaf_fliq  (ipy) )
         else
            cgrid%dmean_leaf_temp (ipy) = cgrid%dmean_can_temp (ipy)
            if (cgrid%dmean_can_temp(ipy) > t00) then
               cgrid%dmean_leaf_fliq(ipy) = 1.0
            elseif (cgrid%dmean_can_temp(ipy) == t00) then
               cgrid%dmean_leaf_fliq(ipy) = 0.5
            else
               cgrid%dmean_leaf_fliq(ipy) = 0.0
            end if
         end if
         !----- Wood. ---------------------------------------------------------------------!
         if (cgrid%dmean_wood_hcap(ipy) > 0.) then
            call uextcm2tl( cgrid%dmean_wood_energy(ipy)                                   &
                          , cgrid%dmean_wood_water (ipy)                                   &
                          , cgrid%dmean_wood_hcap  (ipy)                                   &
                          , cgrid%dmean_wood_temp  (ipy)                                   &
                          , cgrid%dmean_wood_fliq  (ipy) )
         else
            cgrid%dmean_wood_temp(ipy) = cgrid%dmean_can_temp(ipy)
            if (cgrid%dmean_can_temp(ipy) > t00) then
               cgrid%dmean_wood_fliq(ipy) = 1.0
            elseif (cgrid%dmean_can_temp(ipy) == t00) then
               cgrid%dmean_wood_fliq(ipy) = 0.5
            else
               cgrid%dmean_wood_fliq(ipy) = 0.0
            end if
         end if
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !    Normalise the "intensive" properties.  The weight was either the LAI, WAI,   !
         ! or plant density.  In case none of the cohorts qualified to contribute, then we !
         ! assign either the canopy air space property, or a default number.               !
         !---------------------------------------------------------------------------------!
         if (poly_lai > 0.) then
            cgrid%dmean_fs_open   (ipy) = cgrid%dmean_fs_open    (ipy) / poly_lai
            cgrid%dmean_fsw       (ipy) = cgrid%dmean_fsw        (ipy) / poly_lai
            cgrid%dmean_fsn       (ipy) = cgrid%dmean_fsn        (ipy) / poly_lai
         else
            cgrid%dmean_fs_open   (ipy) = 0.5
            cgrid%dmean_fsw       (ipy) = 0.5
            cgrid%dmean_fsn       (ipy) = 0.5
         end if
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
    end subroutine normalize_ed_dmean_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine resets the daily_averages for variables actually used in the       !
   ! integration.                                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine zero_ed_today_vars(cgrid)
      use ed_state_vars, only : edtype       & ! structure
                              , polygontype  & ! structure
                              , sitetype     & ! structure
                              , patchtype    ! ! structure
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target  :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer :: cpoly
      type(sitetype)   , pointer :: csite
      type(patchtype)  , pointer :: cpatch
      integer                    :: ipy
      integer                    :: isi
      integer                    :: ipa
      integer                    :: ico
      !------------------------------------------------------------------------------------!
      do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
               
         do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)
               
               !----- Reset variables stored in sitetype. ---------------------------------!
               csite%today_A_decomp(ipa)  = 0.0
               csite%today_Af_decomp(ipa) = 0.0
               !---------------------------------------------------------------------------!


               !----- Reset variables stored in patchtype. --------------------------------!
               do ico = 1, cpatch%ncohorts
                  cpatch%today_gpp          (ico) = 0.0
                  cpatch%today_nppleaf      (ico) = 0.0
                  cpatch%today_nppfroot     (ico) = 0.0
                  cpatch%today_nppsapwood   (ico) = 0.0
                  cpatch%today_nppcroot     (ico) = 0.0
                  cpatch%today_nppseeds     (ico) = 0.0
                  cpatch%today_nppwood      (ico) = 0.0
                  cpatch%today_nppdaily     (ico) = 0.0
                  cpatch%today_gpp_pot      (ico) = 0.0
                  cpatch%today_gpp_lightmax (ico) = 0.0
                  cpatch%today_gpp_moistmax (ico) = 0.0
                  cpatch%today_leaf_resp    (ico) = 0.0
                  cpatch%today_root_resp    (ico) = 0.0
               end do
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!
      return
   end subroutine zero_ed_today_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine resets the daily averages once the daily average file has been     !
   ! written and used to compute the monthly mean (in case the latter was requested).      !
   !---------------------------------------------------------------------------------------!
   subroutine zero_ed_dmean_vars(cgrid)
      use ed_state_vars, only : edtype        & ! structure
                              , polygontype   & ! structure
                              , sitetype      & ! structure
                              , patchtype     ! ! structure
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target  :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer :: cpoly
      type(sitetype)   , pointer :: csite
      type(patchtype)  , pointer :: cpatch
      integer                    :: ipy
      integer                    :: isi
      integer                    :: ipa
      integer                    :: ico
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !       Loop over polygons.                                                          !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)
         
         !----- Variables stored in edtype. -----------------------------------------------!
         cgrid%dmean_nppleaf            (ipy) = 0.0
         cgrid%dmean_nppfroot           (ipy) = 0.0
         cgrid%dmean_nppsapwood         (ipy) = 0.0
         cgrid%dmean_nppcroot           (ipy) = 0.0
         cgrid%dmean_nppseeds           (ipy) = 0.0
         cgrid%dmean_nppwood            (ipy) = 0.0
         cgrid%dmean_nppdaily           (ipy) = 0.0
         cgrid%dmean_A_decomp           (ipy) = 0.0
         cgrid%dmean_Af_decomp          (ipy) = 0.0
         cgrid%dmean_co2_residual       (ipy) = 0.0
         cgrid%dmean_energy_residual    (ipy) = 0.0
         cgrid%dmean_water_residual     (ipy) = 0.0
         cgrid%dmean_gpp                (ipy) = 0.0
         cgrid%dmean_npp                (ipy) = 0.0
         cgrid%dmean_leaf_resp          (ipy) = 0.0
         cgrid%dmean_root_resp          (ipy) = 0.0
         cgrid%dmean_growth_resp        (ipy) = 0.0
         cgrid%dmean_storage_resp       (ipy) = 0.0
         cgrid%dmean_vleaf_resp         (ipy) = 0.0
         cgrid%dmean_plresp             (ipy) = 0.0
         cgrid%dmean_leaf_energy        (ipy) = 0.0
         cgrid%dmean_leaf_water         (ipy) = 0.0
         cgrid%dmean_leaf_hcap          (ipy) = 0.0
         cgrid%dmean_leaf_vpdef         (ipy) = 0.0
         cgrid%dmean_leaf_temp          (ipy) = 0.0
         cgrid%dmean_leaf_fliq          (ipy) = 0.0
         cgrid%dmean_leaf_gsw           (ipy) = 0.0
         cgrid%dmean_leaf_gbw           (ipy) = 0.0
         cgrid%dmean_wood_energy        (ipy) = 0.0
         cgrid%dmean_wood_water         (ipy) = 0.0
         cgrid%dmean_wood_hcap          (ipy) = 0.0
         cgrid%dmean_wood_temp          (ipy) = 0.0
         cgrid%dmean_wood_fliq          (ipy) = 0.0
         cgrid%dmean_wood_gbw           (ipy) = 0.0
         cgrid%dmean_fs_open            (ipy) = 0.0
         cgrid%dmean_fsw                (ipy) = 0.0
         cgrid%dmean_fsn                (ipy) = 0.0
         cgrid%dmean_a_light            (ipy) = 0.0
         cgrid%dmean_a_rubp             (ipy) = 0.0
         cgrid%dmean_a_co2              (ipy) = 0.0
         cgrid%dmean_psi_open           (ipy) = 0.0
         cgrid%dmean_psi_closed         (ipy) = 0.0
         cgrid%dmean_water_supply       (ipy) = 0.0
         cgrid%dmean_par_l              (ipy) = 0.0
         cgrid%dmean_par_l_beam         (ipy) = 0.0
         cgrid%dmean_par_l_diff         (ipy) = 0.0
         cgrid%dmean_rshort_l           (ipy) = 0.0
         cgrid%dmean_rlong_l            (ipy) = 0.0
         cgrid%dmean_sensible_lc        (ipy) = 0.0
         cgrid%dmean_vapor_lc           (ipy) = 0.0
         cgrid%dmean_transp             (ipy) = 0.0
         cgrid%dmean_intercepted_al     (ipy) = 0.0
         cgrid%dmean_wshed_lg           (ipy) = 0.0
         cgrid%dmean_rshort_w           (ipy) = 0.0
         cgrid%dmean_rlong_w            (ipy) = 0.0
         cgrid%dmean_sensible_wc        (ipy) = 0.0
         cgrid%dmean_vapor_wc           (ipy) = 0.0
         cgrid%dmean_intercepted_aw     (ipy) = 0.0
         cgrid%dmean_wshed_wg           (ipy) = 0.0
         cgrid%dmean_rh                 (ipy) = 0.0
         cgrid%dmean_cwd_rh             (ipy) = 0.0
         cgrid%dmean_nep                (ipy) = 0.0
         cgrid%dmean_rk4step            (ipy) = 0.0
         cgrid%dmean_available_water    (ipy) = 0.0
         cgrid%dmean_can_theiv          (ipy) = 0.0
         cgrid%dmean_can_theta          (ipy) = 0.0
         cgrid%dmean_can_vpdef          (ipy) = 0.0
         cgrid%dmean_can_temp           (ipy) = 0.0
         cgrid%dmean_can_shv            (ipy) = 0.0
         cgrid%dmean_can_co2            (ipy) = 0.0
         cgrid%dmean_can_rhos           (ipy) = 0.0
         cgrid%dmean_can_prss           (ipy) = 0.0
         cgrid%dmean_gnd_temp           (ipy) = 0.0
         cgrid%dmean_gnd_shv            (ipy) = 0.0
         cgrid%dmean_can_ggnd           (ipy) = 0.0
         cgrid%dmean_sfcw_depth         (ipy) = 0.0
         cgrid%dmean_sfcw_energy        (ipy) = 0.0
         cgrid%dmean_sfcw_mass          (ipy) = 0.0
         cgrid%dmean_sfcw_temp          (ipy) = 0.0
         cgrid%dmean_sfcw_fliq          (ipy) = 0.0
         cgrid%dmean_soil_energy      (:,ipy) = 0.0
         cgrid%dmean_soil_mstpot      (:,ipy) = 0.0
         cgrid%dmean_soil_water       (:,ipy) = 0.0
         cgrid%dmean_soil_temp        (:,ipy) = 0.0
         cgrid%dmean_soil_fliq        (:,ipy) = 0.0
         cgrid%dmean_rshort_gnd         (ipy) = 0.0
         cgrid%dmean_par_gnd            (ipy) = 0.0
         cgrid%dmean_rlong_gnd          (ipy) = 0.0
         cgrid%dmean_rlongup            (ipy) = 0.0
         cgrid%dmean_parup              (ipy) = 0.0
         cgrid%dmean_nirup              (ipy) = 0.0
         cgrid%dmean_rshortup           (ipy) = 0.0
         cgrid%dmean_rnet               (ipy) = 0.0
         cgrid%dmean_albedo             (ipy) = 0.0
         cgrid%dmean_albedo_par         (ipy) = 0.0
         cgrid%dmean_albedo_nir         (ipy) = 0.0
         cgrid%dmean_rlong_albedo       (ipy) = 0.0
         cgrid%dmean_ustar              (ipy) = 0.0
         cgrid%dmean_tstar              (ipy) = 0.0
         cgrid%dmean_qstar              (ipy) = 0.0
         cgrid%dmean_cstar              (ipy) = 0.0
         cgrid%dmean_carbon_ac          (ipy) = 0.0
         cgrid%dmean_carbon_st          (ipy) = 0.0
         cgrid%dmean_vapor_gc           (ipy) = 0.0
         cgrid%dmean_vapor_ac           (ipy) = 0.0
         cgrid%dmean_smoist_gg        (:,ipy) = 0.0
         cgrid%dmean_throughfall        (ipy) = 0.0
         cgrid%dmean_transloss        (:,ipy) = 0.0
         cgrid%dmean_runoff             (ipy) = 0.0
         cgrid%dmean_drainage           (ipy) = 0.0
         cgrid%dmean_sensible_gc        (ipy) = 0.0
         cgrid%dmean_sensible_ac        (ipy) = 0.0
         cgrid%dmean_sensible_gg      (:,ipy) = 0.0
         cgrid%dmean_qthroughfall       (ipy) = 0.0
         cgrid%dmean_qrunoff            (ipy) = 0.0
         cgrid%dmean_qdrainage          (ipy) = 0.0
         cgrid%dmean_atm_theiv          (ipy) = 0.0
         cgrid%dmean_atm_theta          (ipy) = 0.0
         cgrid%dmean_atm_temp           (ipy) = 0.0
         cgrid%dmean_atm_vpdef          (ipy) = 0.0
         cgrid%dmean_atm_shv            (ipy) = 0.0
         cgrid%dmean_atm_rshort         (ipy) = 0.0
         cgrid%dmean_atm_rshort_diff    (ipy) = 0.0
         cgrid%dmean_atm_par            (ipy) = 0.0
         cgrid%dmean_atm_par_diff       (ipy) = 0.0
         cgrid%dmean_atm_rlong          (ipy) = 0.0
         cgrid%dmean_atm_vels           (ipy) = 0.0
         cgrid%dmean_atm_rhos           (ipy) = 0.0
         cgrid%dmean_atm_prss           (ipy) = 0.0
         cgrid%dmean_atm_co2            (ipy) = 0.0
         cgrid%dmean_pcpg               (ipy) = 0.0
         cgrid%dmean_qpcpg              (ipy) = 0.0
         cgrid%dmean_dpcpg              (ipy) = 0.0

         !---------------------------------------------------------------------------------!
         !       Loop over sites.                                                          !
         !---------------------------------------------------------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            cpoly%daylight             (isi) = 0.0
            cpoly%dmean_atm_theiv      (isi) = 0.0
            cpoly%dmean_atm_theta      (isi) = 0.0
            cpoly%dmean_atm_temp       (isi) = 0.0
            cpoly%dmean_atm_vpdef      (isi) = 0.0
            cpoly%dmean_atm_shv        (isi) = 0.0
            cpoly%dmean_atm_rshort     (isi) = 0.0
            cpoly%dmean_atm_rshort_diff(isi) = 0.0
            cpoly%dmean_atm_par        (isi) = 0.0
            cpoly%dmean_atm_par_diff   (isi) = 0.0
            cpoly%dmean_atm_rlong      (isi) = 0.0
            cpoly%dmean_atm_vels       (isi) = 0.0
            cpoly%dmean_atm_rhos       (isi) = 0.0
            cpoly%dmean_atm_prss       (isi) = 0.0
            cpoly%dmean_atm_co2        (isi) = 0.0
            cpoly%dmean_pcpg           (isi) = 0.0
            cpoly%dmean_qpcpg          (isi) = 0.0
            cpoly%dmean_dpcpg          (isi) = 0.0

            !------------------------------------------------------------------------------!
            !       Loop over sites.                                                       !
            !------------------------------------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               csite%dmean_A_decomp         (ipa) = 0.0
               csite%dmean_Af_decomp        (ipa) = 0.0
               csite%dmean_co2_residual     (ipa) = 0.0
               csite%dmean_energy_residual  (ipa) = 0.0
               csite%dmean_water_residual   (ipa) = 0.0
               csite%dmean_rh               (ipa) = 0.0
               csite%dmean_cwd_rh           (ipa) = 0.0
               csite%dmean_nep              (ipa) = 0.0
               csite%dmean_rk4step          (ipa) = 0.0
               csite%dmean_available_water  (ipa) = 0.0
               csite%dmean_can_theiv        (ipa) = 0.0
               csite%dmean_can_theta        (ipa) = 0.0
               csite%dmean_can_vpdef        (ipa) = 0.0
               csite%dmean_can_temp         (ipa) = 0.0
               csite%dmean_can_shv          (ipa) = 0.0
               csite%dmean_can_co2          (ipa) = 0.0
               csite%dmean_can_rhos         (ipa) = 0.0
               csite%dmean_can_prss         (ipa) = 0.0
               csite%dmean_gnd_temp         (ipa) = 0.0
               csite%dmean_gnd_shv          (ipa) = 0.0
               csite%dmean_can_ggnd         (ipa) = 0.0
               csite%dmean_sfcw_depth       (ipa) = 0.0
               csite%dmean_sfcw_energy      (ipa) = 0.0
               csite%dmean_sfcw_mass        (ipa) = 0.0
               csite%dmean_sfcw_temp        (ipa) = 0.0
               csite%dmean_sfcw_fliq        (ipa) = 0.0
               csite%dmean_soil_energy    (:,ipa) = 0.0
               csite%dmean_soil_mstpot    (:,ipa) = 0.0
               csite%dmean_soil_water     (:,ipa) = 0.0
               csite%dmean_soil_temp      (:,ipa) = 0.0
               csite%dmean_soil_fliq      (:,ipa) = 0.0
               csite%dmean_rshort_gnd       (ipa) = 0.0
               csite%dmean_par_gnd          (ipa) = 0.0
               csite%dmean_rlong_gnd        (ipa) = 0.0
               csite%dmean_rlongup          (ipa) = 0.0
               csite%dmean_parup            (ipa) = 0.0
               csite%dmean_nirup            (ipa) = 0.0
               csite%dmean_rshortup         (ipa) = 0.0
               csite%dmean_rnet             (ipa) = 0.0
               csite%dmean_albedo           (ipa) = 0.0
               csite%dmean_albedo_par       (ipa) = 0.0
               csite%dmean_albedo_nir       (ipa) = 0.0
               csite%dmean_rlong_albedo     (ipa) = 0.0
               csite%dmean_ustar            (ipa) = 0.0
               csite%dmean_tstar            (ipa) = 0.0
               csite%dmean_qstar            (ipa) = 0.0
               csite%dmean_cstar            (ipa) = 0.0
               csite%dmean_carbon_ac        (ipa) = 0.0
               csite%dmean_carbon_st        (ipa) = 0.0
               csite%dmean_vapor_gc         (ipa) = 0.0
               csite%dmean_vapor_ac         (ipa) = 0.0
               csite%dmean_smoist_gg      (:,ipa) = 0.0
               csite%dmean_throughfall      (ipa) = 0.0
               csite%dmean_transloss      (:,ipa) = 0.0
               csite%dmean_runoff           (ipa) = 0.0
               csite%dmean_drainage         (ipa) = 0.0
               csite%dmean_sensible_gc      (ipa) = 0.0
               csite%dmean_sensible_ac      (ipa) = 0.0
               csite%dmean_sensible_gg    (:,ipa) = 0.0
               csite%dmean_qthroughfall     (ipa) = 0.0
               csite%dmean_qrunoff          (ipa) = 0.0
               csite%dmean_qdrainage        (ipa) = 0.0



               !---------------------------------------------------------------------------!
               !       Loop over cohorts.                                                  !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1, cpatch%ncohorts
                  cpatch%dmean_nppleaf           (ico) = 0.0
                  cpatch%dmean_nppfroot          (ico) = 0.0
                  cpatch%dmean_nppsapwood        (ico) = 0.0
                  cpatch%dmean_nppcroot          (ico) = 0.0
                  cpatch%dmean_nppseeds          (ico) = 0.0
                  cpatch%dmean_nppwood           (ico) = 0.0
                  cpatch%dmean_nppdaily          (ico) = 0.0
                  cpatch%dmean_gpp               (ico) = 0.0
                  cpatch%dmean_npp               (ico) = 0.0
                  cpatch%dmean_leaf_resp         (ico) = 0.0
                  cpatch%dmean_root_resp         (ico) = 0.0
                  cpatch%dmean_growth_resp       (ico) = 0.0
                  cpatch%dmean_storage_resp      (ico) = 0.0
                  cpatch%dmean_vleaf_resp        (ico) = 0.0
                  cpatch%dmean_plresp            (ico) = 0.0
                  cpatch%dmean_leaf_energy       (ico) = 0.0
                  cpatch%dmean_leaf_water        (ico) = 0.0
                  cpatch%dmean_leaf_hcap         (ico) = 0.0
                  cpatch%dmean_leaf_vpdef        (ico) = 0.0
                  cpatch%dmean_leaf_temp         (ico) = 0.0
                  cpatch%dmean_leaf_fliq         (ico) = 0.0
                  cpatch%dmean_leaf_gsw          (ico) = 0.0
                  cpatch%dmean_leaf_gbw          (ico) = 0.0
                  cpatch%dmean_wood_energy       (ico) = 0.0
                  cpatch%dmean_wood_water        (ico) = 0.0
                  cpatch%dmean_wood_hcap         (ico) = 0.0
                  cpatch%dmean_wood_temp         (ico) = 0.0
                  cpatch%dmean_wood_fliq         (ico) = 0.0
                  cpatch%dmean_wood_gbw          (ico) = 0.0
                  cpatch%dmean_fs_open           (ico) = 0.0
                  cpatch%dmean_fsw               (ico) = 0.0
                  cpatch%dmean_fsn               (ico) = 0.0
                  cpatch%dmean_a_light           (ico) = 0.0
                  cpatch%dmean_a_rubp            (ico) = 0.0
                  cpatch%dmean_a_co2             (ico) = 0.0
                  cpatch%dmean_psi_open          (ico) = 0.0
                  cpatch%dmean_psi_closed        (ico) = 0.0
                  cpatch%dmean_water_supply      (ico) = 0.0
                  cpatch%dmean_light_level       (ico) = 0.0
                  cpatch%dmean_light_level_beam  (ico) = 0.0
                  cpatch%dmean_light_level_diff  (ico) = 0.0

                  cpatch%dmean_par_level_beam    (ico) = 0.0
                  cpatch%dmean_par_level_diffd   (ico) = 0.0
                  cpatch%dmean_par_level_diffu   (ico) = 0.0

                  cpatch%dmean_par_l             (ico) = 0.0
                  cpatch%dmean_par_l_beam        (ico) = 0.0
                  cpatch%dmean_par_l_diff        (ico) = 0.0
                  cpatch%dmean_rshort_l          (ico) = 0.0
                  cpatch%dmean_rlong_l           (ico) = 0.0
                  cpatch%dmean_sensible_lc       (ico) = 0.0
                  cpatch%dmean_vapor_lc          (ico) = 0.0
                  cpatch%dmean_transp            (ico) = 0.0
                  cpatch%dmean_intercepted_al    (ico) = 0.0
                  cpatch%dmean_wshed_lg          (ico) = 0.0
                  cpatch%dmean_rshort_w          (ico) = 0.0
                  cpatch%dmean_rlong_w           (ico) = 0.0
                  cpatch%dmean_rad_profile     (:,ico) = 0.0
                  cpatch%dmean_sensible_wc       (ico) = 0.0
                  cpatch%dmean_vapor_wc          (ico) = 0.0
                  cpatch%dmean_intercepted_aw    (ico) = 0.0
                  cpatch%dmean_wshed_wg          (ico) = 0.0
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine zero_ed_dmean_vars
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !                            |---------------------------------|                        !
   !                            |** MONTHLY AVERAGE SUBROUTINES **|                        !
   !                            |---------------------------------|                        !
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine integrates most of the monthly mean variables.  This sub-routine   !
   ! is called after the dmean variables have been normalised, so we can take advantage of !
   ! their values.                                                                         !
   !    A few variables are _NOT_ integrated here: quantities such as temperature and      !
   ! liquid fraction are found after the monthly mean of the extensive quantities have     !
   ! been normalised. Also, polygon-level variables that were not integrated to the daily  !
   ! means are not integrated here, they are found after the monthly mean at the native    !
   ! level has been found.                                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine integrate_ed_mmean_vars(cgrid)
      use ed_state_vars, only : edtype        & ! structure
                              , polygontype   & ! structure
                              , sitetype      & ! structure
                              , patchtype     ! ! structure
      use ed_max_dims  , only : n_dbh         & ! intent(in)
                              , n_pft         ! ! intent(in)
      use consts_coms  , only : yr_day        ! ! intent(in)
      use ed_misc_coms , only : current_time  & ! intent(in)
                              , simtime       ! ! structure
      implicit none
      !----- Argument. --------------------------------------------------------------------!
      type(edtype)      , target    :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype) , pointer   :: cpoly
      type(sitetype)    , pointer   :: csite
      type(patchtype)   , pointer   :: cpatch
      type(simtime)                 :: daybefore
      integer                       :: ipy
      integer                       :: isi
      integer                       :: ipa
      integer                       :: ico
      real                          :: ndaysi
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find which day we have just integrated, we will use it to determine the right  !
      ! scaling factor.                                                                    !
      !------------------------------------------------------------------------------------!
      call yesterday_info(current_time,daybefore,ndaysi)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Loop over polygons.                                                           !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         !---------------------------------------------------------------------------------!
         !    The following variables don't have daily means, but the instantaneous value  !
         ! is fine because they are updated only once a day.                               !
         !---------------------------------------------------------------------------------!
         cgrid%mmean_lai             (:,:,ipy) = cgrid%mmean_lai              (:,:,ipy)    &
                                               + cgrid%lai                    (:,:,ipy)    &
                                               * ndaysi
         cgrid%mmean_bleaf           (:,:,ipy) = cgrid%mmean_bleaf            (:,:,ipy)    &
                                               + cgrid%bleaf                  (:,:,ipy)    &
                                               * ndaysi
         cgrid%mmean_broot           (:,:,ipy) = cgrid%mmean_broot            (:,:,ipy)    &
                                               + cgrid%broot                  (:,:,ipy)    &
                                               * ndaysi
         cgrid%mmean_bstorage        (:,:,ipy) = cgrid%mmean_bstorage         (:,:,ipy)    &
                                               + cgrid%bstorage               (:,:,ipy)    &
                                               * ndaysi
         cgrid%mmean_bleaf_n         (:,:,ipy) = cgrid%mmean_bleaf_n          (:,:,ipy)    &
                                               + cgrid%bleaf_n                (:,:,ipy)    &
                                               * ndaysi
         cgrid%mmean_broot_n         (:,:,ipy) = cgrid%mmean_broot_n          (:,:,ipy)    &
                                               + cgrid%broot_n                (:,:,ipy)    &
                                               * ndaysi
         cgrid%mmean_bstorage_n      (:,:,ipy) = cgrid%mmean_bstorage_n       (:,:,ipy)    &
                                               + cgrid%bstorage_n             (:,:,ipy)    &
                                               * ndaysi
         cgrid%mmean_leaf_maintenance(:,:,ipy) = cgrid%mmean_leaf_maintenance (:,:,ipy)    &
                                               + cgrid%leaf_maintenance       (:,:,ipy)    &
                                               * ndaysi
         cgrid%mmean_root_maintenance(:,:,ipy) = cgrid%mmean_root_maintenance (:,:,ipy)    &
                                               + cgrid%root_maintenance       (:,:,ipy)    &
                                               * ndaysi
         cgrid%mmean_leaf_drop       (:,:,ipy) = cgrid%mmean_leaf_drop        (:,:,ipy)    &
                                               + cgrid%leaf_drop              (:,:,ipy)    &
                                               * ndaysi
         cgrid%mmean_fast_soil_c         (ipy) = cgrid%mmean_fast_soil_c          (ipy)    &
                                               + cgrid%fast_soil_c                (ipy)    &
                                               * ndaysi
         cgrid%mmean_slow_soil_c         (ipy) = cgrid%mmean_slow_soil_c          (ipy)    &
                                               + cgrid%slow_soil_c                (ipy)    &
                                               * ndaysi
         cgrid%mmean_struct_soil_c       (ipy) = cgrid%mmean_struct_soil_c        (ipy)    &
                                               + cgrid%struct_soil_c              (ipy)    &
                                               * ndaysi
         cgrid%mmean_struct_soil_l       (ipy) = cgrid%mmean_struct_soil_l        (ipy)    &
                                               + cgrid%struct_soil_l              (ipy)    &
                                               * ndaysi
         cgrid%mmean_cwd_c               (ipy) = cgrid%mmean_cwd_c                (ipy)    &
                                               + cgrid%cwd_c                      (ipy)    &
                                               * ndaysi
         cgrid%mmean_fast_soil_n         (ipy) = cgrid%mmean_fast_soil_n          (ipy)    &
                                               + cgrid%fast_soil_n                (ipy)    &
                                               * ndaysi
         cgrid%mmean_mineral_soil_n      (ipy) = cgrid%mmean_mineral_soil_n       (ipy)    &
                                               + cgrid%mineral_soil_n             (ipy)    &
                                               * ndaysi
         cgrid%mmean_cwd_n               (ipy) = cgrid%mmean_cwd_n                (ipy)    &
                                               + cgrid%cwd_n                      (ipy)    &
                                               * ndaysi
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Integrate the other polygon-level variables.                                 !
         !---------------------------------------------------------------------------------!
         cgrid%mmean_gpp              (ipy) = cgrid%mmean_gpp              (ipy)           &
                                            + cgrid%dmean_gpp              (ipy)           &
                                            * ndaysi
         cgrid%mmean_npp              (ipy) = cgrid%mmean_npp              (ipy)           &
                                            + cgrid%dmean_npp              (ipy)           &
                                            * ndaysi
         cgrid%mmean_leaf_resp        (ipy) = cgrid%mmean_leaf_resp        (ipy)           &
                                            + cgrid%dmean_leaf_resp        (ipy)           &
                                            * ndaysi
         cgrid%mmean_root_resp        (ipy) = cgrid%mmean_root_resp        (ipy)           &
                                            + cgrid%dmean_root_resp        (ipy)           &
                                            * ndaysi
         cgrid%mmean_growth_resp      (ipy) = cgrid%mmean_growth_resp      (ipy)           &
                                            + cgrid%dmean_growth_resp      (ipy)           &
                                            * ndaysi
         cgrid%mmean_storage_resp     (ipy) = cgrid%mmean_storage_resp     (ipy)           &
                                            + cgrid%dmean_storage_resp     (ipy)           &
                                            * ndaysi
         cgrid%mmean_vleaf_resp       (ipy) = cgrid%mmean_vleaf_resp       (ipy)           &
                                            + cgrid%dmean_vleaf_resp       (ipy)           &
                                            * ndaysi
         cgrid%mmean_plresp           (ipy) = cgrid%mmean_plresp           (ipy)           &
                                            + cgrid%dmean_plresp           (ipy)           &
                                            * ndaysi
         cgrid%mmean_leaf_energy      (ipy) = cgrid%mmean_leaf_energy      (ipy)           &
                                            + cgrid%dmean_leaf_energy      (ipy)           &
                                            * ndaysi
         cgrid%mmean_leaf_water       (ipy) = cgrid%mmean_leaf_water       (ipy)           &
                                            + cgrid%dmean_leaf_water       (ipy)           &
                                            * ndaysi
         cgrid%mmean_leaf_hcap        (ipy) = cgrid%mmean_leaf_hcap        (ipy)           &
                                            + cgrid%dmean_leaf_hcap        (ipy)           &
                                            * ndaysi
         cgrid%mmean_leaf_vpdef       (ipy) = cgrid%mmean_leaf_vpdef       (ipy)           &
                                            + cgrid%dmean_leaf_vpdef       (ipy)           &
                                            * ndaysi
         cgrid%mmean_leaf_gsw         (ipy) = cgrid%mmean_leaf_gsw         (ipy)           &
                                            + cgrid%dmean_leaf_gsw         (ipy)           &
                                            * ndaysi
         cgrid%mmean_leaf_gbw         (ipy) = cgrid%mmean_leaf_gbw         (ipy)           &
                                            + cgrid%dmean_leaf_gbw         (ipy)           &
                                            * ndaysi
         cgrid%mmean_wood_energy      (ipy) = cgrid%mmean_wood_energy      (ipy)           &
                                            + cgrid%dmean_wood_energy      (ipy)           &
                                            * ndaysi
         cgrid%mmean_wood_water       (ipy) = cgrid%mmean_wood_water       (ipy)           &
                                            + cgrid%dmean_wood_water       (ipy)           &
                                            * ndaysi
         cgrid%mmean_wood_hcap        (ipy) = cgrid%mmean_wood_hcap        (ipy)           &
                                            + cgrid%dmean_wood_hcap        (ipy)           &
                                            * ndaysi
         cgrid%mmean_wood_gbw         (ipy) = cgrid%mmean_wood_gbw         (ipy)           &
                                            + cgrid%dmean_wood_gbw         (ipy)           &
                                            * ndaysi
         cgrid%mmean_fs_open          (ipy) = cgrid%mmean_fs_open          (ipy)           &
                                            + cgrid%dmean_fs_open          (ipy)           &
                                            * ndaysi
         cgrid%mmean_fsw              (ipy) = cgrid%mmean_fsw              (ipy)           &
                                            + cgrid%dmean_fsw              (ipy)           &
                                            * ndaysi
         cgrid%mmean_fsn              (ipy) = cgrid%mmean_fsn              (ipy)           &
                                            + cgrid%dmean_fsn              (ipy)           &
                                            * ndaysi
         cgrid%mmean_a_light          (ipy) = cgrid%mmean_a_light          (ipy)           &
                                            + cgrid%dmean_a_light          (ipy)           &
                                            * ndaysi
         cgrid%mmean_a_rubp           (ipy) = cgrid%mmean_a_rubp           (ipy)           &
                                            + cgrid%dmean_a_rubp           (ipy)           &
                                            * ndaysi
         cgrid%mmean_a_co2            (ipy) = cgrid%mmean_a_co2            (ipy)           &
                                            + cgrid%dmean_a_co2            (ipy)           &
                                            * ndaysi
         cgrid%mmean_psi_open         (ipy) = cgrid%mmean_psi_open         (ipy)           &
                                            + cgrid%dmean_psi_open         (ipy)           &
                                            * ndaysi
         cgrid%mmean_psi_closed       (ipy) = cgrid%mmean_psi_closed       (ipy)           &
                                            + cgrid%dmean_psi_closed       (ipy)           &
                                            * ndaysi
         cgrid%mmean_water_supply     (ipy) = cgrid%mmean_water_supply     (ipy)           &
                                            + cgrid%dmean_water_supply     (ipy)           &
                                            * ndaysi
         cgrid%mmean_par_l            (ipy) = cgrid%mmean_par_l            (ipy)           &
                                            + cgrid%dmean_par_l            (ipy)           &
                                            * ndaysi
         cgrid%mmean_par_l_beam       (ipy) = cgrid%mmean_par_l_beam       (ipy)           &
                                            + cgrid%dmean_par_l_beam       (ipy)           &
                                            * ndaysi
         cgrid%mmean_par_l_diff       (ipy) = cgrid%mmean_par_l_diff       (ipy)           &
                                            + cgrid%dmean_par_l_diff       (ipy)           &
                                            * ndaysi
         cgrid%mmean_rshort_l         (ipy) = cgrid%mmean_rshort_l         (ipy)           &
                                            + cgrid%dmean_rshort_l         (ipy)           &
                                            * ndaysi
         cgrid%mmean_rlong_l          (ipy) = cgrid%mmean_rlong_l          (ipy)           &
                                            + cgrid%dmean_rlong_l          (ipy)           &
                                            * ndaysi
         cgrid%mmean_sensible_lc      (ipy) = cgrid%mmean_sensible_lc      (ipy)           &
                                            + cgrid%dmean_sensible_lc      (ipy)           &
                                            * ndaysi
         cgrid%mmean_vapor_lc         (ipy) = cgrid%mmean_vapor_lc         (ipy)           &
                                            + cgrid%dmean_vapor_lc         (ipy)           &
                                            * ndaysi
         cgrid%mmean_transp           (ipy) = cgrid%mmean_transp           (ipy)           &
                                            + cgrid%dmean_transp           (ipy)           &
                                            * ndaysi
         cgrid%mmean_intercepted_al   (ipy) = cgrid%mmean_intercepted_al   (ipy)           &
                                            + cgrid%dmean_intercepted_al   (ipy)           &
                                            * ndaysi
         cgrid%mmean_wshed_lg         (ipy) = cgrid%mmean_wshed_lg         (ipy)           &
                                            + cgrid%dmean_wshed_lg         (ipy)           &
                                            * ndaysi
         cgrid%mmean_rshort_w         (ipy) = cgrid%mmean_rshort_w         (ipy)           &
                                            + cgrid%dmean_rshort_w         (ipy)           &
                                            * ndaysi
         cgrid%mmean_rlong_w          (ipy) = cgrid%mmean_rlong_w          (ipy)           &
                                            + cgrid%dmean_rlong_w          (ipy)           &
                                            * ndaysi
         cgrid%mmean_sensible_wc      (ipy) = cgrid%mmean_sensible_wc      (ipy)           &
                                            + cgrid%dmean_sensible_wc      (ipy)           &
                                            * ndaysi
         cgrid%mmean_vapor_wc         (ipy) = cgrid%mmean_vapor_wc         (ipy)           &
                                            + cgrid%dmean_vapor_wc         (ipy)           &
                                            * ndaysi
         cgrid%mmean_intercepted_aw   (ipy) = cgrid%mmean_intercepted_aw   (ipy)           &
                                            + cgrid%dmean_intercepted_aw   (ipy)           &
                                            * ndaysi
         cgrid%mmean_wshed_wg         (ipy) = cgrid%mmean_wshed_wg         (ipy)           &
                                            + cgrid%dmean_wshed_wg         (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppleaf          (ipy) = cgrid%mmean_nppleaf          (ipy)           &
                                            + cgrid%dmean_nppleaf          (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppfroot         (ipy) = cgrid%mmean_nppfroot         (ipy)           &
                                            + cgrid%dmean_nppfroot         (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppsapwood       (ipy) = cgrid%mmean_nppsapwood       (ipy)           &
                                            + cgrid%dmean_nppsapwood       (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppcroot         (ipy) = cgrid%mmean_nppcroot         (ipy)           &
                                            + cgrid%dmean_nppcroot         (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppseeds         (ipy) = cgrid%mmean_nppseeds         (ipy)           &
                                            + cgrid%dmean_nppseeds         (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppwood          (ipy) = cgrid%mmean_nppwood          (ipy)           &
                                            + cgrid%dmean_nppwood          (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppdaily         (ipy) = cgrid%mmean_nppdaily         (ipy)           &
                                            + cgrid%dmean_nppdaily         (ipy)           &
                                            * ndaysi
         cgrid%mmean_rh               (ipy) = cgrid%mmean_rh               (ipy)           &
                                            + cgrid%dmean_rh               (ipy)           &
                                            * ndaysi
         cgrid%mmean_cwd_rh           (ipy) = cgrid%mmean_cwd_rh           (ipy)           &
                                            + cgrid%dmean_cwd_rh           (ipy)           &
                                            * ndaysi
         cgrid%mmean_nep              (ipy) = cgrid%mmean_nep              (ipy)           &
                                            + cgrid%dmean_nep              (ipy)           &
                                            * ndaysi
         cgrid%mmean_rk4step          (ipy) = cgrid%mmean_rk4step          (ipy)           &
                                            + cgrid%dmean_rk4step          (ipy)           &
                                            * ndaysi
         cgrid%mmean_available_water  (ipy) = cgrid%mmean_available_water  (ipy)           &
                                            + cgrid%dmean_available_water  (ipy)           &
                                            * ndaysi
         cgrid%mmean_can_theiv        (ipy) = cgrid%mmean_can_theiv        (ipy)           &
                                            + cgrid%dmean_can_theiv        (ipy)           &
                                            * ndaysi
         cgrid%mmean_can_theta        (ipy) = cgrid%mmean_can_theta        (ipy)           &
                                            + cgrid%dmean_can_theta        (ipy)           &
                                            * ndaysi
         cgrid%mmean_can_vpdef        (ipy) = cgrid%mmean_can_vpdef        (ipy)           &
                                            + cgrid%dmean_can_vpdef        (ipy)           &
                                            * ndaysi
         cgrid%mmean_can_shv          (ipy) = cgrid%mmean_can_shv          (ipy)           &
                                            + cgrid%dmean_can_shv          (ipy)           &
                                            * ndaysi
         cgrid%mmean_can_co2          (ipy) = cgrid%mmean_can_co2          (ipy)           &
                                            + cgrid%dmean_can_co2          (ipy)           &
                                            * ndaysi
         cgrid%mmean_can_prss         (ipy) = cgrid%mmean_can_prss         (ipy)           &
                                            + cgrid%dmean_can_prss         (ipy)           &
                                            * ndaysi
         cgrid%mmean_gnd_temp         (ipy) = cgrid%mmean_gnd_temp         (ipy)           &
                                            + cgrid%dmean_gnd_temp         (ipy)           &
                                            * ndaysi
         cgrid%mmean_gnd_shv          (ipy) = cgrid%mmean_gnd_shv          (ipy)           &
                                            + cgrid%dmean_gnd_shv          (ipy)           &
                                            * ndaysi
         cgrid%mmean_can_ggnd         (ipy) = cgrid%mmean_can_ggnd         (ipy)           &
                                            + cgrid%dmean_can_ggnd         (ipy)           &
                                            * ndaysi
         cgrid%mmean_sfcw_depth       (ipy) = cgrid%mmean_sfcw_depth       (ipy)           &
                                            + cgrid%dmean_sfcw_depth       (ipy)           &
                                            * ndaysi
         !----- Temporarily convert energy to extensive [J/m2]. ---------------------------!
         cgrid%mmean_sfcw_energy      (ipy) = cgrid%mmean_sfcw_energy      (ipy)           &
                                            + cgrid%dmean_sfcw_energy      (ipy)           &
                                            * cgrid%dmean_sfcw_mass        (ipy)           &
                                            * ndaysi
         cgrid%mmean_sfcw_mass        (ipy) = cgrid%mmean_sfcw_mass        (ipy)           &
                                            + cgrid%dmean_sfcw_mass        (ipy)           &
                                            * ndaysi
         cgrid%mmean_soil_energy    (:,ipy) = cgrid%mmean_soil_energy    (:,ipy)           &
                                            + cgrid%dmean_soil_energy    (:,ipy)           &
                                            * ndaysi
         cgrid%mmean_soil_mstpot    (:,ipy) = cgrid%mmean_soil_mstpot    (:,ipy)           &
                                            + cgrid%dmean_soil_mstpot    (:,ipy)           &
                                            * ndaysi
         cgrid%mmean_soil_water     (:,ipy) = cgrid%mmean_soil_water     (:,ipy)           &
                                            + cgrid%dmean_soil_water     (:,ipy)           &
                                            * ndaysi
         cgrid%mmean_rshort_gnd       (ipy) = cgrid%mmean_rshort_gnd       (ipy)           &
                                            + cgrid%dmean_rshort_gnd       (ipy)           &
                                            * ndaysi
         cgrid%mmean_par_gnd          (ipy) = cgrid%mmean_par_gnd          (ipy)           &
                                            + cgrid%dmean_par_gnd          (ipy)           &
                                            * ndaysi
         cgrid%mmean_rlong_gnd        (ipy) = cgrid%mmean_rlong_gnd        (ipy)           &
                                            + cgrid%dmean_rlong_gnd        (ipy)           &
                                            * ndaysi
         cgrid%mmean_rlongup          (ipy) = cgrid%mmean_rlongup          (ipy)           &
                                            + cgrid%dmean_rlongup          (ipy)           &
                                            * ndaysi
         cgrid%mmean_parup            (ipy) = cgrid%mmean_parup            (ipy)           &
                                            + cgrid%dmean_parup            (ipy)           &
                                            * ndaysi
         cgrid%mmean_nirup            (ipy) = cgrid%mmean_nirup            (ipy)           &
                                            + cgrid%dmean_nirup            (ipy)           &
                                            * ndaysi
         cgrid%mmean_rshortup         (ipy) = cgrid%mmean_rshortup         (ipy)           &
                                            + cgrid%dmean_rshortup         (ipy)           &
                                            * ndaysi
         cgrid%mmean_rnet             (ipy) = cgrid%mmean_rnet             (ipy)           &
                                            + cgrid%dmean_rnet             (ipy)           &
                                            * ndaysi
         cgrid%mmean_albedo           (ipy) = cgrid%mmean_albedo           (ipy)           &
                                            + cgrid%dmean_albedo           (ipy)           &
                                            * ndaysi
         cgrid%mmean_albedo_par       (ipy) = cgrid%mmean_albedo_par       (ipy)           &
                                            + cgrid%dmean_albedo_par       (ipy)           &
                                            * ndaysi
         cgrid%mmean_albedo_nir       (ipy) = cgrid%mmean_albedo_nir       (ipy)           &
                                            + cgrid%dmean_albedo_nir       (ipy)           &
                                            * ndaysi
         cgrid%mmean_rlong_albedo     (ipy) = cgrid%mmean_rlong_albedo     (ipy)           &
                                            + cgrid%dmean_rlong_albedo     (ipy)           &
                                            * ndaysi
         cgrid%mmean_ustar            (ipy) = cgrid%mmean_ustar            (ipy)           &
                                            + cgrid%dmean_ustar            (ipy)           &
                                            * ndaysi
         cgrid%mmean_tstar            (ipy) = cgrid%mmean_tstar            (ipy)           &
                                            + cgrid%dmean_tstar            (ipy)           &
                                            * ndaysi
         cgrid%mmean_qstar            (ipy) = cgrid%mmean_qstar            (ipy)           &
                                            + cgrid%dmean_qstar            (ipy)           &
                                            * ndaysi
         cgrid%mmean_cstar            (ipy) = cgrid%mmean_cstar            (ipy)           &
                                            + cgrid%dmean_cstar            (ipy)           &
                                            * ndaysi
         cgrid%mmean_carbon_ac        (ipy) = cgrid%mmean_carbon_ac        (ipy)           &
                                            + cgrid%dmean_carbon_ac        (ipy)           &
                                            * ndaysi
         cgrid%mmean_carbon_st        (ipy) = cgrid%mmean_carbon_st        (ipy)           &
                                            + cgrid%dmean_carbon_st        (ipy)           &
                                            * ndaysi
         cgrid%mmean_vapor_gc         (ipy) = cgrid%mmean_vapor_gc         (ipy)           &
                                            + cgrid%dmean_vapor_gc         (ipy)           &
                                            * ndaysi
         cgrid%mmean_vapor_ac         (ipy) = cgrid%mmean_vapor_ac         (ipy)           &
                                            + cgrid%dmean_vapor_ac         (ipy)           &
                                            * ndaysi
         cgrid%mmean_smoist_gg      (:,ipy) = cgrid%mmean_smoist_gg      (:,ipy)           &
                                            + cgrid%dmean_smoist_gg      (:,ipy)           &
                                            * ndaysi
         cgrid%mmean_throughfall      (ipy) = cgrid%mmean_throughfall      (ipy)           &
                                            + cgrid%dmean_throughfall      (ipy)           &
                                            * ndaysi
         cgrid%mmean_transloss      (:,ipy) = cgrid%mmean_transloss      (:,ipy)           &
                                            + cgrid%dmean_transloss      (:,ipy)           &
                                            * ndaysi
         cgrid%mmean_runoff           (ipy) = cgrid%mmean_runoff           (ipy)           &
                                            + cgrid%dmean_runoff           (ipy)           &
                                            * ndaysi
         cgrid%mmean_drainage         (ipy) = cgrid%mmean_drainage         (ipy)           &
                                            + cgrid%dmean_drainage         (ipy)           &
                                            * ndaysi
         cgrid%mmean_sensible_gc      (ipy) = cgrid%mmean_sensible_gc      (ipy)           &
                                            + cgrid%dmean_sensible_gc      (ipy)           &
                                            * ndaysi
         cgrid%mmean_sensible_ac      (ipy) = cgrid%mmean_sensible_ac      (ipy)           &
                                            + cgrid%dmean_sensible_ac      (ipy)           &
                                            * ndaysi
         cgrid%mmean_sensible_gg    (:,ipy) = cgrid%mmean_sensible_gg    (:,ipy)           &
                                            + cgrid%dmean_sensible_gg    (:,ipy)           &
                                            * ndaysi
         cgrid%mmean_qthroughfall     (ipy) = cgrid%mmean_qthroughfall     (ipy)           &
                                            + cgrid%dmean_qthroughfall     (ipy)           &
                                            * ndaysi
         cgrid%mmean_qrunoff          (ipy) = cgrid%mmean_qrunoff          (ipy)           &
                                            + cgrid%dmean_qrunoff          (ipy)           &
                                            * ndaysi
         cgrid%mmean_qdrainage        (ipy) = cgrid%mmean_qdrainage        (ipy)           &
                                            + cgrid%dmean_qdrainage        (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppleaf          (ipy) = cgrid%mmean_nppleaf          (ipy)           &
                                            + cgrid%dmean_nppleaf          (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppfroot         (ipy) = cgrid%mmean_nppfroot         (ipy)           &
                                            + cgrid%dmean_nppfroot         (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppsapwood       (ipy) = cgrid%mmean_nppsapwood       (ipy)           &
                                            + cgrid%dmean_nppsapwood       (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppcroot         (ipy) = cgrid%mmean_nppcroot         (ipy)           &
                                            + cgrid%dmean_nppcroot         (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppseeds         (ipy) = cgrid%mmean_nppseeds         (ipy)           &
                                            + cgrid%dmean_nppseeds         (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppwood          (ipy) = cgrid%mmean_nppwood          (ipy)           &
                                            + cgrid%dmean_nppwood          (ipy)           &
                                            * ndaysi
         cgrid%mmean_nppdaily         (ipy) = cgrid%mmean_nppdaily         (ipy)           &
                                            + cgrid%dmean_nppdaily         (ipy)           &
                                            * ndaysi
         cgrid%mmean_A_decomp         (ipy) = cgrid%mmean_A_decomp         (ipy)           &
                                            + cgrid%dmean_A_decomp         (ipy)           &
                                            * ndaysi
         cgrid%mmean_Af_decomp        (ipy) = cgrid%mmean_Af_decomp        (ipy)           &
                                            + cgrid%dmean_Af_decomp        (ipy)           &
                                            * ndaysi
         cgrid%mmean_co2_residual     (ipy) = cgrid%mmean_co2_residual     (ipy)           &
                                            + cgrid%dmean_co2_residual     (ipy)           &
                                            * ndaysi
         cgrid%mmean_energy_residual  (ipy) = cgrid%mmean_energy_residual  (ipy)           &
                                            + cgrid%dmean_energy_residual  (ipy)           &
                                            * ndaysi
         cgrid%mmean_water_residual   (ipy) = cgrid%mmean_water_residual   (ipy)           &
                                            + cgrid%dmean_water_residual   (ipy)           &
                                            * ndaysi
         cgrid%mmean_atm_theiv        (ipy) = cgrid%mmean_atm_theiv        (ipy)           &
                                            + cgrid%dmean_atm_theiv        (ipy)           &
                                            * ndaysi
         cgrid%mmean_atm_theta        (ipy) = cgrid%mmean_atm_theta        (ipy)           &
                                            + cgrid%dmean_atm_theta        (ipy)           &
                                            * ndaysi
         cgrid%mmean_atm_vpdef        (ipy) = cgrid%mmean_atm_vpdef        (ipy)           &
                                            + cgrid%dmean_atm_vpdef        (ipy)           &
                                            * ndaysi
         cgrid%mmean_atm_shv          (ipy) = cgrid%mmean_atm_shv          (ipy)           &
                                            + cgrid%dmean_atm_shv          (ipy)           &
                                            * ndaysi
         cgrid%mmean_atm_rshort       (ipy) = cgrid%mmean_atm_rshort       (ipy)           &
                                            + cgrid%dmean_atm_rshort       (ipy)           &
                                            * ndaysi
         cgrid%mmean_atm_rshort_diff  (ipy) = cgrid%mmean_atm_rshort_diff  (ipy)           &
                                            + cgrid%dmean_atm_rshort_diff  (ipy)           &
                                            * ndaysi
         cgrid%mmean_atm_par          (ipy) = cgrid%mmean_atm_par          (ipy)           &
                                            + cgrid%dmean_atm_par          (ipy)           &
                                            * ndaysi
         cgrid%mmean_atm_par_diff     (ipy) = cgrid%mmean_atm_par_diff     (ipy)           &
                                            + cgrid%dmean_atm_par_diff     (ipy)           &
                                            * ndaysi
         cgrid%mmean_atm_rlong        (ipy) = cgrid%mmean_atm_rlong        (ipy)           &
                                            + cgrid%dmean_atm_rlong        (ipy)           &
                                            * ndaysi
         cgrid%mmean_atm_vels         (ipy) = cgrid%mmean_atm_vels         (ipy)           &
                                            + cgrid%dmean_atm_vels         (ipy)           &
                                            * ndaysi
         cgrid%mmean_atm_prss         (ipy) = cgrid%mmean_atm_prss         (ipy)           &
                                            + cgrid%dmean_atm_prss         (ipy)           &
                                            * ndaysi
         cgrid%mmean_atm_co2          (ipy) = cgrid%mmean_atm_co2          (ipy)           &
                                            + cgrid%dmean_atm_co2          (ipy)           &
                                            * ndaysi
         cgrid%mmean_pcpg             (ipy) = cgrid%mmean_pcpg             (ipy)           &
                                            + cgrid%dmean_pcpg             (ipy)           &
                                            * ndaysi
         cgrid%mmean_qpcpg            (ipy) = cgrid%mmean_qpcpg            (ipy)           &
                                            + cgrid%dmean_qpcpg            (ipy)           &
                                            * ndaysi
         cgrid%mmean_dpcpg            (ipy) = cgrid%mmean_dpcpg            (ipy)           &
                                            + cgrid%dmean_dpcpg            (ipy)           &
                                            * ndaysi
         !---------------------------------------------------------------------------------!
         !     Mean sum of squares.  Use double precision to integrating term, then        !
         ! convert the term back to single precision.  This step is needed to avoid under- !
         ! flows.                                                                          !
         !---------------------------------------------------------------------------------!
         cgrid%mmsqu_gpp              (ipy) = cgrid%mmsqu_gpp                 (ipy)        &
                                            + isqu_ftz(cgrid%dmean_gpp        (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_npp              (ipy) = cgrid%mmsqu_npp                 (ipy)        &
                                            + isqu_ftz(cgrid%dmean_npp        (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_plresp           (ipy) = cgrid%mmsqu_plresp              (ipy)        &
                                            + isqu_ftz(cgrid%dmean_plresp     (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_sensible_lc      (ipy) = cgrid%mmsqu_sensible_lc         (ipy)        &
                                            + isqu_ftz(cgrid%dmean_sensible_lc(ipy))       &
                                            * ndaysi
         cgrid%mmsqu_vapor_lc         (ipy) = cgrid%mmsqu_vapor_lc            (ipy)        &
                                            + isqu_ftz(cgrid%dmean_vapor_lc   (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_transp           (ipy) = cgrid%mmsqu_transp              (ipy)        &
                                            + isqu_ftz(cgrid%dmean_transp     (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_sensible_wc      (ipy) = cgrid%mmsqu_sensible_wc         (ipy)        &
                                            + isqu_ftz(cgrid%dmean_sensible_wc(ipy))       &
                                            * ndaysi
         cgrid%mmsqu_vapor_wc         (ipy) = cgrid%mmsqu_vapor_wc            (ipy)        &
                                            + isqu_ftz(cgrid%dmean_vapor_wc   (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_rh               (ipy) = cgrid%mmsqu_rh                  (ipy)        &
                                            + isqu_ftz(cgrid%dmean_rh         (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_cwd_rh           (ipy) = cgrid%mmsqu_cwd_rh              (ipy)        &
                                            + isqu_ftz(cgrid%dmean_cwd_rh     (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_nep              (ipy) = cgrid%mmsqu_nep                 (ipy)        &
                                            + isqu_ftz(cgrid%dmean_nep        (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_rlongup          (ipy) = cgrid%mmsqu_rlongup             (ipy)        &
                                            + isqu_ftz(cgrid%dmean_rlongup    (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_parup            (ipy) = cgrid%mmsqu_parup               (ipy)        &
                                            + isqu_ftz(cgrid%dmean_parup      (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_nirup            (ipy) = cgrid%mmsqu_nirup               (ipy)        &
                                            + isqu_ftz(cgrid%dmean_nirup      (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_rshortup         (ipy) = cgrid%mmsqu_rshortup            (ipy)        &
                                            + isqu_ftz(cgrid%dmean_rshortup   (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_rnet             (ipy) = cgrid%mmsqu_rnet                (ipy)        &
                                            + isqu_ftz(cgrid%dmean_rnet       (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_albedo           (ipy) = cgrid%mmsqu_albedo              (ipy)        &
                                            + isqu_ftz(cgrid%dmean_albedo     (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_ustar            (ipy) = cgrid%mmsqu_ustar               (ipy)        &
                                            + isqu_ftz(cgrid%dmean_ustar      (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_carbon_ac        (ipy) = cgrid%mmsqu_carbon_ac           (ipy)        &
                                            + isqu_ftz(cgrid%dmean_carbon_ac  (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_carbon_st        (ipy) = cgrid%mmsqu_carbon_st           (ipy)        &
                                            + isqu_ftz(cgrid%dmean_carbon_st  (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_vapor_gc         (ipy) = cgrid%mmsqu_vapor_gc            (ipy)        &
                                            + isqu_ftz(cgrid%dmean_vapor_gc   (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_vapor_ac         (ipy) = cgrid%mmsqu_vapor_ac            (ipy)        &
                                            + isqu_ftz(cgrid%dmean_vapor_ac   (ipy))       &
                                            * ndaysi
         cgrid%mmsqu_sensible_gc      (ipy) = cgrid%mmsqu_sensible_gc         (ipy)        &
                                            + isqu_ftz(cgrid%dmean_sensible_gc(ipy))       &
                                            * ndaysi
         cgrid%mmsqu_sensible_ac      (ipy) = cgrid%mmsqu_sensible_ac         (ipy)        &
                                            + isqu_ftz(cgrid%dmean_sensible_ac(ipy))       &
                                            * ndaysi
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Site loop.                                                                 !
         !---------------------------------------------------------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            !------------------------------------------------------------------------------!
            !    Integrate site-level variables.                                           !
            !------------------------------------------------------------------------------!
            cpoly%mmean_atm_theiv      (isi) = cpoly%mmean_atm_theiv      (isi)            &
                                             + cpoly%dmean_atm_theiv      (isi)            &
                                             * ndaysi
            cpoly%mmean_atm_theta      (isi) = cpoly%mmean_atm_theta      (isi)            &
                                             + cpoly%dmean_atm_theta      (isi)            &
                                             * ndaysi
            cpoly%mmean_atm_vpdef      (isi) = cpoly%mmean_atm_vpdef      (isi)            &
                                             + cpoly%dmean_atm_vpdef      (isi)            &
                                             * ndaysi
            cpoly%mmean_atm_shv        (isi) = cpoly%mmean_atm_shv        (isi)            &
                                             + cpoly%dmean_atm_shv        (isi)            &
                                             * ndaysi
            cpoly%mmean_atm_rshort     (isi) = cpoly%mmean_atm_rshort     (isi)            &
                                             + cpoly%dmean_atm_rshort     (isi)            &
                                             * ndaysi
            cpoly%mmean_atm_rshort_diff(isi) = cpoly%mmean_atm_rshort_diff(isi)            &
                                             + cpoly%dmean_atm_rshort_diff(isi)            &
                                             * ndaysi
            cpoly%mmean_atm_par        (isi) = cpoly%mmean_atm_par        (isi)            &
                                             + cpoly%dmean_atm_par        (isi)            &
                                             * ndaysi
            cpoly%mmean_atm_par_diff   (isi) = cpoly%mmean_atm_par_diff   (isi)            &
                                             + cpoly%dmean_atm_par_diff   (isi)            &
                                             * ndaysi
            cpoly%mmean_atm_rlong      (isi) = cpoly%mmean_atm_rlong      (isi)            &
                                             + cpoly%dmean_atm_rlong      (isi)            &
                                             * ndaysi
            cpoly%mmean_atm_vels       (isi) = cpoly%mmean_atm_vels       (isi)            &
                                             + cpoly%dmean_atm_vels       (isi)            &
                                             * ndaysi
            cpoly%mmean_atm_prss       (isi) = cpoly%mmean_atm_prss       (isi)            &
                                             + cpoly%dmean_atm_prss       (isi)            &
                                             * ndaysi
            cpoly%mmean_atm_co2        (isi) = cpoly%mmean_atm_co2        (isi)            &
                                             + cpoly%dmean_atm_co2        (isi)            &
                                             * ndaysi
            cpoly%mmean_pcpg           (isi) = cpoly%mmean_pcpg           (isi)            &
                                             + cpoly%dmean_pcpg           (isi)            &
                                             * ndaysi
            cpoly%mmean_qpcpg          (isi) = cpoly%mmean_qpcpg          (isi)            &
                                             + cpoly%dmean_qpcpg          (isi)            &
                                             * ndaysi
            cpoly%mmean_dpcpg          (isi) = cpoly%mmean_dpcpg          (isi)            &
                                             + cpoly%dmean_dpcpg          (isi)            &
                                             * ndaysi
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Patch loop.                                                             !
            !------------------------------------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               !---------------------------------------------------------------------------!
               !      Integrate the variables that don't have daily means because their    !
               ! time step is one day.                                                     !
               !---------------------------------------------------------------------------!
               csite%mmean_fast_soil_c      (ipa) = csite%mmean_fast_soil_c      (ipa)     &
                                                  + csite%fast_soil_c            (ipa)     &
                                                  * ndaysi
               csite%mmean_slow_soil_c      (ipa) = csite%mmean_slow_soil_c      (ipa)     &
                                                  + csite%slow_soil_c            (ipa)     &
                                                  * ndaysi
               csite%mmean_struct_soil_c    (ipa) = csite%mmean_struct_soil_c    (ipa)     &
                                                  + csite%structural_soil_c      (ipa)     &
                                                  * ndaysi
               csite%mmean_struct_soil_l    (ipa) = csite%mmean_struct_soil_l    (ipa)     &
                                                  + csite%structural_soil_l      (ipa)     &
                                                  * ndaysi
               csite%mmean_fast_soil_n      (ipa) = csite%mmean_fast_soil_n      (ipa)     &
                                                  + csite%fast_soil_n            (ipa)     &
                                                  * ndaysi
               csite%mmean_mineral_soil_n   (ipa) = csite%mmean_mineral_soil_n   (ipa)     &
                                                  + csite%mineralized_soil_n     (ipa)     &
                                                  * ndaysi
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !      Integrate patch-level variables.                                     !
               !---------------------------------------------------------------------------!
               csite%mmean_co2_residual     (ipa) = csite%mmean_co2_residual     (ipa)     &
                                                  + csite%dmean_co2_residual     (ipa)     &
                                                  * ndaysi
               csite%mmean_energy_residual  (ipa) = csite%mmean_energy_residual  (ipa)     &
                                                  + csite%dmean_energy_residual  (ipa)     &
                                                  * ndaysi
               csite%mmean_water_residual   (ipa) = csite%mmean_water_residual   (ipa)     &
                                                  + csite%dmean_water_residual   (ipa)     &
                                                  * ndaysi
               csite%mmean_rh               (ipa) = csite%mmean_rh               (ipa)     &
                                                  + csite%dmean_rh               (ipa)     &
                                                  * ndaysi
               csite%mmean_cwd_rh           (ipa) = csite%mmean_cwd_rh           (ipa)     &
                                                  + csite%dmean_cwd_rh           (ipa)     &
                                                  * ndaysi
               csite%mmean_nep              (ipa) = csite%mmean_nep              (ipa)     &
                                                  + csite%dmean_nep              (ipa)     &
                                                  * ndaysi
               csite%mmean_A_decomp         (ipa) = csite%mmean_A_decomp         (ipa)     &
                                                  + csite%dmean_A_decomp         (ipa)     &
                                                  * ndaysi
               csite%mmean_Af_decomp        (ipa) = csite%mmean_Af_decomp        (ipa)     &
                                                  + csite%dmean_Af_decomp        (ipa)     &
                                                  * ndaysi
               csite%mmean_rk4step          (ipa) = csite%mmean_rk4step          (ipa)     &
                                                  + csite%dmean_rk4step          (ipa)     &
                                                  * ndaysi
               csite%mmean_available_water  (ipa) = csite%mmean_available_water  (ipa)     &
                                                  + csite%dmean_available_water  (ipa)     &
                                                  * ndaysi
               csite%mmean_can_theiv        (ipa) = csite%mmean_can_theiv        (ipa)     &
                                                  + csite%dmean_can_theiv        (ipa)     &
                                                  * ndaysi
               csite%mmean_can_theta        (ipa) = csite%mmean_can_theta        (ipa)     &
                                                  + csite%dmean_can_theta        (ipa)     &
                                                  * ndaysi
               csite%mmean_can_vpdef        (ipa) = csite%mmean_can_vpdef        (ipa)     &
                                                  + csite%dmean_can_vpdef        (ipa)     &
                                                  * ndaysi
               csite%mmean_can_shv          (ipa) = csite%mmean_can_shv          (ipa)     &
                                                  + csite%dmean_can_shv          (ipa)     &
                                                  * ndaysi
               csite%mmean_can_co2          (ipa) = csite%mmean_can_co2          (ipa)     &
                                                  + csite%dmean_can_co2          (ipa)     &
                                                  * ndaysi
               csite%mmean_can_prss         (ipa) = csite%mmean_can_prss         (ipa)     &
                                                  + csite%dmean_can_prss         (ipa)     &
                                                  * ndaysi
               csite%mmean_gnd_temp         (ipa) = csite%mmean_gnd_temp         (ipa)     &
                                                  + csite%dmean_gnd_temp         (ipa)     &
                                                  * ndaysi
               csite%mmean_gnd_shv          (ipa) = csite%mmean_gnd_shv          (ipa)     &
                                                  + csite%dmean_gnd_shv          (ipa)     &
                                                  * ndaysi
               csite%mmean_can_ggnd         (ipa) = csite%mmean_can_ggnd         (ipa)     &
                                                  + csite%dmean_can_ggnd         (ipa)     &
                                                  * ndaysi
               csite%mmean_sfcw_depth       (ipa) = csite%mmean_sfcw_depth       (ipa)     &
                                                  + csite%dmean_sfcw_depth       (ipa)     &
                                                  * ndaysi
               !----- Temporarily make pounding energy extensive [J/m2]. ------------------!
               csite%mmean_sfcw_energy      (ipa) = csite%mmean_sfcw_energy      (ipa)     &
                                                  + csite%dmean_sfcw_energy      (ipa)     &
                                                  * csite%dmean_sfcw_mass        (ipa)     &
                                                  * ndaysi
               csite%mmean_sfcw_mass        (ipa) = csite%mmean_sfcw_mass        (ipa)     &
                                                  + csite%dmean_sfcw_mass        (ipa)     &
                                                  * ndaysi
               csite%mmean_soil_energy    (:,ipa) = csite%mmean_soil_energy    (:,ipa)     &
                                                  + csite%dmean_soil_energy    (:,ipa)     &
                                                  * ndaysi
               csite%mmean_soil_mstpot    (:,ipa) = csite%mmean_soil_mstpot    (:,ipa)     &
                                                  + csite%dmean_soil_mstpot    (:,ipa)     &
                                                  * ndaysi
               csite%mmean_soil_water     (:,ipa) = csite%mmean_soil_water     (:,ipa)     &
                                                  + csite%dmean_soil_water     (:,ipa)     &
                                                  * ndaysi
               csite%mmean_rshort_gnd       (ipa) = csite%mmean_rshort_gnd       (ipa)     &
                                                  + csite%dmean_rshort_gnd       (ipa)     &
                                                  * ndaysi
               csite%mmean_par_gnd          (ipa) = csite%mmean_par_gnd          (ipa)     &
                                                  + csite%dmean_par_gnd          (ipa)     &
                                                  * ndaysi
               csite%mmean_rlong_gnd        (ipa) = csite%mmean_rlong_gnd        (ipa)     &
                                                  + csite%dmean_rlong_gnd        (ipa)     &
                                                  * ndaysi
               csite%mmean_rlongup          (ipa) = csite%mmean_rlongup          (ipa)     &
                                                  + csite%dmean_rlongup          (ipa)     &
                                                  * ndaysi
               csite%mmean_parup            (ipa) = csite%mmean_parup            (ipa)     &
                                                  + csite%dmean_parup            (ipa)     &
                                                  * ndaysi
               csite%mmean_nirup            (ipa) = csite%mmean_nirup            (ipa)     &
                                                  + csite%dmean_nirup            (ipa)     &
                                                  * ndaysi
               csite%mmean_rshortup         (ipa) = csite%mmean_rshortup         (ipa)     &
                                                  + csite%dmean_rshortup         (ipa)     &
                                                  * ndaysi
               csite%mmean_rnet             (ipa) = csite%mmean_rnet             (ipa)     &
                                                  + csite%dmean_rnet             (ipa)     &
                                                  * ndaysi
               csite%mmean_albedo           (ipa) = csite%mmean_albedo           (ipa)     &
                                                  + csite%dmean_albedo           (ipa)     &
                                                  * ndaysi
               csite%mmean_albedo_par       (ipa) = csite%mmean_albedo_par       (ipa)     &
                                                  + csite%dmean_albedo_par       (ipa)     &
                                                  * ndaysi
               csite%mmean_albedo_nir       (ipa) = csite%mmean_albedo_nir       (ipa)     &
                                                  + csite%dmean_albedo_nir       (ipa)     &
                                                  * ndaysi
               csite%mmean_rlong_albedo     (ipa) = csite%mmean_rlong_albedo     (ipa)     &
                                                  + csite%dmean_rlong_albedo     (ipa)     &
                                                  * ndaysi
               csite%mmean_ustar            (ipa) = csite%mmean_ustar            (ipa)     &
                                                  + csite%dmean_ustar            (ipa)     &
                                                  * ndaysi
               csite%mmean_tstar            (ipa) = csite%mmean_tstar            (ipa)     &
                                                  + csite%dmean_tstar            (ipa)     &
                                                  * ndaysi
               csite%mmean_qstar            (ipa) = csite%mmean_qstar            (ipa)     &
                                                  + csite%dmean_qstar            (ipa)     &
                                                  * ndaysi
               csite%mmean_cstar            (ipa) = csite%mmean_cstar            (ipa)     &
                                                  + csite%dmean_cstar            (ipa)     &
                                                  * ndaysi
               csite%mmean_carbon_ac        (ipa) = csite%mmean_carbon_ac        (ipa)     &
                                                  + csite%dmean_carbon_ac        (ipa)     &
                                                  * ndaysi
               csite%mmean_carbon_st        (ipa) = csite%mmean_carbon_st        (ipa)     &
                                                  + csite%dmean_carbon_st        (ipa)     &
                                                  * ndaysi
               csite%mmean_vapor_gc         (ipa) = csite%mmean_vapor_gc         (ipa)     &
                                                  + csite%dmean_vapor_gc         (ipa)     &
                                                  * ndaysi
               csite%mmean_vapor_ac         (ipa) = csite%mmean_vapor_ac         (ipa)     &
                                                  + csite%dmean_vapor_ac         (ipa)     &
                                                  * ndaysi
               csite%mmean_smoist_gg      (:,ipa) = csite%mmean_smoist_gg      (:,ipa)     &
                                                  + csite%dmean_smoist_gg      (:,ipa)     &
                                                  * ndaysi
               csite%mmean_throughfall      (ipa) = csite%mmean_throughfall      (ipa)     &
                                                  + csite%dmean_throughfall      (ipa)     &
                                                  * ndaysi
               csite%mmean_transloss      (:,ipa) = csite%mmean_transloss      (:,ipa)     &
                                                  + csite%dmean_transloss      (:,ipa)     &
                                                  * ndaysi
               csite%mmean_runoff           (ipa) = csite%mmean_runoff           (ipa)     &
                                                  + csite%dmean_runoff           (ipa)     &
                                                  * ndaysi
               csite%mmean_drainage         (ipa) = csite%mmean_drainage         (ipa)     &
                                                  + csite%dmean_drainage         (ipa)     &
                                                  * ndaysi
               csite%mmean_sensible_gc      (ipa) = csite%mmean_sensible_gc      (ipa)     &
                                                  + csite%dmean_sensible_gc      (ipa)     &
                                                  * ndaysi
               csite%mmean_sensible_ac      (ipa) = csite%mmean_sensible_ac      (ipa)     &
                                                  + csite%dmean_sensible_ac      (ipa)     &
                                                  * ndaysi
               csite%mmean_sensible_gg    (:,ipa) = csite%mmean_sensible_gg    (:,ipa)     &
                                                  + csite%dmean_sensible_gg    (:,ipa)     &
                                                  * ndaysi
               csite%mmean_qthroughfall     (ipa) = csite%mmean_qthroughfall     (ipa)     &
                                                  + csite%dmean_qthroughfall     (ipa)     &
                                                  * ndaysi
               csite%mmean_qrunoff          (ipa) = csite%mmean_qrunoff          (ipa)     &
                                                  + csite%dmean_qrunoff          (ipa)     &
                                                  * ndaysi
               csite%mmean_qdrainage        (ipa) = csite%mmean_qdrainage        (ipa)     &
                                                  + csite%dmean_qdrainage        (ipa)     &
                                                  * ndaysi
               csite%mmean_A_decomp         (ipa) = csite%mmean_A_decomp         (ipa)     &
                                                  + csite%dmean_A_decomp         (ipa)     &
                                                  * ndaysi
               csite%mmean_Af_decomp        (ipa) = csite%mmean_Af_decomp        (ipa)     &
                                                  + csite%dmean_Af_decomp        (ipa)     &
                                                  * ndaysi
               csite%mmean_co2_residual     (ipa) = csite%mmean_co2_residual     (ipa)     &
                                                  + csite%dmean_co2_residual     (ipa)     &
                                                  * ndaysi
               csite%mmean_energy_residual  (ipa) = csite%mmean_energy_residual  (ipa)     &
                                                  + csite%dmean_energy_residual  (ipa)     &
                                                  * ndaysi
               csite%mmean_water_residual   (ipa) = csite%mmean_water_residual   (ipa)     &
                                                  + csite%dmean_water_residual   (ipa)     &
                                                  * ndaysi
               !---------------------------------------------------------------------------!
               !     Integrate the sum of squares.                                         !
               !---------------------------------------------------------------------------!
               csite%mmsqu_rh               (ipa) = csite%mmsqu_rh                  (ipa)  &
                                                  + isqu_ftz(csite%dmean_rh         (ipa)) &
                                                  * ndaysi
               csite%mmsqu_cwd_rh           (ipa) = csite%mmsqu_cwd_rh              (ipa)  &
                                                  + isqu_ftz(csite%dmean_cwd_rh     (ipa)) &
                                                  * ndaysi
               csite%mmsqu_nep              (ipa) = csite%mmsqu_nep                 (ipa)  &
                                                  + isqu_ftz(csite%dmean_nep        (ipa)) &
                                                  * ndaysi
               csite%mmsqu_rlongup          (ipa) = csite%mmsqu_rlongup             (ipa)  &
                                                  + isqu_ftz(csite%dmean_rlongup    (ipa)) &
                                                  * ndaysi
               csite%mmsqu_parup            (ipa) = csite%mmsqu_parup               (ipa)  &
                                                  + isqu_ftz(csite%dmean_parup      (ipa)) &
                                                  * ndaysi
               csite%mmsqu_nirup            (ipa) = csite%mmsqu_nirup               (ipa)  &
                                                  + isqu_ftz(csite%dmean_nirup      (ipa)) &
                                                  * ndaysi
               csite%mmsqu_rshortup         (ipa) = csite%mmsqu_rshortup            (ipa)  &
                                                  + isqu_ftz(csite%dmean_rshortup   (ipa)) &
                                                  * ndaysi
               csite%mmsqu_rnet             (ipa) = csite%mmsqu_rnet                (ipa)  &
                                                  + isqu_ftz(csite%dmean_rnet       (ipa)  &
                                                  * csite%dmean_rnet                (ipa)) &
                                                  * ndaysi
               csite%mmsqu_albedo           (ipa) = csite%mmsqu_albedo              (ipa)  &
                                                  + isqu_ftz(csite%dmean_albedo     (ipa)  &
                                                  * csite%dmean_albedo              (ipa)) &
                                                  * ndaysi
               csite%mmsqu_ustar            (ipa) = csite%mmsqu_ustar               (ipa)  &
                                                  + isqu_ftz(csite%dmean_ustar      (ipa)  &
                                                  * csite%dmean_ustar               (ipa)) &
                                                  * ndaysi
               csite%mmsqu_carbon_ac        (ipa) = csite%mmsqu_carbon_ac           (ipa)  &
                                                  + isqu_ftz(csite%dmean_carbon_ac  (ipa)) &
                                                  * ndaysi
               csite%mmsqu_carbon_st        (ipa) = csite%mmsqu_carbon_st           (ipa)  &
                                                  + isqu_ftz(csite%dmean_carbon_st  (ipa)) &
                                                  * ndaysi
               csite%mmsqu_vapor_gc         (ipa) = csite%mmsqu_vapor_gc            (ipa)  &
                                                  + isqu_ftz(csite%dmean_vapor_gc   (ipa)) &
                                                  * ndaysi
               csite%mmsqu_vapor_ac         (ipa) = csite%mmsqu_vapor_ac            (ipa)  &
                                                  + isqu_ftz(csite%dmean_vapor_ac   (ipa)) &
                                                  * ndaysi
               csite%mmsqu_sensible_gc      (ipa) = csite%mmsqu_sensible_gc         (ipa)  &
                                                  + isqu_ftz(csite%dmean_sensible_gc(ipa)) &
                                                  * ndaysi
               csite%mmsqu_sensible_ac      (ipa) = csite%mmsqu_sensible_ac         (ipa)  &
                                                  + isqu_ftz(csite%dmean_sensible_ac(ipa)) &
                                                  * ndaysi
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Patch loop.                                                          !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts
                  !------------------------------------------------------------------------!
                  !      Integrate the cohort-level variables that have no daily means     !
                  ! because their time step is one day.                                    !
                  !------------------------------------------------------------------------!
                  cpatch%mmean_lai             (ico) = cpatch%mmean_lai             (ico)  &
                                                     + cpatch%lai                   (ico)  &
                                                     * ndaysi
                  cpatch%mmean_bleaf           (ico) = cpatch%mmean_bleaf           (ico)  &
                                                     + cpatch%bleaf                 (ico)  &
                                                     * ndaysi
                  cpatch%mmean_broot           (ico) = cpatch%mmean_broot           (ico)  &
                                                     + cpatch%broot                 (ico)  &
                                                     * ndaysi
                  cpatch%mmean_bstorage        (ico) = cpatch%mmean_bstorage        (ico)  &
                                                     + cpatch%bstorage              (ico)  &
                                                     * ndaysi
                  cpatch%mmean_mort_rate     (:,ico) = cpatch%mmean_mort_rate     (:,ico)  &
                                                     + cpatch%mort_rate           (:,ico)  &
                                                     * ndaysi
                  cpatch%mmean_leaf_maintenance(ico) = cpatch%mmean_leaf_maintenance(ico)  &
                                                     + cpatch%leaf_maintenance      (ico)  &
                                                     * ndaysi
                  cpatch%mmean_root_maintenance(ico) = cpatch%mmean_root_maintenance(ico)  &
                                                     + cpatch%root_maintenance      (ico)  &
                                                     * ndaysi
                  cpatch%mmean_leaf_drop       (ico) = cpatch%mmean_leaf_drop       (ico)  &
                                                     + cpatch%leaf_drop             (ico)  &
                                                     * ndaysi
                  cpatch%mmean_cb              (ico) = cpatch%mmean_cb              (ico)  &
                                                     + ( cpatch%dmean_gpp           (ico)  &
                                                       - cpatch%dmean_plresp        (ico)  &
                                                       - cpatch%leaf_maintenance    (ico)  &
                                                       - cpatch%root_maintenance    (ico)  &
                                                       - cpatch%leaf_drop           (ico)  &
                                                       + cpatch%bstorage            (ico)  &
                                                       / yr_day )                          &
                                                     * ndaysi
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      Integrate the other cohort-level variables.                       !
                  !------------------------------------------------------------------------!
                  cpatch%mmean_gpp             (ico) = cpatch%mmean_gpp             (ico)  &
                                                     + cpatch%dmean_gpp             (ico)  &
                                                     * ndaysi
                  cpatch%mmean_npp             (ico) = cpatch%mmean_npp             (ico)  &
                                                     + cpatch%dmean_npp             (ico)  &
                                                     * ndaysi
                  cpatch%mmean_leaf_resp       (ico) = cpatch%mmean_leaf_resp       (ico)  &
                                                     + cpatch%dmean_leaf_resp       (ico)  &
                                                     * ndaysi
                  cpatch%mmean_root_resp       (ico) = cpatch%mmean_root_resp       (ico)  &
                                                     + cpatch%dmean_root_resp       (ico)  &
                                                     * ndaysi
                  cpatch%mmean_growth_resp     (ico) = cpatch%mmean_growth_resp     (ico)  &
                                                     + cpatch%dmean_growth_resp     (ico)  &
                                                     * ndaysi
                  cpatch%mmean_storage_resp    (ico) = cpatch%mmean_storage_resp    (ico)  &
                                                     + cpatch%dmean_storage_resp    (ico)  &
                                                     * ndaysi
                  cpatch%mmean_vleaf_resp      (ico) = cpatch%mmean_vleaf_resp      (ico)  &
                                                     + cpatch%dmean_vleaf_resp      (ico)  &
                                                     * ndaysi
                  cpatch%mmean_plresp          (ico) = cpatch%mmean_plresp          (ico)  &
                                                     + cpatch%dmean_plresp          (ico)  &
                                                     * ndaysi
                  cpatch%mmean_leaf_energy     (ico) = cpatch%mmean_leaf_energy     (ico)  &
                                                     + cpatch%dmean_leaf_energy     (ico)  &
                                                     * ndaysi
                  cpatch%mmean_leaf_water      (ico) = cpatch%mmean_leaf_water      (ico)  &
                                                     + cpatch%dmean_leaf_water      (ico)  &
                                                     * ndaysi
                  cpatch%mmean_leaf_hcap       (ico) = cpatch%mmean_leaf_hcap       (ico)  &
                                                     + cpatch%dmean_leaf_hcap       (ico)  &
                                                     * ndaysi
                  cpatch%mmean_leaf_vpdef      (ico) = cpatch%mmean_leaf_vpdef      (ico)  &
                                                     + cpatch%dmean_leaf_vpdef      (ico)  &
                                                     * ndaysi
                  cpatch%mmean_leaf_temp       (ico) = cpatch%mmean_leaf_temp       (ico)  &
                                                     + cpatch%dmean_leaf_temp       (ico)  &
                                                     * ndaysi
                  cpatch%mmean_leaf_fliq       (ico) = cpatch%mmean_leaf_fliq       (ico)  &
                                                     + cpatch%dmean_leaf_fliq       (ico)  &
                                                     * ndaysi
                  cpatch%mmean_leaf_gsw        (ico) = cpatch%mmean_leaf_gsw        (ico)  &
                                                     + cpatch%dmean_leaf_gsw        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_leaf_gbw        (ico) = cpatch%mmean_leaf_gbw        (ico)  &
                                                     + cpatch%dmean_leaf_gbw        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_wood_energy     (ico) = cpatch%mmean_wood_energy     (ico)  &
                                                     + cpatch%dmean_wood_energy     (ico)  &
                                                     * ndaysi
                  cpatch%mmean_wood_water      (ico) = cpatch%mmean_wood_water      (ico)  &
                                                     + cpatch%dmean_wood_water      (ico)  &
                                                     * ndaysi
                  cpatch%mmean_wood_hcap       (ico) = cpatch%mmean_wood_hcap       (ico)  &
                                                     + cpatch%dmean_wood_hcap       (ico)  &
                                                     * ndaysi
                  cpatch%mmean_wood_temp       (ico) = cpatch%mmean_wood_temp       (ico)  &
                                                     + cpatch%dmean_wood_temp       (ico)  &
                                                     * ndaysi
                  cpatch%mmean_wood_fliq       (ico) = cpatch%mmean_wood_fliq       (ico)  &
                                                     + cpatch%dmean_wood_fliq       (ico)  &
                                                     * ndaysi
                  cpatch%mmean_wood_gbw        (ico) = cpatch%mmean_wood_gbw        (ico)  &
                                                     + cpatch%dmean_wood_gbw        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_fs_open         (ico) = cpatch%mmean_fs_open         (ico)  &
                                                     + cpatch%dmean_fs_open         (ico)  &
                                                     * ndaysi
                  cpatch%mmean_fsw             (ico) = cpatch%mmean_fsw             (ico)  &
                                                     + cpatch%dmean_fsw             (ico)  &
                                                     * ndaysi
                  cpatch%mmean_fsn             (ico) = cpatch%mmean_fsn             (ico)  &
                                                     + cpatch%dmean_fsn             (ico)  &
                                                     * ndaysi
                  cpatch%mmean_a_light         (ico) = cpatch%mmean_a_light         (ico)  &
                                                     + cpatch%dmean_a_light         (ico)  &
                                                     * ndaysi
                  cpatch%mmean_a_rubp          (ico) = cpatch%mmean_a_rubp          (ico)  &
                                                     + cpatch%dmean_a_rubp          (ico)  &
                                                     * ndaysi
                  cpatch%mmean_a_co2           (ico) = cpatch%mmean_a_co2           (ico)  &
                                                     + cpatch%dmean_a_co2           (ico)  &
                                                     * ndaysi
                  cpatch%mmean_psi_open        (ico) = cpatch%mmean_psi_open        (ico)  &
                                                     + cpatch%dmean_psi_open        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_psi_closed      (ico) = cpatch%mmean_psi_closed      (ico)  &
                                                     + cpatch%dmean_psi_closed      (ico)  &
                                                     * ndaysi
                  cpatch%mmean_water_supply    (ico) = cpatch%mmean_water_supply    (ico)  &
                                                     + cpatch%dmean_water_supply    (ico)  &
                                                     * ndaysi

                  cpatch%mmean_par_level_beam  (ico) = cpatch%mmean_par_level_beam  (ico)  &
                                                     + cpatch%dmean_par_level_beam  (ico)  &
                                                     * ndaysi

                  cpatch%mmean_par_level_diffd (ico) = cpatch%mmean_par_level_diffd  (ico)  &
                                                     + cpatch%dmean_par_level_diffd  (ico)  &
                                                     * ndaysi

                  cpatch%mmean_par_level_diffu (ico) = cpatch%mmean_par_level_diffu  (ico)  &
                                                     + cpatch%dmean_par_level_diffu  (ico)  &
                                                     * ndaysi

                  cpatch%mmean_light_level     (ico) = cpatch%mmean_light_level     (ico)  &
                                                     + cpatch%dmean_light_level     (ico)  &
                                                     * ndaysi
                  cpatch%mmean_light_level_beam(ico) = cpatch%mmean_light_level_beam(ico)  &
                                                     + cpatch%dmean_light_level_beam(ico)  &
                                                     * ndaysi
                  cpatch%mmean_light_level_diff(ico) = cpatch%mmean_light_level_diff(ico)  &
                                                     + cpatch%dmean_light_level_diff(ico)  &
                                                     * ndaysi


                  cpatch%mmean_par_l           (ico) = cpatch%mmean_par_l           (ico)  &
                                                     + cpatch%dmean_par_l           (ico)  &
                                                     * ndaysi
                  cpatch%mmean_par_l_beam      (ico) = cpatch%mmean_par_l_beam      (ico)  &
                                                     + cpatch%dmean_par_l_beam      (ico)  &
                                                     * ndaysi
                  cpatch%mmean_par_l_diff      (ico) = cpatch%mmean_par_l_diff      (ico)  &
                                                     + cpatch%dmean_par_l_diff      (ico)  &
                                                     * ndaysi
                  cpatch%mmean_rshort_l        (ico) = cpatch%mmean_rshort_l        (ico)  &
                                                     + cpatch%dmean_rshort_l        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_rlong_l         (ico) = cpatch%mmean_rlong_l         (ico)  &
                                                     + cpatch%dmean_rlong_l         (ico)  &
                                                     * ndaysi
                  cpatch%mmean_sensible_lc     (ico) = cpatch%mmean_sensible_lc     (ico)  &
                                                     + cpatch%dmean_sensible_lc     (ico)  &
                                                     * ndaysi
                  cpatch%mmean_vapor_lc        (ico) = cpatch%mmean_vapor_lc        (ico)  &
                                                     + cpatch%dmean_vapor_lc        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_transp          (ico) = cpatch%mmean_transp          (ico)  &
                                                     + cpatch%dmean_transp          (ico)  &
                                                     * ndaysi
                  cpatch%mmean_intercepted_al  (ico) = cpatch%mmean_intercepted_al  (ico)  &
                                                     + cpatch%dmean_intercepted_al  (ico)  &
                                                     * ndaysi
                  cpatch%mmean_wshed_lg        (ico) = cpatch%mmean_wshed_lg        (ico)  &
                                                     + cpatch%dmean_wshed_lg        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_rshort_w        (ico) = cpatch%mmean_rshort_w        (ico)  &
                                                     + cpatch%dmean_rshort_w        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_rlong_w         (ico) = cpatch%mmean_rlong_w         (ico)  &
                                                     + cpatch%dmean_rlong_w         (ico)  &
                                                     * ndaysi
                  cpatch%mmean_rad_profile   (:,ico) = cpatch%mmean_rad_profile   (:,ico)  &
                                                     + cpatch%dmean_rad_profile   (:,ico)  &
                                                     * ndaysi
                  cpatch%mmean_sensible_wc     (ico) = cpatch%mmean_sensible_wc     (ico)  &
                                                     + cpatch%dmean_sensible_wc     (ico)  &
                                                     * ndaysi
                  cpatch%mmean_vapor_wc        (ico) = cpatch%mmean_vapor_wc        (ico)  &
                                                     + cpatch%dmean_vapor_wc        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_intercepted_aw  (ico) = cpatch%mmean_intercepted_aw  (ico)  &
                                                     + cpatch%dmean_intercepted_aw  (ico)  &
                                                     * ndaysi
                  cpatch%mmean_wshed_wg        (ico) = cpatch%mmean_wshed_wg        (ico)  &
                                                     + cpatch%dmean_wshed_wg        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_nppleaf         (ico) = cpatch%mmean_nppleaf         (ico)  &
                                                     + cpatch%dmean_nppleaf         (ico)  &
                                                     * ndaysi
                  cpatch%mmean_nppfroot        (ico) = cpatch%mmean_nppfroot        (ico)  &
                                                     + cpatch%dmean_nppfroot        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_nppsapwood      (ico) = cpatch%mmean_nppsapwood      (ico)  &
                                                     + cpatch%dmean_nppsapwood      (ico)  &
                                                     * ndaysi
                  cpatch%mmean_nppcroot        (ico) = cpatch%mmean_nppcroot        (ico)  &
                                                     + cpatch%dmean_nppcroot        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_nppseeds        (ico) = cpatch%mmean_nppseeds        (ico)  &
                                                     + cpatch%dmean_nppseeds        (ico)  &
                                                     * ndaysi
                  cpatch%mmean_nppwood         (ico) = cpatch%mmean_nppwood         (ico)  &
                                                     + cpatch%dmean_nppwood         (ico)  &
                                                     * ndaysi
                  cpatch%mmean_nppdaily        (ico) = cpatch%mmean_nppdaily        (ico)  &
                                                     + cpatch%dmean_nppdaily        (ico)  &
                                                     * ndaysi
                  !----- Integrate mean sum of squares. -----------------------------------!
                  cpatch%mmsqu_gpp        (ico) = cpatch%mmsqu_gpp                  (ico)  &
                                                + isqu_ftz(cpatch%dmean_gpp         (ico)) &
                                                * ndaysi
                  cpatch%mmsqu_npp        (ico) = cpatch%mmsqu_npp                  (ico)  &
                                                + isqu_ftz(cpatch%dmean_npp         (ico)) &
                                                * ndaysi
                  cpatch%mmsqu_plresp     (ico) = cpatch%mmsqu_plresp               (ico)  &
                                                + isqu_ftz(cpatch%dmean_plresp      (ico)) &
                                                * ndaysi
                  cpatch%mmsqu_sensible_lc(ico) = cpatch%mmsqu_sensible_lc          (ico)  &
                                                + isqu_ftz(cpatch%dmean_sensible_lc (ico)) &
                                                * ndaysi
                  cpatch%mmsqu_vapor_lc   (ico) = cpatch%mmsqu_vapor_lc             (ico)  &
                                                + isqu_ftz(cpatch%dmean_vapor_lc    (ico)) &
                                                * ndaysi
                  cpatch%mmsqu_transp     (ico) = cpatch%mmsqu_transp               (ico)  &
                                                + isqu_ftz(cpatch%dmean_transp      (ico)) &
                                                * ndaysi
                  cpatch%mmsqu_sensible_wc(ico) = cpatch%mmsqu_sensible_wc          (ico)  &
                                                + isqu_ftz(cpatch%dmean_sensible_wc (ico)) &
                                                * ndaysi
                  cpatch%mmsqu_vapor_wc   (ico) = cpatch%mmsqu_vapor_wc             (ico)  &
                                                + isqu_ftz(cpatch%dmean_vapor_wc    (ico)) &
                                                * ndaysi
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine integrate_ed_mmean_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine normalises the daily mean variables of those variables that could  !
   ! not be integrated directly.  This includes temperatures, liquid fraction, and soil    !
   ! matric potential.                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine normalize_ed_mmean_vars(cgrid)
      use ed_state_vars        , only : edtype             & ! structure
                                      , polygontype        & ! structure
                                      , sitetype           & ! structure
                                      , patchtype          ! ! structure
      use grid_coms            , only : nzg                ! ! intent(in)
      use ed_misc_coms         , only : dtlsm              ! ! intent(in)
      use therm_lib            , only : press2exner        & ! function
                                      , extheta2temp       & ! function
                                      , uextcm2tl          & ! subroutine
                                      , uint2tl            & ! subroutine
                                      , idealdenssh        ! ! function
      use soil_coms            , only : tiny_sfcwater_mass & ! intent(in)
                                      , soil               ! ! intent(in)
      use consts_coms          , only : t00                & ! intent(in)
                                      , wdns               ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)                          , target  :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype)                     , pointer :: cpoly
      type(sitetype)                        , pointer :: csite
      type(patchtype)                       , pointer :: cpatch
      real             , dimension(nzg)               :: cgrid_mmean_soil_hcap
      integer                                         :: ipy
      integer                                         :: isi
      integer                                         :: ipa
      integer                                         :: ico
      integer                                         :: k
      integer                                         :: nsoil
      real                                            :: can_exner
      real                                            :: atm_exner
      real                                            :: site_area_i
      real                                            :: poly_area_i
      real                                            :: site_wgt
      real                                            :: patch_wgt
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Loop over polygons.                                                            !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)


         !----- Inverse of this polygon area (it should be always 1.) ---------------------!
         poly_area_i = 1./sum(cpoly%area)
         !---------------------------------------------------------------------------------!


         !----- Re-set some support variables. --------------------------------------------!
         cgrid_mmean_soil_hcap(:) = 0.0
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Loop over sites.                                                            !
         !---------------------------------------------------------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            !----- Inverse of this site area (it should be always 1.) ---------------------!
            site_area_i = 1./sum(csite%area)
            !------------------------------------------------------------------------------!


            !----- Site weight. -----------------------------------------------------------!
            site_wgt = cpoly%area(isi) * poly_area_i
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Find the derived properties for the air above canopy.                   !
            !------------------------------------------------------------------------------!
            atm_exner                 = press2exner (cpoly%mmean_atm_prss(isi))
            cpoly%mmean_atm_temp(isi) = extheta2temp(atm_exner,cpoly%mmean_atm_theta(isi))
            cpoly%mmean_atm_rhos(isi) = idealdenssh ( cpoly%mmean_atm_prss  (isi)          &
                                                    , cpoly%mmean_atm_temp  (isi)          &
                                                    , cpoly%mmean_atm_shv   (isi) )
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Loop over patches.                                                       !
            !------------------------------------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)


               !----- Site weight. --------------------------------------------------------!
               patch_wgt = csite%area(ipa) * site_area_i * site_wgt
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !      Now we find the derived properties for the canopy air space.         !
               !---------------------------------------------------------------------------!
               can_exner                 = press2exner (csite%mmean_can_prss(ipa))
               csite%mmean_can_temp(ipa) = extheta2temp( can_exner                         &
                                                       , csite%mmean_can_theta(ipa))
               csite%mmean_can_rhos(ipa) = idealdenssh ( csite%mmean_can_prss  (ipa)       &
                                                       , csite%mmean_can_temp  (ipa)       &
                                                       , csite%mmean_can_shv   (ipa)       )
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !     Soil matric potential, temperature, and liquid water.                 !
               !---------------------------------------------------------------------------!
               do k=1,nzg
                  nsoil = cpoly%ntext_soil(k,isi)
                  call uextcm2tl( csite%mmean_soil_energy(k,ipa)                           &
                                , csite%mmean_soil_water (k,ipa) * wdns                    &
                                , soil(nsoil)%slcpd                                        &
                                , csite%mmean_soil_temp  (k,ipa)                           &
                                , csite%mmean_soil_fliq  (k,ipa))

                  cgrid_mmean_soil_hcap   (k)     = cgrid_mmean_soil_hcap(k)               &
                                                  + soil(nsoil)%slcpd * patch_wgt
                  !------------------------------------------------------------------------!

               end do
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !   If the patch had some temporary snow/pounding layer, convert the mean   !
               ! energy to J/kg, then find the mean temperature and liquid fraction.       !
               ! Otherwise, set them to either zero or default values.                     !
               !---------------------------------------------------------------------------!
               if (csite%mmean_sfcw_mass(ipa) > tiny_sfcwater_mass) then
                  csite%mmean_sfcw_energy(ipa) = csite%mmean_sfcw_energy(ipa)              &
                                               / csite%mmean_sfcw_mass  (ipa)
                  call uint2tl( csite%mmean_sfcw_energy(ipa), csite%mmean_sfcw_temp(ipa)   &
                              , csite%mmean_sfcw_fliq  (ipa))
               else
                  csite%mmean_sfcw_mass  (ipa)  = 0.
                  csite%mmean_sfcw_depth (ipa)  = 0.
                  csite%mmean_sfcw_energy(ipa)  = 0.
                  csite%mmean_sfcw_temp  (ipa)  = csite%mmean_soil_temp(nzg,ipa)
                  csite%mmean_sfcw_fliq  (ipa)  = csite%mmean_soil_fliq(nzg,ipa)
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Loop over the cohorts.                                               !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts



                  !------------------------------------------------------------------------!
                  !     Find the vegetation temperature and liquid fraction.               !
                  !------------------------------------------------------------------------!
                  !----- Leaf. ------------------------------------------------------------!
                  if (cpatch%mmean_leaf_hcap(ico) > 0.) then
                     call uextcm2tl( cpatch%mmean_leaf_energy(ico)                         &
                                   , cpatch%mmean_leaf_water (ico)                         &
                                   , cpatch%mmean_leaf_hcap  (ico)                         &
                                   , cpatch%mmean_leaf_temp  (ico)                         &
                                   , cpatch%mmean_leaf_fliq  (ico) )
                  else
                     cpatch%mmean_leaf_vpdef(ico) = csite%mmean_can_vpdef(ipa)
                     cpatch%mmean_leaf_temp (ico) = csite%mmean_can_temp (ipa)
                     if (csite%mmean_can_temp(ipa) > t00) then
                        cpatch%mmean_leaf_fliq(ico) = 1.0
                     elseif (csite%mmean_can_temp(ipa) == t00) then
                        cpatch%mmean_leaf_fliq(ico) = 0.5
                     else
                        cpatch%mmean_leaf_fliq(ico) = 0.0
                     end if
                  end if
                  !----- Wood. ------------------------------------------------------------!
                  if (cpatch%mmean_wood_hcap(ico) > 0.) then
                     call uextcm2tl( cpatch%mmean_wood_energy(ico)                         &
                                   , cpatch%mmean_wood_water (ico)                         &
                                   , cpatch%mmean_wood_hcap  (ico)                         &
                                   , cpatch%mmean_wood_temp  (ico)                         &
                                   , cpatch%mmean_wood_fliq  (ico) )
                  else
                     cpatch%mmean_wood_temp(ico) = csite%mmean_can_temp(ipa)
                     if (csite%mmean_can_temp(ipa) > t00) then
                        cpatch%mmean_wood_fliq(ico) = 1.0
                     elseif (csite%mmean_can_temp(ipa) == t00) then
                        cpatch%mmean_wood_fliq(ico) = 0.5
                     else
                        cpatch%mmean_wood_fliq(ico) = 0.0
                     end if
                  end if
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!






         !---------------------------------------------------------------------------------!
         !      Find the derived properties for the air above canopy.                      !
         !---------------------------------------------------------------------------------!
         atm_exner                 = press2exner (cgrid%mmean_atm_prss(ipy))
         cgrid%mmean_atm_temp(ipy) = extheta2temp(atm_exner,cgrid%mmean_atm_theta(ipy))
         cgrid%mmean_atm_rhos(ipy) = idealdenssh ( cgrid%mmean_atm_prss  (ipy)             &
                                                 , cgrid%mmean_atm_temp  (ipy)             &
                                                 , cgrid%mmean_atm_shv   (ipy) )
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !      Find the derived properties for the canopy air space.                      !
         !---------------------------------------------------------------------------------!
         can_exner                 = press2exner (cgrid%mmean_can_prss(ipy))
         cgrid%mmean_can_temp(ipy) = extheta2temp(can_exner,cgrid%mmean_can_theta(ipy))
         cgrid%mmean_can_rhos(ipy) = idealdenssh ( cgrid%mmean_can_prss  (ipy)             &
                                                 , cgrid%mmean_can_temp  (ipy)             &
                                                 , cgrid%mmean_can_shv   (ipy) )
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !   If the patch had some temporary snow/pounding layer, convert the mean energy  !
         ! to J/kg, then find the mean temperature and liquid fraction.  Otherwise, set    !
         ! them to either zero or default values.                                          !
         !---------------------------------------------------------------------------------!
         if (cgrid%mmean_sfcw_mass(ipy) > tiny_sfcwater_mass) then
            cgrid%mmean_sfcw_energy(ipy) = cgrid%mmean_sfcw_energy(ipy)                    &
                                         / cgrid%mmean_sfcw_mass(ipy)
            call uint2tl(cgrid%mmean_sfcw_energy(ipy),cgrid%mmean_sfcw_temp(ipy)           &
                        ,cgrid%mmean_sfcw_fliq(ipy))
         else
            cgrid%mmean_sfcw_mass  (ipy)  = 0.
            cgrid%mmean_sfcw_depth (ipy)  = 0.
            cgrid%mmean_sfcw_energy(ipy)  = 0.
            cgrid%mmean_sfcw_temp  (ipy)  = cgrid%mmean_soil_temp(nzg,ipy)
            cgrid%mmean_sfcw_fliq  (ipy)  = cgrid%mmean_soil_fliq(nzg,ipy)
         end if
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Find the temperature and the fraction of liquid water.                      !
         !---------------------------------------------------------------------------------!
         do k=1,nzg
            call uextcm2tl( cgrid%mmean_soil_energy(k,ipy)                                 &
                          , cgrid%mmean_soil_water (k,ipy) * wdns                          &
                          , cgrid_mmean_soil_hcap  (k)                                     &
                          , cgrid%mmean_soil_temp  (k,ipy)                                 &
                          , cgrid%mmean_soil_fliq  (k,ipy) )
         end do
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !     Find the vegetation temperature and liquid fraction.                        !
         !---------------------------------------------------------------------------------!
         !----- Leaf. ---------------------------------------------------------------------!
         if (cgrid%mmean_leaf_hcap(ipy) > 0.) then
            call uextcm2tl( cgrid%mmean_leaf_energy(ipy), cgrid%mmean_leaf_water (ipy)     &
                          , cgrid%mmean_leaf_hcap  (ipy), cgrid%mmean_leaf_temp  (ipy)     &
                          , cgrid%mmean_leaf_fliq  (ipy) )
         else
            cgrid%mmean_leaf_temp (ipy) = cgrid%mmean_can_temp (ipy)
            if (cgrid%mmean_can_temp(ipy) > t00) then
               cgrid%mmean_leaf_fliq(ipy) = 1.0
            elseif (cgrid%mmean_can_temp(ipy) == t00) then
               cgrid%mmean_leaf_fliq(ipy) = 0.5
            else
               cgrid%mmean_leaf_fliq(ipy) = 0.0
            end if
         end if
         !----- Wood. ---------------------------------------------------------------------!
         if (cgrid%mmean_wood_hcap(ipy) > 0.) then
            call uextcm2tl( cgrid%mmean_wood_energy(ipy)                                   &
                          , cgrid%mmean_wood_water (ipy)                                   &
                          , cgrid%mmean_wood_hcap  (ipy)                                   &
                          , cgrid%mmean_wood_temp  (ipy)                                   &
                          , cgrid%mmean_wood_fliq  (ipy) )
         else
            cgrid%mmean_wood_temp(ipy) = cgrid%mmean_can_temp(ipy)
            if (cgrid%mmean_can_temp(ipy) > t00) then
               cgrid%mmean_wood_fliq(ipy) = 1.0
            elseif (cgrid%mmean_can_temp(ipy) == t00) then
               cgrid%mmean_wood_fliq(ipy) = 0.5
            else
               cgrid%mmean_wood_fliq(ipy) = 0.0
            end if
         end if
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine normalize_ed_mmean_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine resets the monthly averages for variables actually used in the     !
   ! integration.                                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine zero_ed_mmean_vars(cgrid)
      use ed_state_vars , only : edtype        & ! structure
                               , polygontype   & ! structure
                               , sitetype      & ! structure
                               , patchtype     ! ! structure
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target  :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer :: cpoly
      type(sitetype)   , pointer :: csite
      type(patchtype)  , pointer :: cpatch
      integer                    :: ipy
      integer                    :: isi
      integer                    :: ipa
      integer                    :: ico
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !       Loop over polygons.                                                          !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         cgrid%mmean_lai             (:,:,ipy) = 0.0
         cgrid%mmean_bleaf           (:,:,ipy) = 0.0
         cgrid%mmean_broot           (:,:,ipy) = 0.0
         cgrid%mmean_bstorage        (:,:,ipy) = 0.0
         cgrid%mmean_bleaf_n         (:,:,ipy) = 0.0
         cgrid%mmean_broot_n         (:,:,ipy) = 0.0
         cgrid%mmean_bstorage_n      (:,:,ipy) = 0.0
         cgrid%mmean_leaf_maintenance(:,:,ipy) = 0.0
         cgrid%mmean_root_maintenance(:,:,ipy) = 0.0
         cgrid%mmean_leaf_drop       (:,:,ipy) = 0.0
         cgrid%mmean_fast_soil_c         (ipy) = 0.0 
         cgrid%mmean_slow_soil_c         (ipy) = 0.0 
         cgrid%mmean_struct_soil_c       (ipy) = 0.0 
         cgrid%mmean_struct_soil_l       (ipy) = 0.0 
         cgrid%mmean_cwd_c               (ipy) = 0.0 
         cgrid%mmean_fast_soil_n         (ipy) = 0.0 
         cgrid%mmean_mineral_soil_n      (ipy) = 0.0 
         cgrid%mmean_cwd_n               (ipy) = 0.0 
         cgrid%mmean_gpp                 (ipy) = 0.0 
         cgrid%mmean_npp                 (ipy) = 0.0 
         cgrid%mmean_leaf_resp           (ipy) = 0.0 
         cgrid%mmean_root_resp           (ipy) = 0.0 
         cgrid%mmean_growth_resp         (ipy) = 0.0 
         cgrid%mmean_storage_resp        (ipy) = 0.0 
         cgrid%mmean_vleaf_resp          (ipy) = 0.0 
         cgrid%mmean_plresp              (ipy) = 0.0 
         cgrid%mmean_leaf_energy         (ipy) = 0.0 
         cgrid%mmean_leaf_water          (ipy) = 0.0 
         cgrid%mmean_leaf_hcap           (ipy) = 0.0 
         cgrid%mmean_leaf_vpdef          (ipy) = 0.0 
         cgrid%mmean_leaf_temp           (ipy) = 0.0 
         cgrid%mmean_leaf_fliq           (ipy) = 0.0 
         cgrid%mmean_leaf_gsw            (ipy) = 0.0 
         cgrid%mmean_leaf_gbw            (ipy) = 0.0 
         cgrid%mmean_wood_energy         (ipy) = 0.0 
         cgrid%mmean_wood_water          (ipy) = 0.0 
         cgrid%mmean_wood_hcap           (ipy) = 0.0 
         cgrid%mmean_wood_temp           (ipy) = 0.0 
         cgrid%mmean_wood_fliq           (ipy) = 0.0 
         cgrid%mmean_wood_gbw            (ipy) = 0.0 
         cgrid%mmean_fs_open             (ipy) = 0.0 
         cgrid%mmean_fsw                 (ipy) = 0.0 
         cgrid%mmean_fsn                 (ipy) = 0.0 
         cgrid%mmean_a_light             (ipy) = 0.0
         cgrid%mmean_a_rubp              (ipy) = 0.0
         cgrid%mmean_a_co2               (ipy) = 0.0
         cgrid%mmean_psi_open            (ipy) = 0.0 
         cgrid%mmean_psi_closed          (ipy) = 0.0 
         cgrid%mmean_water_supply        (ipy) = 0.0 
         cgrid%mmean_par_l               (ipy) = 0.0 
         cgrid%mmean_par_l_beam          (ipy) = 0.0 
         cgrid%mmean_par_l_diff          (ipy) = 0.0 
         cgrid%mmean_rshort_l            (ipy) = 0.0 
         cgrid%mmean_rlong_l             (ipy) = 0.0 
         cgrid%mmean_sensible_lc         (ipy) = 0.0 
         cgrid%mmean_vapor_lc            (ipy) = 0.0 
         cgrid%mmean_transp              (ipy) = 0.0 
         cgrid%mmean_intercepted_al      (ipy) = 0.0 
         cgrid%mmean_wshed_lg            (ipy) = 0.0 
         cgrid%mmean_rshort_w            (ipy) = 0.0 
         cgrid%mmean_rlong_w             (ipy) = 0.0 
         cgrid%mmean_sensible_wc         (ipy) = 0.0 
         cgrid%mmean_vapor_wc            (ipy) = 0.0 
         cgrid%mmean_intercepted_aw      (ipy) = 0.0 
         cgrid%mmean_wshed_wg            (ipy) = 0.0 
         cgrid%mmean_nppleaf             (ipy) = 0.0 
         cgrid%mmean_nppfroot            (ipy) = 0.0 
         cgrid%mmean_nppsapwood          (ipy) = 0.0 
         cgrid%mmean_nppcroot            (ipy) = 0.0 
         cgrid%mmean_nppseeds            (ipy) = 0.0 
         cgrid%mmean_nppwood             (ipy) = 0.0 
         cgrid%mmean_nppdaily            (ipy) = 0.0 
         cgrid%mmean_rh                  (ipy) = 0.0 
         cgrid%mmean_cwd_rh              (ipy) = 0.0 
         cgrid%mmean_nep                 (ipy) = 0.0 
         cgrid%mmean_rk4step             (ipy) = 0.0 
         cgrid%mmean_available_water     (ipy) = 0.0 
         cgrid%mmean_can_theiv           (ipy) = 0.0 
         cgrid%mmean_can_theta           (ipy) = 0.0 
         cgrid%mmean_can_vpdef           (ipy) = 0.0 
         cgrid%mmean_can_temp            (ipy) = 0.0 
         cgrid%mmean_can_shv             (ipy) = 0.0 
         cgrid%mmean_can_co2             (ipy) = 0.0 
         cgrid%mmean_can_rhos            (ipy) = 0.0 
         cgrid%mmean_can_prss            (ipy) = 0.0 
         cgrid%mmean_gnd_temp            (ipy) = 0.0 
         cgrid%mmean_gnd_shv             (ipy) = 0.0 
         cgrid%mmean_can_ggnd            (ipy) = 0.0 
         cgrid%mmean_sfcw_depth          (ipy) = 0.0 
         cgrid%mmean_sfcw_energy         (ipy) = 0.0 
         cgrid%mmean_sfcw_mass           (ipy) = 0.0 
         cgrid%mmean_sfcw_temp           (ipy) = 0.0 
         cgrid%mmean_sfcw_fliq           (ipy) = 0.0 
         cgrid%mmean_soil_energy       (:,ipy) = 0.0 
         cgrid%mmean_soil_mstpot       (:,ipy) = 0.0 
         cgrid%mmean_soil_water        (:,ipy) = 0.0 
         cgrid%mmean_soil_temp         (:,ipy) = 0.0 
         cgrid%mmean_soil_fliq         (:,ipy) = 0.0 
         cgrid%mmean_rshort_gnd          (ipy) = 0.0 
         cgrid%mmean_par_gnd             (ipy) = 0.0 
         cgrid%mmean_rlong_gnd           (ipy) = 0.0 
         cgrid%mmean_rlongup             (ipy) = 0.0 
         cgrid%mmean_parup               (ipy) = 0.0 
         cgrid%mmean_nirup               (ipy) = 0.0 
         cgrid%mmean_rshortup            (ipy) = 0.0 
         cgrid%mmean_rnet                (ipy) = 0.0 
         cgrid%mmean_albedo              (ipy) = 0.0 
         cgrid%mmean_albedo_par          (ipy) = 0.0 
         cgrid%mmean_albedo_nir          (ipy) = 0.0 
         cgrid%mmean_rlong_albedo        (ipy) = 0.0 
         cgrid%mmean_ustar               (ipy) = 0.0 
         cgrid%mmean_tstar               (ipy) = 0.0 
         cgrid%mmean_qstar               (ipy) = 0.0 
         cgrid%mmean_cstar               (ipy) = 0.0 
         cgrid%mmean_carbon_ac           (ipy) = 0.0 
         cgrid%mmean_carbon_st           (ipy) = 0.0 
         cgrid%mmean_vapor_gc            (ipy) = 0.0 
         cgrid%mmean_vapor_ac            (ipy) = 0.0 
         cgrid%mmean_smoist_gg         (:,ipy) = 0.0 
         cgrid%mmean_throughfall         (ipy) = 0.0 
         cgrid%mmean_transloss         (:,ipy) = 0.0 
         cgrid%mmean_runoff              (ipy) = 0.0 
         cgrid%mmean_drainage            (ipy) = 0.0 
         cgrid%mmean_sensible_gc         (ipy) = 0.0 
         cgrid%mmean_sensible_ac         (ipy) = 0.0 
         cgrid%mmean_sensible_gg       (:,ipy) = 0.0 
         cgrid%mmean_qthroughfall        (ipy) = 0.0 
         cgrid%mmean_qrunoff             (ipy) = 0.0 
         cgrid%mmean_qdrainage           (ipy) = 0.0 
         cgrid%mmean_nppleaf             (ipy) = 0.0
         cgrid%mmean_nppfroot            (ipy) = 0.0
         cgrid%mmean_nppsapwood          (ipy) = 0.0
         cgrid%mmean_nppcroot            (ipy) = 0.0
         cgrid%mmean_nppseeds            (ipy) = 0.0
         cgrid%mmean_nppwood             (ipy) = 0.0
         cgrid%mmean_nppdaily            (ipy) = 0.0
         cgrid%mmean_A_decomp            (ipy) = 0.0 
         cgrid%mmean_Af_decomp           (ipy) = 0.0 
         cgrid%mmean_co2_residual        (ipy) = 0.0 
         cgrid%mmean_energy_residual     (ipy) = 0.0 
         cgrid%mmean_water_residual      (ipy) = 0.0 
         cgrid%mmean_atm_theiv           (ipy) = 0.0 
         cgrid%mmean_atm_theta           (ipy) = 0.0 
         cgrid%mmean_atm_temp            (ipy) = 0.0 
         cgrid%mmean_atm_vpdef           (ipy) = 0.0 
         cgrid%mmean_atm_shv             (ipy) = 0.0 
         cgrid%mmean_atm_rshort          (ipy) = 0.0 
         cgrid%mmean_atm_rshort_diff     (ipy) = 0.0 
         cgrid%mmean_atm_par             (ipy) = 0.0 
         cgrid%mmean_atm_par_diff        (ipy) = 0.0 
         cgrid%mmean_atm_rlong           (ipy) = 0.0 
         cgrid%mmean_atm_vels            (ipy) = 0.0 
         cgrid%mmean_atm_rhos            (ipy) = 0.0 
         cgrid%mmean_atm_prss            (ipy) = 0.0 
         cgrid%mmean_atm_co2             (ipy) = 0.0 
         cgrid%mmean_pcpg                (ipy) = 0.0 
         cgrid%mmean_qpcpg               (ipy) = 0.0 
         cgrid%mmean_dpcpg               (ipy) = 0.0 
         cgrid%mmsqu_gpp                 (ipy) = 0.0 
         cgrid%mmsqu_npp                 (ipy) = 0.0 
         cgrid%mmsqu_plresp              (ipy) = 0.0 
         cgrid%mmsqu_sensible_lc         (ipy) = 0.0 
         cgrid%mmsqu_vapor_lc            (ipy) = 0.0 
         cgrid%mmsqu_transp              (ipy) = 0.0 
         cgrid%mmsqu_sensible_wc         (ipy) = 0.0 
         cgrid%mmsqu_vapor_wc            (ipy) = 0.0 
         cgrid%mmsqu_rh                  (ipy) = 0.0 
         cgrid%mmsqu_cwd_rh              (ipy) = 0.0 
         cgrid%mmsqu_nep                 (ipy) = 0.0 
         cgrid%mmsqu_rlongup             (ipy) = 0.0 
         cgrid%mmsqu_parup               (ipy) = 0.0 
         cgrid%mmsqu_nirup               (ipy) = 0.0 
         cgrid%mmsqu_rshortup            (ipy) = 0.0 
         cgrid%mmsqu_rnet                (ipy) = 0.0 
         cgrid%mmsqu_albedo              (ipy) = 0.0 
         cgrid%mmsqu_ustar               (ipy) = 0.0 
         cgrid%mmsqu_carbon_ac           (ipy) = 0.0 
         cgrid%mmsqu_carbon_st           (ipy) = 0.0 
         cgrid%mmsqu_vapor_gc            (ipy) = 0.0 
         cgrid%mmsqu_vapor_ac            (ipy) = 0.0 
         cgrid%mmsqu_sensible_gc         (ipy) = 0.0 
         cgrid%mmsqu_sensible_ac         (ipy) = 0.0 


         !---------------------------------------------------------------------------------!
         !       Loop over sites.                                                          !
         !---------------------------------------------------------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            cpoly%mmean_atm_theiv      (isi) = 0.0
            cpoly%mmean_atm_theta      (isi) = 0.0
            cpoly%mmean_atm_temp       (isi) = 0.0
            cpoly%mmean_atm_vpdef      (isi) = 0.0
            cpoly%mmean_atm_shv        (isi) = 0.0
            cpoly%mmean_atm_rshort     (isi) = 0.0
            cpoly%mmean_atm_rshort_diff(isi) = 0.0
            cpoly%mmean_atm_par        (isi) = 0.0
            cpoly%mmean_atm_par_diff   (isi) = 0.0
            cpoly%mmean_atm_rlong      (isi) = 0.0
            cpoly%mmean_atm_vels       (isi) = 0.0
            cpoly%mmean_atm_rhos       (isi) = 0.0
            cpoly%mmean_atm_prss       (isi) = 0.0
            cpoly%mmean_atm_co2        (isi) = 0.0
            cpoly%mmean_pcpg           (isi) = 0.0
            cpoly%mmean_qpcpg          (isi) = 0.0
            cpoly%mmean_dpcpg          (isi) = 0.0


            !------------------------------------------------------------------------------!
            !       Loop over sites.                                                       !
            !------------------------------------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch=> csite%patch(ipa)

               csite%mmean_fast_soil_c      (ipa) = 0.0
               csite%mmean_slow_soil_c      (ipa) = 0.0
               csite%mmean_struct_soil_c    (ipa) = 0.0
               csite%mmean_struct_soil_l    (ipa) = 0.0
               csite%mmean_fast_soil_n      (ipa) = 0.0
               csite%mmean_mineral_soil_n   (ipa) = 0.0
               csite%mmean_co2_residual     (ipa) = 0.0
               csite%mmean_energy_residual  (ipa) = 0.0
               csite%mmean_water_residual   (ipa) = 0.0
               csite%mmean_rh               (ipa) = 0.0
               csite%mmean_cwd_rh           (ipa) = 0.0
               csite%mmean_nep              (ipa) = 0.0
               csite%mmean_A_decomp         (ipa) = 0.0
               csite%mmean_Af_decomp        (ipa) = 0.0
               csite%mmean_rk4step          (ipa) = 0.0
               csite%mmean_available_water  (ipa) = 0.0
               csite%mmean_can_theiv        (ipa) = 0.0
               csite%mmean_can_theta        (ipa) = 0.0
               csite%mmean_can_vpdef        (ipa) = 0.0
               csite%mmean_can_temp         (ipa) = 0.0
               csite%mmean_can_shv          (ipa) = 0.0
               csite%mmean_can_co2          (ipa) = 0.0
               csite%mmean_can_rhos         (ipa) = 0.0
               csite%mmean_can_prss         (ipa) = 0.0
               csite%mmean_gnd_temp         (ipa) = 0.0
               csite%mmean_gnd_shv          (ipa) = 0.0
               csite%mmean_can_ggnd         (ipa) = 0.0
               csite%mmean_sfcw_depth       (ipa) = 0.0
               csite%mmean_sfcw_energy      (ipa) = 0.0
               csite%mmean_sfcw_mass        (ipa) = 0.0
               csite%mmean_sfcw_temp        (ipa) = 0.0
               csite%mmean_sfcw_fliq        (ipa) = 0.0
               csite%mmean_soil_energy    (:,ipa) = 0.0
               csite%mmean_soil_mstpot    (:,ipa) = 0.0
               csite%mmean_soil_water     (:,ipa) = 0.0
               csite%mmean_soil_temp      (:,ipa) = 0.0
               csite%mmean_soil_fliq      (:,ipa) = 0.0
               csite%mmean_rshort_gnd       (ipa) = 0.0
               csite%mmean_par_gnd          (ipa) = 0.0
               csite%mmean_rlong_gnd        (ipa) = 0.0
               csite%mmean_rlongup          (ipa) = 0.0
               csite%mmean_parup            (ipa) = 0.0
               csite%mmean_nirup            (ipa) = 0.0
               csite%mmean_rshortup         (ipa) = 0.0
               csite%mmean_rnet             (ipa) = 0.0
               csite%mmean_albedo           (ipa) = 0.0
               csite%mmean_albedo_par       (ipa) = 0.0
               csite%mmean_albedo_nir       (ipa) = 0.0
               csite%mmean_rlong_albedo     (ipa) = 0.0
               csite%mmean_ustar            (ipa) = 0.0
               csite%mmean_tstar            (ipa) = 0.0
               csite%mmean_qstar            (ipa) = 0.0
               csite%mmean_cstar            (ipa) = 0.0
               csite%mmean_carbon_ac        (ipa) = 0.0
               csite%mmean_carbon_st        (ipa) = 0.0
               csite%mmean_vapor_gc         (ipa) = 0.0
               csite%mmean_vapor_ac         (ipa) = 0.0
               csite%mmean_smoist_gg      (:,ipa) = 0.0
               csite%mmean_throughfall      (ipa) = 0.0
               csite%mmean_transloss      (:,ipa) = 0.0
               csite%mmean_runoff           (ipa) = 0.0
               csite%mmean_drainage         (ipa) = 0.0
               csite%mmean_sensible_gc      (ipa) = 0.0
               csite%mmean_sensible_ac      (ipa) = 0.0
               csite%mmean_sensible_gg    (:,ipa) = 0.0
               csite%mmean_qthroughfall     (ipa) = 0.0
               csite%mmean_qrunoff          (ipa) = 0.0
               csite%mmean_qdrainage        (ipa) = 0.0
               csite%mmean_A_decomp         (ipa) = 0.0
               csite%mmean_Af_decomp        (ipa) = 0.0
               csite%mmean_co2_residual     (ipa) = 0.0
               csite%mmean_energy_residual  (ipa) = 0.0
               csite%mmean_water_residual   (ipa) = 0.0
               csite%mmsqu_rh               (ipa) = 0.0
               csite%mmsqu_cwd_rh           (ipa) = 0.0
               csite%mmsqu_nep              (ipa) = 0.0
               csite%mmsqu_rlongup          (ipa) = 0.0
               csite%mmsqu_parup            (ipa) = 0.0
               csite%mmsqu_nirup            (ipa) = 0.0
               csite%mmsqu_rshortup         (ipa) = 0.0
               csite%mmsqu_rnet             (ipa) = 0.0
               csite%mmsqu_albedo           (ipa) = 0.0
               csite%mmsqu_ustar            (ipa) = 0.0
               csite%mmsqu_carbon_ac        (ipa) = 0.0
               csite%mmsqu_carbon_st        (ipa) = 0.0
               csite%mmsqu_vapor_gc         (ipa) = 0.0
               csite%mmsqu_vapor_ac         (ipa) = 0.0
               csite%mmsqu_sensible_gc      (ipa) = 0.0
               csite%mmsqu_sensible_ac      (ipa) = 0.0


               !---------------------------------------------------------------------------!
               !       Loop over cohorts.                                                  !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts
                  cpatch%mmean_lai               (ico) = 0.0
                  cpatch%mmean_bleaf             (ico) = 0.0
                  cpatch%mmean_broot             (ico) = 0.0
                  cpatch%mmean_bstorage          (ico) = 0.0
                  cpatch%mmean_mort_rate       (:,ico) = 0.0
                  cpatch%mmean_leaf_maintenance  (ico) = 0.0
                  cpatch%mmean_root_maintenance  (ico) = 0.0
                  cpatch%mmean_leaf_drop         (ico) = 0.0
                  cpatch%mmean_cb                (ico) = 0.0
                  cpatch%mmean_gpp               (ico) = 0.0
                  cpatch%mmean_npp               (ico) = 0.0
                  cpatch%mmean_leaf_resp         (ico) = 0.0
                  cpatch%mmean_root_resp         (ico) = 0.0
                  cpatch%mmean_growth_resp       (ico) = 0.0
                  cpatch%mmean_storage_resp      (ico) = 0.0
                  cpatch%mmean_vleaf_resp        (ico) = 0.0
                  cpatch%mmean_plresp            (ico) = 0.0
                  cpatch%mmean_leaf_energy       (ico) = 0.0
                  cpatch%mmean_leaf_water        (ico) = 0.0
                  cpatch%mmean_leaf_hcap         (ico) = 0.0
                  cpatch%mmean_leaf_vpdef        (ico) = 0.0
                  cpatch%mmean_leaf_temp         (ico) = 0.0
                  cpatch%mmean_leaf_fliq         (ico) = 0.0
                  cpatch%mmean_leaf_gsw          (ico) = 0.0
                  cpatch%mmean_leaf_gbw          (ico) = 0.0
                  cpatch%mmean_wood_energy       (ico) = 0.0
                  cpatch%mmean_wood_water        (ico) = 0.0
                  cpatch%mmean_wood_hcap         (ico) = 0.0
                  cpatch%mmean_wood_temp         (ico) = 0.0
                  cpatch%mmean_wood_fliq         (ico) = 0.0
                  cpatch%mmean_wood_gbw          (ico) = 0.0
                  cpatch%mmean_fs_open           (ico) = 0.0
                  cpatch%mmean_fsw               (ico) = 0.0
                  cpatch%mmean_fsn               (ico) = 0.0
                  cpatch%mmean_a_light           (ico) = 0.0
                  cpatch%mmean_a_rubp            (ico) = 0.0
                  cpatch%mmean_a_co2             (ico) = 0.0
                  cpatch%mmean_psi_open          (ico) = 0.0
                  cpatch%mmean_psi_closed        (ico) = 0.0
                  cpatch%mmean_water_supply      (ico) = 0.0
                  cpatch%mmean_light_level       (ico) = 0.0
                  cpatch%mmean_light_level_beam  (ico) = 0.0
                  cpatch%mmean_light_level_diff  (ico) = 0.0

                  cpatch%mmean_par_level_beam    (ico) = 0.0
                  cpatch%mmean_par_level_diffd   (ico) = 0.0
                  cpatch%mmean_par_level_diffu   (ico) = 0.0

                  cpatch%mmean_par_l             (ico) = 0.0
                  cpatch%mmean_par_l_beam        (ico) = 0.0
                  cpatch%mmean_par_l_diff        (ico) = 0.0
                  cpatch%mmean_rshort_l          (ico) = 0.0
                  cpatch%mmean_rlong_l           (ico) = 0.0
                  cpatch%mmean_sensible_lc       (ico) = 0.0
                  cpatch%mmean_vapor_lc          (ico) = 0.0
                  cpatch%mmean_transp            (ico) = 0.0
                  cpatch%mmean_intercepted_al    (ico) = 0.0
                  cpatch%mmean_wshed_lg          (ico) = 0.0
                  cpatch%mmean_rshort_w          (ico) = 0.0
                  cpatch%mmean_rlong_w           (ico) = 0.0
                  cpatch%mmean_rad_profile     (:,ico) = 0.0
                  cpatch%mmean_sensible_wc       (ico) = 0.0
                  cpatch%mmean_vapor_wc          (ico) = 0.0
                  cpatch%mmean_intercepted_aw    (ico) = 0.0
                  cpatch%mmean_wshed_wg          (ico) = 0.0
                  cpatch%mmean_nppleaf           (ico) = 0.0
                  cpatch%mmean_nppfroot          (ico) = 0.0
                  cpatch%mmean_nppsapwood        (ico) = 0.0
                  cpatch%mmean_nppcroot          (ico) = 0.0
                  cpatch%mmean_nppseeds          (ico) = 0.0
                  cpatch%mmean_nppwood           (ico) = 0.0
                  cpatch%mmean_nppdaily          (ico) = 0.0
                  cpatch%mmsqu_gpp               (ico) = 0.0
                  cpatch%mmsqu_npp               (ico) = 0.0
                  cpatch%mmsqu_plresp            (ico) = 0.0
                  cpatch%mmsqu_sensible_lc       (ico) = 0.0
                  cpatch%mmsqu_vapor_lc          (ico) = 0.0
                  cpatch%mmsqu_transp            (ico) = 0.0
                  cpatch%mmsqu_sensible_wc       (ico) = 0.0
                  cpatch%mmsqu_vapor_wc          (ico) = 0.0
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine zero_ed_mmean_vars
   !=======================================================================================!
   !=======================================================================================!










   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !                             |-------------------------------|                         !
   !                             |**** MEAN DIEL SUBROUTINES ****|                         !
   !                             |-------------------------------|                         !
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine integrates most of the mean diel variables.  This subroutine is    !
   ! called after the fmean variables have been normalised, so we can take advantage of    !
   ! their values.                                                                         !
   !                                                                                       !
   !    A few variables are _NOT_ integrated here: quantities such as temperature and      !
   ! liquid fraction are found after the mean diel of the extensive variables have been    !
   ! normalised.                                                                           !
   !---------------------------------------------------------------------------------------!
   subroutine integrate_ed_qmean_vars(cgrid)
      use ed_state_vars, only : edtype              & ! structure
                              , polygontype         & ! structure
                              , sitetype            & ! structure
                              , patchtype           ! ! structure
      use ed_misc_coms , only : frqfast             & ! intent(in)
                              , ndcycle             & ! intent(in)
                              , current_time        & ! intent(in)
                              , simtime             ! ! structure
      use consts_coms  , only : day_sec             ! ! intent(in)
      implicit none
      !----- Argument ---------------------------------------------------------------------!
      type(edtype)      , target  :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype) , pointer :: cpoly
      type(sitetype)    , pointer :: csite
      type(patchtype)   , pointer :: cpatch
      type(simtime)               :: daybefore
      integer                     :: ipy
      integer                     :: isi
      integer                     :: ipa
      integer                     :: ico
      integer                     :: t
      real                        :: ndaysi
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the index corresponding to this time of the day for the mean diurnal      !
      ! cycle averages.                                                                    !
      !------------------------------------------------------------------------------------!
      t = ceiling(mod(current_time%time,day_sec)/frqfast)
      if (t == 0) t = nint(day_sec/frqfast)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find which day we have just integrated, we will use it to determine the right  !
      ! scaling factor.                                                                    !
      !------------------------------------------------------------------------------------!
      call yesterday_info(current_time,daybefore,ndaysi)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Use the variables that have been already aggregated.                          !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)


         !---------------------------------------------------------------------------------!
         !    Integrate polygon-level variables.                                           !
         !---------------------------------------------------------------------------------!
         cgrid%qmean_gpp              (t,ipy) = cgrid%qmean_gpp              (t,ipy)       &
                                              + cgrid%fmean_gpp                (ipy)       &
                                              * ndaysi
         cgrid%qmean_npp              (t,ipy) = cgrid%qmean_npp              (t,ipy)       &
                                              + cgrid%fmean_npp                (ipy)       &
                                              * ndaysi
         cgrid%qmean_leaf_resp        (t,ipy) = cgrid%qmean_leaf_resp        (t,ipy)       &
                                              + cgrid%fmean_leaf_resp          (ipy)       &
                                              * ndaysi
         cgrid%qmean_root_resp        (t,ipy) = cgrid%qmean_root_resp        (t,ipy)       &
                                              + cgrid%fmean_root_resp          (ipy)       &
                                              * ndaysi
         cgrid%qmean_growth_resp      (t,ipy) = cgrid%qmean_growth_resp      (t,ipy)       &
                                              + cgrid%fmean_growth_resp        (ipy)       &
                                              * ndaysi
         cgrid%qmean_storage_resp     (t,ipy) = cgrid%qmean_storage_resp     (t,ipy)       &
                                              + cgrid%fmean_storage_resp       (ipy)       &
                                              * ndaysi
         cgrid%qmean_vleaf_resp       (t,ipy) = cgrid%qmean_vleaf_resp       (t,ipy)       &
                                              + cgrid%fmean_vleaf_resp         (ipy)       &
                                              * ndaysi
         cgrid%qmean_plresp           (t,ipy) = cgrid%qmean_plresp           (t,ipy)       &
                                              + cgrid%fmean_plresp             (ipy)       &
                                              * ndaysi
         cgrid%qmean_leaf_energy      (t,ipy) = cgrid%qmean_leaf_energy      (t,ipy)       &
                                              + cgrid%fmean_leaf_energy        (ipy)       &
                                              * ndaysi
         cgrid%qmean_leaf_water       (t,ipy) = cgrid%qmean_leaf_water       (t,ipy)       &
                                              + cgrid%fmean_leaf_water         (ipy)       &
                                              * ndaysi
         cgrid%qmean_leaf_hcap        (t,ipy) = cgrid%qmean_leaf_hcap        (t,ipy)       &
                                              + cgrid%fmean_leaf_hcap          (ipy)       &
                                              * ndaysi
         cgrid%qmean_leaf_vpdef       (t,ipy) = cgrid%qmean_leaf_vpdef       (t,ipy)       &
                                              + cgrid%fmean_leaf_vpdef         (ipy)       &
                                              * ndaysi
         cgrid%qmean_leaf_gsw         (t,ipy) = cgrid%qmean_leaf_gsw         (t,ipy)       &
                                              + cgrid%fmean_leaf_gsw           (ipy)       &
                                              * ndaysi
         cgrid%qmean_leaf_gbw         (t,ipy) = cgrid%qmean_leaf_gbw         (t,ipy)       &
                                              + cgrid%fmean_leaf_gbw           (ipy)       &
                                              * ndaysi
         cgrid%qmean_wood_energy      (t,ipy) = cgrid%qmean_wood_energy      (t,ipy)       &
                                              + cgrid%fmean_wood_energy        (ipy)       &
                                              * ndaysi
         cgrid%qmean_wood_water       (t,ipy) = cgrid%qmean_wood_water       (t,ipy)       &
                                              + cgrid%fmean_wood_water         (ipy)       &
                                              * ndaysi
         cgrid%qmean_wood_hcap        (t,ipy) = cgrid%qmean_wood_hcap        (t,ipy)       &
                                              + cgrid%fmean_wood_hcap          (ipy)       &
                                              * ndaysi
         cgrid%qmean_wood_gbw         (t,ipy) = cgrid%qmean_wood_gbw         (t,ipy)       &
                                              + cgrid%fmean_wood_gbw           (ipy)       &
                                              * ndaysi
         cgrid%qmean_fs_open          (t,ipy) = cgrid%qmean_fs_open          (t,ipy)       &
                                              + cgrid%fmean_fs_open            (ipy)       &
                                              * ndaysi
         cgrid%qmean_fsw              (t,ipy) = cgrid%qmean_fsw              (t,ipy)       &
                                              + cgrid%fmean_fsw                (ipy)       &
                                              * ndaysi
         cgrid%qmean_fsn              (t,ipy) = cgrid%qmean_fsn              (t,ipy)       &
                                              + cgrid%fmean_fsn                (ipy)       &
                                              * ndaysi
         cgrid%qmean_a_light          (t,ipy) = cgrid%qmean_a_light          (t,ipy)       &
                                              + cgrid%fmean_a_light            (ipy)       &
                                              * ndaysi
         cgrid%qmean_a_rubp           (t,ipy) = cgrid%qmean_a_rubp           (t,ipy)       &
                                              + cgrid%fmean_a_rubp             (ipy)       &
                                              * ndaysi
         cgrid%qmean_a_co2            (t,ipy) = cgrid%qmean_a_co2            (t,ipy)       &
                                              + cgrid%fmean_a_co2              (ipy)       &
                                              * ndaysi
         cgrid%qmean_psi_open         (t,ipy) = cgrid%qmean_psi_open         (t,ipy)       &
                                              + cgrid%fmean_psi_open           (ipy)       &
                                              * ndaysi
         cgrid%qmean_psi_closed       (t,ipy) = cgrid%qmean_psi_closed       (t,ipy)       &
                                              + cgrid%fmean_psi_closed         (ipy)       &
                                              * ndaysi
         cgrid%qmean_water_supply     (t,ipy) = cgrid%qmean_water_supply     (t,ipy)       &
                                              + cgrid%fmean_water_supply       (ipy)       &
                                              * ndaysi
         cgrid%qmean_par_l            (t,ipy) = cgrid%qmean_par_l            (t,ipy)       &
                                              + cgrid%fmean_par_l              (ipy)       &
                                              * ndaysi
         cgrid%qmean_par_l_beam       (t,ipy) = cgrid%qmean_par_l_beam       (t,ipy)       &
                                              + cgrid%fmean_par_l_beam         (ipy)       &
                                              * ndaysi
         cgrid%qmean_par_l_diff       (t,ipy) = cgrid%qmean_par_l_diff       (t,ipy)       &
                                              + cgrid%fmean_par_l_diff         (ipy)       &
                                              * ndaysi
         cgrid%qmean_rshort_l         (t,ipy) = cgrid%qmean_rshort_l         (t,ipy)       &
                                              + cgrid%fmean_rshort_l           (ipy)       &
                                              * ndaysi
         cgrid%qmean_rlong_l          (t,ipy) = cgrid%qmean_rlong_l          (t,ipy)       &
                                              + cgrid%fmean_rlong_l            (ipy)       &
                                              * ndaysi
         cgrid%qmean_sensible_lc      (t,ipy) = cgrid%qmean_sensible_lc      (t,ipy)       &
                                              + cgrid%fmean_sensible_lc        (ipy)       &
                                              * ndaysi
         cgrid%qmean_vapor_lc         (t,ipy) = cgrid%qmean_vapor_lc         (t,ipy)       &
                                              + cgrid%fmean_vapor_lc           (ipy)       &
                                              * ndaysi
         cgrid%qmean_transp           (t,ipy) = cgrid%qmean_transp           (t,ipy)       &
                                              + cgrid%fmean_transp             (ipy)       &
                                              * ndaysi
         cgrid%qmean_intercepted_al   (t,ipy) = cgrid%qmean_intercepted_al   (t,ipy)       &
                                              + cgrid%fmean_intercepted_al     (ipy)       &
                                              * ndaysi
         cgrid%qmean_wshed_lg         (t,ipy) = cgrid%qmean_wshed_lg         (t,ipy)       &
                                              + cgrid%fmean_wshed_lg           (ipy)       &
                                              * ndaysi
         cgrid%qmean_rshort_w         (t,ipy) = cgrid%qmean_rshort_w         (t,ipy)       &
                                              + cgrid%fmean_rshort_w           (ipy)       &
                                              * ndaysi
         cgrid%qmean_rlong_w          (t,ipy) = cgrid%qmean_rlong_w          (t,ipy)       &
                                              + cgrid%fmean_rlong_w            (ipy)       &
                                              * ndaysi
         cgrid%qmean_sensible_wc      (t,ipy) = cgrid%qmean_sensible_wc      (t,ipy)       &
                                              + cgrid%fmean_sensible_wc        (ipy)       &
                                              * ndaysi
         cgrid%qmean_vapor_wc         (t,ipy) = cgrid%qmean_vapor_wc         (t,ipy)       &
                                              + cgrid%fmean_vapor_wc           (ipy)       &
                                              * ndaysi
         cgrid%qmean_intercepted_aw   (t,ipy) = cgrid%qmean_intercepted_aw   (t,ipy)       &
                                              + cgrid%fmean_intercepted_aw     (ipy)       &
                                              * ndaysi
         cgrid%qmean_wshed_wg         (t,ipy) = cgrid%qmean_wshed_wg         (t,ipy)       &
                                              + cgrid%fmean_wshed_wg           (ipy)       &
                                              * ndaysi
         cgrid%qmean_rh               (t,ipy) = cgrid%qmean_rh               (t,ipy)       &
                                              + cgrid%fmean_rh                 (ipy)       &
                                              * ndaysi
         cgrid%qmean_cwd_rh           (t,ipy) = cgrid%qmean_cwd_rh           (t,ipy)       &
                                              + cgrid%fmean_cwd_rh             (ipy)       &
                                              * ndaysi
         cgrid%qmean_nep              (t,ipy) = cgrid%qmean_nep              (t,ipy)       &
                                              + cgrid%fmean_nep                (ipy)       &
                                              * ndaysi
         cgrid%qmean_rk4step          (t,ipy) = cgrid%qmean_rk4step          (t,ipy)       &
                                              + cgrid%fmean_rk4step            (ipy)       &
                                              * ndaysi
         cgrid%qmean_available_water  (t,ipy) = cgrid%qmean_available_water  (t,ipy)       &
                                              + cgrid%fmean_available_water    (ipy)       &
                                              * ndaysi
         cgrid%qmean_can_theiv        (t,ipy) = cgrid%qmean_can_theiv        (t,ipy)       &
                                              + cgrid%fmean_can_theiv          (ipy)       &
                                              * ndaysi
         cgrid%qmean_can_theta        (t,ipy) = cgrid%qmean_can_theta        (t,ipy)       &
                                              + cgrid%fmean_can_theta          (ipy)       &
                                              * ndaysi
         cgrid%qmean_can_vpdef        (t,ipy) = cgrid%qmean_can_vpdef        (t,ipy)       &
                                              + cgrid%fmean_can_vpdef          (ipy)       &
                                              * ndaysi
         cgrid%qmean_can_shv          (t,ipy) = cgrid%qmean_can_shv          (t,ipy)       &
                                              + cgrid%fmean_can_shv            (ipy)       &
                                              * ndaysi
         cgrid%qmean_can_co2          (t,ipy) = cgrid%qmean_can_co2          (t,ipy)       &
                                              + cgrid%fmean_can_co2            (ipy)       &
                                              * ndaysi
         cgrid%qmean_can_prss         (t,ipy) = cgrid%qmean_can_prss         (t,ipy)       &
                                              + cgrid%fmean_can_prss           (ipy)       &
                                              * ndaysi
         cgrid%qmean_gnd_temp         (t,ipy) = cgrid%qmean_gnd_temp         (t,ipy)       &
                                              + cgrid%fmean_gnd_temp           (ipy)       &
                                              * ndaysi
         cgrid%qmean_gnd_shv          (t,ipy) = cgrid%qmean_gnd_shv          (t,ipy)       &
                                              + cgrid%fmean_gnd_shv            (ipy)       &
                                              * ndaysi
         cgrid%qmean_can_ggnd         (t,ipy) = cgrid%qmean_can_ggnd         (t,ipy)       &
                                              + cgrid%fmean_can_ggnd           (ipy)       &
                                              * ndaysi
         cgrid%qmean_sfcw_depth       (t,ipy) = cgrid%qmean_sfcw_depth       (t,ipy)       &
                                              + cgrid%fmean_sfcw_depth         (ipy)       &
                                              * ndaysi
         !----- During the integration, pounding internal energy must be in J/m2. ---------!
         cgrid%qmean_sfcw_energy      (t,ipy) = cgrid%qmean_sfcw_energy      (t,ipy)       &
                                              + cgrid%fmean_sfcw_energy        (ipy)       &
                                              * cgrid%fmean_sfcw_mass          (ipy)       &
                                              * ndaysi
         cgrid%qmean_sfcw_mass        (t,ipy) = cgrid%qmean_sfcw_mass        (t,ipy)       &
                                              + cgrid%fmean_sfcw_mass          (ipy)       &
                                              * ndaysi
         cgrid%qmean_soil_energy    (:,t,ipy) = cgrid%qmean_soil_energy    (:,t,ipy)       &
                                              + cgrid%fmean_soil_energy      (:,ipy)       &
                                              * ndaysi
         cgrid%qmean_soil_mstpot    (:,t,ipy) = cgrid%qmean_soil_mstpot    (:,t,ipy)       &
                                              + cgrid%fmean_soil_mstpot      (:,ipy)       &
                                              * ndaysi
         cgrid%qmean_soil_water     (:,t,ipy) = cgrid%qmean_soil_water     (:,t,ipy)       &
                                              + cgrid%fmean_soil_water       (:,ipy)       &
                                              * ndaysi
         cgrid%qmean_rshort_gnd       (t,ipy) = cgrid%qmean_rshort_gnd       (t,ipy)       &
                                              + cgrid%fmean_rshort_gnd         (ipy)       &
                                              * ndaysi
         cgrid%qmean_par_gnd          (t,ipy) = cgrid%qmean_par_gnd          (t,ipy)       &
                                              + cgrid%fmean_par_gnd            (ipy)       &
                                              * ndaysi
         cgrid%qmean_rlong_gnd        (t,ipy) = cgrid%qmean_rlong_gnd        (t,ipy)       &
                                              + cgrid%fmean_rlong_gnd          (ipy)       &
                                              * ndaysi
         cgrid%qmean_rlongup          (t,ipy) = cgrid%qmean_rlongup          (t,ipy)       &
                                              + cgrid%fmean_rlongup            (ipy)       &
                                              * ndaysi
         cgrid%qmean_parup            (t,ipy) = cgrid%qmean_parup            (t,ipy)       &
                                              + cgrid%fmean_parup              (ipy)       &
                                              * ndaysi
         cgrid%qmean_nirup            (t,ipy) = cgrid%qmean_nirup            (t,ipy)       &
                                              + cgrid%fmean_nirup              (ipy)       &
                                              * ndaysi
         cgrid%qmean_rshortup         (t,ipy) = cgrid%qmean_rshortup         (t,ipy)       &
                                              + cgrid%fmean_rshortup           (ipy)       &
                                              * ndaysi
         cgrid%qmean_rnet             (t,ipy) = cgrid%qmean_rnet             (t,ipy)       &
                                              + cgrid%fmean_rnet               (ipy)       &
                                              * ndaysi
         cgrid%qmean_albedo           (t,ipy) = cgrid%qmean_albedo           (t,ipy)       &
                                              + cgrid%fmean_albedo             (ipy)       &
                                              * ndaysi
         cgrid%qmean_albedo_par       (t,ipy) = cgrid%qmean_albedo_par       (t,ipy)       &
                                              + cgrid%fmean_albedo_par         (ipy)       &
                                              * ndaysi
         cgrid%qmean_albedo_nir       (t,ipy) = cgrid%qmean_albedo_nir       (t,ipy)       &
                                              + cgrid%fmean_albedo_nir         (ipy)       &
                                              * ndaysi
         cgrid%qmean_rlong_albedo     (t,ipy) = cgrid%qmean_rlong_albedo     (t,ipy)       &
                                              + cgrid%fmean_rlong_albedo       (ipy)       &
                                              * ndaysi
         cgrid%qmean_ustar            (t,ipy) = cgrid%qmean_ustar            (t,ipy)       &
                                              + cgrid%fmean_ustar              (ipy)       &
                                              * ndaysi
         cgrid%qmean_tstar            (t,ipy) = cgrid%qmean_tstar            (t,ipy)       &
                                              + cgrid%fmean_tstar              (ipy)       &
                                              * ndaysi
         cgrid%qmean_qstar            (t,ipy) = cgrid%qmean_qstar            (t,ipy)       &
                                              + cgrid%fmean_qstar              (ipy)       &
                                              * ndaysi
         cgrid%qmean_cstar            (t,ipy) = cgrid%qmean_cstar            (t,ipy)       &
                                              + cgrid%fmean_cstar              (ipy)       &
                                              * ndaysi
         cgrid%qmean_carbon_ac        (t,ipy) = cgrid%qmean_carbon_ac        (t,ipy)       &
                                              + cgrid%fmean_carbon_ac          (ipy)       &
                                              * ndaysi
         cgrid%qmean_carbon_st        (t,ipy) = cgrid%qmean_carbon_st        (t,ipy)       &
                                              + cgrid%fmean_carbon_st          (ipy)       &
                                              * ndaysi
         cgrid%qmean_vapor_gc         (t,ipy) = cgrid%qmean_vapor_gc         (t,ipy)       &
                                              + cgrid%fmean_vapor_gc           (ipy)       &
                                              * ndaysi
         cgrid%qmean_vapor_ac         (t,ipy) = cgrid%qmean_vapor_ac         (t,ipy)       &
                                              + cgrid%fmean_vapor_ac           (ipy)       &
                                              * ndaysi
         cgrid%qmean_smoist_gg      (:,t,ipy) = cgrid%qmean_smoist_gg      (:,t,ipy)       &
                                              + cgrid%fmean_smoist_gg        (:,ipy)       &
                                              * ndaysi
         cgrid%qmean_throughfall      (t,ipy) = cgrid%qmean_throughfall      (t,ipy)       &
                                              + cgrid%fmean_throughfall        (ipy)       &
                                              * ndaysi
         cgrid%qmean_transloss      (:,t,ipy) = cgrid%qmean_transloss      (:,t,ipy)       &
                                              + cgrid%fmean_transloss        (:,ipy)       &
                                              * ndaysi
         cgrid%qmean_runoff           (t,ipy) = cgrid%qmean_runoff           (t,ipy)       &
                                              + cgrid%fmean_runoff             (ipy)       &
                                              * ndaysi
         cgrid%qmean_drainage         (t,ipy) = cgrid%qmean_drainage         (t,ipy)       &
                                              + cgrid%fmean_drainage           (ipy)       &
                                              * ndaysi
         cgrid%qmean_sensible_gc      (t,ipy) = cgrid%qmean_sensible_gc      (t,ipy)       &
                                              + cgrid%fmean_sensible_gc        (ipy)       &
                                              * ndaysi
         cgrid%qmean_sensible_ac      (t,ipy) = cgrid%qmean_sensible_ac      (t,ipy)       &
                                              + cgrid%fmean_sensible_ac        (ipy)       &
                                              * ndaysi
         cgrid%qmean_sensible_gg    (:,t,ipy) = cgrid%qmean_sensible_gg    (:,t,ipy)       &
                                              + cgrid%fmean_sensible_gg      (:,ipy)       &
                                              * ndaysi
         cgrid%qmean_qthroughfall     (t,ipy) = cgrid%qmean_qthroughfall     (t,ipy)       &
                                              + cgrid%fmean_qthroughfall       (ipy)       &
                                              * ndaysi
         cgrid%qmean_qrunoff          (t,ipy) = cgrid%qmean_qrunoff          (t,ipy)       &
                                              + cgrid%fmean_qrunoff            (ipy)       &
                                              * ndaysi
         cgrid%qmean_qdrainage        (t,ipy) = cgrid%qmean_qdrainage        (t,ipy)       &
                                              + cgrid%fmean_qdrainage          (ipy)       &
                                              * ndaysi
         cgrid%qmean_atm_theiv        (t,ipy) = cgrid%qmean_atm_theiv        (t,ipy)       &
                                              + cgrid%fmean_atm_theiv          (ipy)       &
                                              * ndaysi
         cgrid%qmean_atm_theta        (t,ipy) = cgrid%qmean_atm_theta        (t,ipy)       &
                                              + cgrid%fmean_atm_theta          (ipy)       &
                                              * ndaysi
         cgrid%qmean_atm_vpdef        (t,ipy) = cgrid%qmean_atm_vpdef        (t,ipy)       &
                                              + cgrid%fmean_atm_vpdef          (ipy)       &
                                              * ndaysi
         cgrid%qmean_atm_shv          (t,ipy) = cgrid%qmean_atm_shv          (t,ipy)       &
                                              + cgrid%fmean_atm_shv            (ipy)       &
                                              * ndaysi
         cgrid%qmean_atm_rshort       (t,ipy) = cgrid%qmean_atm_rshort       (t,ipy)       &
                                              + cgrid%fmean_atm_rshort         (ipy)       &
                                              * ndaysi
         cgrid%qmean_atm_rshort_diff  (t,ipy) = cgrid%qmean_atm_rshort_diff  (t,ipy)       &
                                              + cgrid%fmean_atm_rshort_diff    (ipy)       &
                                              * ndaysi
         cgrid%qmean_atm_par          (t,ipy) = cgrid%qmean_atm_par          (t,ipy)       &
                                              + cgrid%fmean_atm_par            (ipy)       &
                                              * ndaysi
         cgrid%qmean_atm_par_diff     (t,ipy) = cgrid%qmean_atm_par_diff     (t,ipy)       &
                                              + cgrid%fmean_atm_par_diff       (ipy)       &
                                              * ndaysi
         cgrid%qmean_atm_rlong        (t,ipy) = cgrid%qmean_atm_rlong        (t,ipy)       &
                                              + cgrid%fmean_atm_rlong          (ipy)       &
                                              * ndaysi
         cgrid%qmean_atm_vels         (t,ipy) = cgrid%qmean_atm_vels         (t,ipy)       &
                                              + cgrid%fmean_atm_vels           (ipy)       &
                                              * ndaysi
         cgrid%qmean_atm_prss         (t,ipy) = cgrid%qmean_atm_prss         (t,ipy)       &
                                              + cgrid%fmean_atm_prss           (ipy)       &
                                              * ndaysi
         cgrid%qmean_atm_co2          (t,ipy) = cgrid%qmean_atm_co2          (t,ipy)       &
                                              + cgrid%fmean_atm_co2            (ipy)       &
                                              * ndaysi
         cgrid%qmean_pcpg             (t,ipy) = cgrid%qmean_pcpg             (t,ipy)       &
                                              + cgrid%fmean_pcpg               (ipy)       &
                                              * ndaysi
         cgrid%qmean_qpcpg            (t,ipy) = cgrid%qmean_qpcpg            (t,ipy)       &
                                              + cgrid%fmean_qpcpg              (ipy)       &
                                              * ndaysi
         cgrid%qmean_dpcpg            (t,ipy) = cgrid%qmean_dpcpg            (t,ipy)       &
                                              + cgrid%fmean_dpcpg              (ipy)       &
                                              * ndaysi
         !----- Mean sum of squares. ------------------------------------------------------#
         cgrid%qmsqu_gpp              (t,ipy) = cgrid%qmsqu_gpp                 (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_gpp          (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_npp              (t,ipy) = cgrid%qmsqu_npp                 (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_npp          (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_plresp           (t,ipy) = cgrid%qmsqu_plresp              (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_plresp       (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_sensible_lc      (t,ipy) = cgrid%qmsqu_sensible_lc         (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_sensible_lc  (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_vapor_lc         (t,ipy) = cgrid%qmsqu_vapor_lc            (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_vapor_lc     (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_transp           (t,ipy) = cgrid%qmsqu_transp              (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_transp       (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_sensible_wc      (t,ipy) = cgrid%qmsqu_sensible_wc         (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_sensible_wc  (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_vapor_wc         (t,ipy) = cgrid%qmsqu_vapor_wc            (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_vapor_wc     (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_rh               (t,ipy) = cgrid%qmsqu_rh                  (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_rh           (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_cwd_rh           (t,ipy) = cgrid%qmsqu_cwd_rh              (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_cwd_rh       (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_nep              (t,ipy) = cgrid%qmsqu_nep                 (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_nep          (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_rlongup          (t,ipy) = cgrid%qmsqu_rlongup             (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_rlongup      (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_parup            (t,ipy) = cgrid%qmsqu_parup               (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_parup        (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_nirup            (t,ipy) = cgrid%qmsqu_nirup               (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_nirup        (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_rshortup         (t,ipy) = cgrid%qmsqu_rshortup            (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_rshortup     (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_rnet             (t,ipy) = cgrid%qmsqu_rnet                (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_rnet         (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_albedo           (t,ipy) = cgrid%qmsqu_albedo              (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_albedo       (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_ustar            (t,ipy) = cgrid%qmsqu_ustar               (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_ustar        (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_carbon_ac        (t,ipy) = cgrid%qmsqu_carbon_ac           (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_carbon_ac    (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_carbon_st        (t,ipy) = cgrid%qmsqu_carbon_st           (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_carbon_st    (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_vapor_gc         (t,ipy) = cgrid%qmsqu_vapor_gc            (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_vapor_gc     (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_vapor_ac         (t,ipy) = cgrid%qmsqu_vapor_ac            (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_vapor_ac     (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_sensible_gc      (t,ipy) = cgrid%qmsqu_sensible_gc         (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_sensible_gc  (ipy))   &
                                              * ndaysi
         cgrid%qmsqu_sensible_ac      (t,ipy) = cgrid%qmsqu_sensible_ac         (t,ipy)    &
                                              + isqu_ftz(cgrid%fmean_sensible_ac  (ipy))   &
                                              * ndaysi
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Site loop.                                                                 !
         !---------------------------------------------------------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            !------------------------------------------------------------------------------!
            !    Integrate site-level variables.                                           !
            !------------------------------------------------------------------------------!
            cpoly%qmean_atm_theiv      (t,isi) = cpoly%qmean_atm_theiv      (t,isi)        &
                                               + cpoly%fmean_atm_theiv        (isi)        &
                                               * ndaysi
            cpoly%qmean_atm_theta      (t,isi) = cpoly%qmean_atm_theta      (t,isi)        &
                                               + cpoly%fmean_atm_theta        (isi)        &
                                               * ndaysi
            cpoly%qmean_atm_vpdef      (t,isi) = cpoly%qmean_atm_vpdef      (t,isi)        &
                                               + cpoly%fmean_atm_vpdef        (isi)        &
                                               * ndaysi
            cpoly%qmean_atm_shv        (t,isi) = cpoly%qmean_atm_shv        (t,isi)        &
                                               + cpoly%fmean_atm_shv          (isi)        &
                                               * ndaysi
            cpoly%qmean_atm_rshort     (t,isi) = cpoly%qmean_atm_rshort     (t,isi)        &
                                               + cpoly%fmean_atm_rshort       (isi)        &
                                               * ndaysi
            cpoly%qmean_atm_rshort_diff(t,isi) = cpoly%qmean_atm_rshort_diff(t,isi)        &
                                               + cpoly%fmean_atm_rshort_diff  (isi)        &
                                               * ndaysi
            cpoly%qmean_atm_par        (t,isi) = cpoly%qmean_atm_par        (t,isi)        &
                                               + cpoly%fmean_atm_par          (isi)        &
                                               * ndaysi
            cpoly%qmean_atm_par_diff   (t,isi) = cpoly%qmean_atm_par_diff   (t,isi)        &
                                               + cpoly%fmean_atm_par_diff     (isi)        &
                                               * ndaysi
            cpoly%qmean_atm_rlong      (t,isi) = cpoly%qmean_atm_rlong      (t,isi)        &
                                               + cpoly%fmean_atm_rlong        (isi)        &
                                               * ndaysi
            cpoly%qmean_atm_vels       (t,isi) = cpoly%qmean_atm_vels       (t,isi)        &
                                               + cpoly%fmean_atm_vels         (isi)        &
                                               * ndaysi
            cpoly%qmean_atm_prss       (t,isi) = cpoly%qmean_atm_prss       (t,isi)        &
                                               + cpoly%fmean_atm_prss         (isi)        &
                                               * ndaysi
            cpoly%qmean_atm_co2        (t,isi) = cpoly%qmean_atm_co2        (t,isi)        &
                                               + cpoly%fmean_atm_co2          (isi)        &
                                               * ndaysi
            cpoly%qmean_pcpg           (t,isi) = cpoly%qmean_pcpg           (t,isi)        &
                                               + cpoly%fmean_pcpg             (isi)        &
                                               * ndaysi
            cpoly%qmean_qpcpg          (t,isi) = cpoly%qmean_qpcpg          (t,isi)        &
                                               + cpoly%fmean_qpcpg            (isi)        &
                                               * ndaysi
            cpoly%qmean_dpcpg          (t,isi) = cpoly%qmean_dpcpg          (t,isi)        &
                                               + cpoly%fmean_dpcpg            (isi)        &
                                               * ndaysi
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Patch loop.                                                             !
            !------------------------------------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               !---------------------------------------------------------------------------!
               !      Integrate patch-level variables.                                     !
               !---------------------------------------------------------------------------!
               csite%qmean_rh               (t,ipa) = csite%qmean_rh               (t,ipa) &
                                                    + csite%fmean_rh                 (ipa) &
                                                    * ndaysi
               csite%qmean_cwd_rh           (t,ipa) = csite%qmean_cwd_rh           (t,ipa) &
                                                    + csite%fmean_cwd_rh             (ipa) &
                                                    * ndaysi
               csite%qmean_nep              (t,ipa) = csite%qmean_nep              (t,ipa) &
                                                    + csite%fmean_nep                (ipa) &
                                                    * ndaysi
               csite%qmean_rk4step          (t,ipa) = csite%qmean_rk4step          (t,ipa) &
                                                    + csite%fmean_rk4step            (ipa) &
                                                    * ndaysi
               csite%qmean_available_water  (t,ipa) = csite%qmean_available_water  (t,ipa) &
                                                    + csite%fmean_available_water    (ipa) &
                                                    * ndaysi
               csite%qmean_can_theiv        (t,ipa) = csite%qmean_can_theiv        (t,ipa) &
                                                    + csite%fmean_can_theiv          (ipa) &
                                                    * ndaysi
               csite%qmean_can_theta        (t,ipa) = csite%qmean_can_theta        (t,ipa) &
                                                    + csite%fmean_can_theta          (ipa) &
                                                    * ndaysi
               csite%qmean_can_vpdef        (t,ipa) = csite%qmean_can_vpdef        (t,ipa) &
                                                    + csite%fmean_can_vpdef          (ipa) &
                                                    * ndaysi
               csite%qmean_can_shv          (t,ipa) = csite%qmean_can_shv          (t,ipa) &
                                                    + csite%fmean_can_shv            (ipa) &
                                                    * ndaysi
               csite%qmean_can_co2          (t,ipa) = csite%qmean_can_co2          (t,ipa) &
                                                    + csite%fmean_can_co2            (ipa) &
                                                    * ndaysi
               csite%qmean_can_prss         (t,ipa) = csite%qmean_can_prss         (t,ipa) &
                                                    + csite%fmean_can_prss           (ipa) &
                                                    * ndaysi
               csite%qmean_gnd_temp         (t,ipa) = csite%qmean_gnd_temp         (t,ipa) &
                                                    + csite%fmean_gnd_temp           (ipa) &
                                                    * ndaysi
               csite%qmean_gnd_shv          (t,ipa) = csite%qmean_gnd_shv          (t,ipa) &
                                                    + csite%fmean_gnd_shv            (ipa) &
                                                    * ndaysi
               csite%qmean_can_ggnd         (t,ipa) = csite%qmean_can_ggnd         (t,ipa) &
                                                    + csite%fmean_can_ggnd           (ipa) &
                                                    * ndaysi
               csite%qmean_sfcw_depth       (t,ipa) = csite%qmean_sfcw_depth       (t,ipa) &
                                                    + csite%fmean_sfcw_depth         (ipa) &
                                                    * ndaysi
               !------ Integrate pounding internal energy in extensive format [J/m2]. -----!
               csite%qmean_sfcw_energy      (t,ipa) = csite%qmean_sfcw_energy      (t,ipa) &
                                                    + csite%fmean_sfcw_energy        (ipa) &
                                                    * csite%fmean_sfcw_mass          (ipa) &
                                                    * ndaysi
               csite%qmean_sfcw_mass        (t,ipa) = csite%qmean_sfcw_mass        (t,ipa) &
                                                    + csite%fmean_sfcw_mass          (ipa) &
                                                    * ndaysi
               csite%qmean_soil_energy    (:,t,ipa) = csite%qmean_soil_energy    (:,t,ipa) &
                                                    + csite%fmean_soil_energy      (:,ipa) &
                                                    * ndaysi
               csite%qmean_soil_mstpot    (:,t,ipa) = csite%qmean_soil_mstpot    (:,t,ipa) &
                                                    + csite%fmean_soil_mstpot      (:,ipa) &
                                                    * ndaysi
               csite%qmean_soil_water     (:,t,ipa) = csite%qmean_soil_water     (:,t,ipa) &
                                                    + csite%fmean_soil_water       (:,ipa) &
                                                    * ndaysi
               csite%qmean_rshort_gnd       (t,ipa) = csite%qmean_rshort_gnd       (t,ipa) &
                                                    + csite%fmean_rshort_gnd         (ipa) &
                                                    * ndaysi
               csite%qmean_par_gnd          (t,ipa) = csite%qmean_par_gnd          (t,ipa) &
                                                    + csite%fmean_par_gnd            (ipa) &
                                                    * ndaysi
               csite%qmean_rlong_gnd        (t,ipa) = csite%qmean_rlong_gnd        (t,ipa) &
                                                    + csite%fmean_rlong_gnd          (ipa) &
                                                    * ndaysi
               csite%qmean_rlongup          (t,ipa) = csite%qmean_rlongup          (t,ipa) &
                                                    + csite%fmean_rlongup            (ipa) &
                                                    * ndaysi
               csite%qmean_parup            (t,ipa) = csite%qmean_parup            (t,ipa) &
                                                    + csite%fmean_parup              (ipa) &
                                                    * ndaysi
               csite%qmean_nirup            (t,ipa) = csite%qmean_nirup            (t,ipa) &
                                                    + csite%fmean_nirup              (ipa) &
                                                    * ndaysi
               csite%qmean_rshortup         (t,ipa) = csite%qmean_rshortup         (t,ipa) &
                                                    + csite%fmean_rshortup           (ipa) &
                                                    * ndaysi
               csite%qmean_rnet             (t,ipa) = csite%qmean_rnet             (t,ipa) &
                                                    + csite%fmean_rnet               (ipa) &
                                                    * ndaysi
               csite%qmean_albedo           (t,ipa) = csite%qmean_albedo           (t,ipa) &
                                                    + csite%fmean_albedo             (ipa) &
                                                    * ndaysi
               csite%qmean_albedo_par       (t,ipa) = csite%qmean_albedo_par       (t,ipa) &
                                                    + csite%fmean_albedo_par         (ipa) &
                                                    * ndaysi
               csite%qmean_albedo_nir       (t,ipa) = csite%qmean_albedo_nir       (t,ipa) &
                                                    + csite%fmean_albedo_nir         (ipa) &
                                                    * ndaysi
               csite%qmean_rlong_albedo     (t,ipa) = csite%qmean_rlong_albedo     (t,ipa) &
                                                    + csite%fmean_rlong_albedo       (ipa) &
                                                    * ndaysi
               csite%qmean_ustar            (t,ipa) = csite%qmean_ustar            (t,ipa) &
                                                    + csite%fmean_ustar              (ipa) &
                                                    * ndaysi
               csite%qmean_tstar            (t,ipa) = csite%qmean_tstar            (t,ipa) &
                                                    + csite%fmean_tstar              (ipa) &
                                                    * ndaysi
               csite%qmean_qstar            (t,ipa) = csite%qmean_qstar            (t,ipa) &
                                                    + csite%fmean_qstar              (ipa) &
                                                    * ndaysi
               csite%qmean_cstar            (t,ipa) = csite%qmean_cstar            (t,ipa) &
                                                    + csite%fmean_cstar              (ipa) &
                                                    * ndaysi
               csite%qmean_carbon_ac        (t,ipa) = csite%qmean_carbon_ac        (t,ipa) &
                                                    + csite%fmean_carbon_ac          (ipa) &
                                                    * ndaysi
               csite%qmean_carbon_st        (t,ipa) = csite%qmean_carbon_st        (t,ipa) &
                                                    + csite%fmean_carbon_st          (ipa) &
                                                    * ndaysi
               csite%qmean_vapor_gc         (t,ipa) = csite%qmean_vapor_gc         (t,ipa) &
                                                    + csite%fmean_vapor_gc           (ipa) &
                                                    * ndaysi
               csite%qmean_vapor_ac         (t,ipa) = csite%qmean_vapor_ac         (t,ipa) &
                                                    + csite%fmean_vapor_ac           (ipa) &
                                                    * ndaysi
               csite%qmean_smoist_gg      (:,t,ipa) = csite%qmean_smoist_gg      (:,t,ipa) &
                                                    + csite%fmean_smoist_gg        (:,ipa) &
                                                    * ndaysi
               csite%qmean_throughfall      (t,ipa) = csite%qmean_throughfall      (t,ipa) &
                                                    + csite%fmean_throughfall        (ipa) &
                                                    * ndaysi
               csite%qmean_transloss      (:,t,ipa) = csite%qmean_transloss      (:,t,ipa) &
                                                    + csite%fmean_transloss        (:,ipa) &
                                                    * ndaysi
               csite%qmean_runoff           (t,ipa) = csite%qmean_runoff           (t,ipa) &
                                                    + csite%fmean_runoff             (ipa) &
                                                    * ndaysi
               csite%qmean_drainage         (t,ipa) = csite%qmean_drainage         (t,ipa) &
                                                    + csite%fmean_drainage           (ipa) &
                                                    * ndaysi
               csite%qmean_sensible_gc      (t,ipa) = csite%qmean_sensible_gc      (t,ipa) &
                                                    + csite%fmean_sensible_gc        (ipa) &
                                                    * ndaysi
               csite%qmean_sensible_ac      (t,ipa) = csite%qmean_sensible_ac      (t,ipa) &
                                                    + csite%fmean_sensible_ac        (ipa) &
                                                    * ndaysi
               csite%qmean_sensible_gg    (:,t,ipa) = csite%qmean_sensible_gg    (:,t,ipa) &
                                                    + csite%fmean_sensible_gg      (:,ipa) &
                                                    * ndaysi
               csite%qmean_qthroughfall     (t,ipa) = csite%qmean_qthroughfall     (t,ipa) &
                                                    + csite%fmean_qthroughfall       (ipa) &
                                                    * ndaysi
               csite%qmean_qrunoff          (t,ipa) = csite%qmean_qrunoff          (t,ipa) &
                                                    + csite%fmean_qrunoff            (ipa) &
                                                    * ndaysi
               csite%qmean_qdrainage        (t,ipa) = csite%qmean_qdrainage        (t,ipa) &
                                                    + csite%fmean_qdrainage          (ipa) &
                                                    * ndaysi
               !------ Integrate the mean sum of squares. ---------------------------------!
               csite%qmsqu_rh           (t,ipa) = csite%qmsqu_rh                  (t,ipa)  &
                                                + isqu_ftz(csite%fmean_rh           (ipa)) &
                                                * ndaysi                                 
               csite%qmsqu_cwd_rh       (t,ipa) = csite%qmsqu_cwd_rh              (t,ipa)  &
                                                + isqu_ftz(csite%fmean_cwd_rh       (ipa)) &
                                                * ndaysi
               csite%qmsqu_nep          (t,ipa) = csite%qmsqu_nep                 (t,ipa)  &
                                                + isqu_ftz(csite%fmean_nep          (ipa)) &
                                                * ndaysi
               csite%qmsqu_rlongup      (t,ipa) = csite%qmsqu_rlongup             (t,ipa)  &
                                                + isqu_ftz(csite%fmean_rlongup      (ipa)) &
                                                * ndaysi
               csite%qmsqu_parup        (t,ipa) = csite%qmsqu_parup               (t,ipa)  &
                                                + isqu_ftz(csite%fmean_parup        (ipa)) &
                                                * ndaysi
               csite%qmsqu_nirup        (t,ipa) = csite%qmsqu_nirup               (t,ipa)  &
                                                + isqu_ftz(csite%fmean_nirup        (ipa)) &
                                                * ndaysi                                 
               csite%qmsqu_rshortup     (t,ipa) = csite%qmsqu_rshortup            (t,ipa)  &
                                                + isqu_ftz(csite%fmean_rshortup     (ipa)) &
                                                * ndaysi
               csite%qmsqu_rnet         (t,ipa) = csite%qmsqu_rnet                (t,ipa)  &
                                                + isqu_ftz(csite%fmean_rnet         (ipa)) &
                                                * ndaysi
               csite%qmsqu_albedo       (t,ipa) = csite%qmsqu_albedo              (t,ipa)  &
                                                + isqu_ftz(csite%fmean_albedo       (ipa)) &
                                                * ndaysi
               csite%qmsqu_ustar        (t,ipa) = csite%qmsqu_ustar               (t,ipa)  &
                                                + isqu_ftz(csite%fmean_ustar        (ipa)) &
                                                * ndaysi
               csite%qmsqu_carbon_ac    (t,ipa) = csite%qmsqu_carbon_ac           (t,ipa)  &
                                                + isqu_ftz(csite%fmean_carbon_ac    (ipa)) &
                                                * ndaysi
               csite%qmsqu_carbon_st    (t,ipa) = csite%qmsqu_carbon_st           (t,ipa)  &
                                                + isqu_ftz(csite%fmean_carbon_st    (ipa)) &
                                                * ndaysi
               csite%qmsqu_vapor_gc     (t,ipa) = csite%qmsqu_vapor_gc            (t,ipa)  &
                                                + isqu_ftz(csite%fmean_vapor_gc     (ipa)) &
                                                * ndaysi
               csite%qmsqu_vapor_ac     (t,ipa) = csite%qmsqu_vapor_ac            (t,ipa)  &
                                                + isqu_ftz(csite%fmean_vapor_ac     (ipa)) &
                                                * ndaysi
               csite%qmsqu_sensible_gc  (t,ipa) = csite%qmsqu_sensible_gc         (t,ipa)  &
                                                + isqu_ftz(csite%fmean_sensible_gc  (ipa)) &
                                                * ndaysi
               csite%qmsqu_sensible_ac  (t,ipa) = csite%qmsqu_sensible_ac         (t,ipa)  &
                                                + isqu_ftz(csite%fmean_sensible_ac  (ipa)) &
                                                * ndaysi
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Patch loop.                                                          !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts
                  !------------------------------------------------------------------------!
                  !      Integrate cohort-level variables.                                 !
                  !------------------------------------------------------------------------!
                  cpatch%qmean_gpp           (t,ico) = cpatch%qmean_gpp           (t,ico)  &
                                                     + cpatch%fmean_gpp             (ico)  &
                                                     * ndaysi
                  cpatch%qmean_npp           (t,ico) = cpatch%qmean_npp           (t,ico)  &
                                                     + cpatch%fmean_npp             (ico)  &
                                                     * ndaysi
                  cpatch%qmean_leaf_resp     (t,ico) = cpatch%qmean_leaf_resp     (t,ico)  &
                                                     + cpatch%fmean_leaf_resp       (ico)  &
                                                     * ndaysi
                  cpatch%qmean_root_resp     (t,ico) = cpatch%qmean_root_resp     (t,ico)  &
                                                     + cpatch%fmean_root_resp       (ico)  &
                                                     * ndaysi
                  cpatch%qmean_growth_resp   (t,ico) = cpatch%qmean_growth_resp   (t,ico)  &
                                                     + cpatch%fmean_growth_resp     (ico)  &
                                                     * ndaysi
                  cpatch%qmean_storage_resp  (t,ico) = cpatch%qmean_storage_resp  (t,ico)  &
                                                     + cpatch%fmean_storage_resp    (ico)  &
                                                     * ndaysi
                  cpatch%qmean_vleaf_resp    (t,ico) = cpatch%qmean_vleaf_resp    (t,ico)  &
                                                     + cpatch%fmean_vleaf_resp      (ico)  &
                                                     * ndaysi
                  cpatch%qmean_plresp        (t,ico) = cpatch%qmean_plresp        (t,ico)  &
                                                     + cpatch%fmean_plresp          (ico)  &
                                                     * ndaysi
                  cpatch%qmean_leaf_energy   (t,ico) = cpatch%qmean_leaf_energy   (t,ico)  &
                                                     + cpatch%fmean_leaf_energy     (ico)  &
                                                     * ndaysi
                  cpatch%qmean_leaf_water    (t,ico) = cpatch%qmean_leaf_water    (t,ico)  &
                                                     + cpatch%fmean_leaf_water      (ico)  &
                                                     * ndaysi
                  cpatch%qmean_leaf_hcap     (t,ico) = cpatch%qmean_leaf_hcap     (t,ico)  &
                                                     + cpatch%fmean_leaf_hcap       (ico)  &
                                                     * ndaysi
                  cpatch%qmean_leaf_vpdef    (t,ico) = cpatch%qmean_leaf_vpdef    (t,ico)  &
                                                     + cpatch%fmean_leaf_vpdef      (ico)  &
                                                     * ndaysi
                  cpatch%qmean_leaf_gsw      (t,ico) = cpatch%qmean_leaf_gsw      (t,ico)  &
                                                     + cpatch%fmean_leaf_gsw        (ico)  &
                                                     * ndaysi
                  cpatch%qmean_leaf_gbw      (t,ico) = cpatch%qmean_leaf_gbw      (t,ico)  &
                                                     + cpatch%fmean_leaf_gbw        (ico)  &
                                                     * ndaysi
                  cpatch%qmean_wood_energy   (t,ico) = cpatch%qmean_wood_energy   (t,ico)  &
                                                     + cpatch%fmean_wood_energy     (ico)  &
                                                     * ndaysi
                  cpatch%qmean_wood_water    (t,ico) = cpatch%qmean_wood_water    (t,ico)  &
                                                     + cpatch%fmean_wood_water      (ico)  &
                                                     * ndaysi
                  cpatch%qmean_wood_hcap     (t,ico) = cpatch%qmean_wood_hcap     (t,ico)  &
                                                     + cpatch%fmean_wood_hcap       (ico)  &
                                                     * ndaysi
                  cpatch%qmean_wood_gbw      (t,ico) = cpatch%qmean_wood_gbw      (t,ico)  &
                                                     + cpatch%fmean_wood_gbw        (ico)  &
                                                     * ndaysi
                  cpatch%qmean_fs_open       (t,ico) = cpatch%qmean_fs_open       (t,ico)  &
                                                     + cpatch%fmean_fs_open         (ico)  &
                                                     * ndaysi
                  cpatch%qmean_fsw           (t,ico) = cpatch%qmean_fsw           (t,ico)  &
                                                     + cpatch%fmean_fsw             (ico)  &
                                                     * ndaysi
                  cpatch%qmean_fsn           (t,ico) = cpatch%qmean_fsn           (t,ico)  &
                                                     + cpatch%fmean_fsn             (ico)  &
                                                     * ndaysi
                  cpatch%qmean_a_light       (t,ico) = cpatch%qmean_a_light       (t,ico)  &
                                                     + cpatch%fmean_a_light         (ico)  &
                                                     * ndaysi
                  cpatch%qmean_a_rubp        (t,ico) = cpatch%qmean_a_rubp        (t,ico)  &
                                                     + cpatch%fmean_a_rubp          (ico)  &
                                                     * ndaysi
                  cpatch%qmean_a_co2         (t,ico) = cpatch%qmean_a_co2         (t,ico)  &
                                                     + cpatch%fmean_a_co2           (ico)  &
                                                     * ndaysi
                  cpatch%qmean_psi_open      (t,ico) = cpatch%qmean_psi_open      (t,ico)  &
                                                     + cpatch%fmean_psi_open        (ico)  &
                                                     * ndaysi
                  cpatch%qmean_psi_closed    (t,ico) = cpatch%qmean_psi_closed    (t,ico)  &
                                                     + cpatch%fmean_psi_closed      (ico)  &
                                                     * ndaysi
                  cpatch%qmean_water_supply  (t,ico) = cpatch%qmean_water_supply  (t,ico)  &
                                                     + cpatch%fmean_water_supply    (ico)  &
                                                     * ndaysi
                  
                  cpatch%qmean_par_level_beam(t,ico) = cpatch%qmean_par_level_beam(t,ico) &
                                                     + cpatch%fmean_par_level_beam(ico)    &
                                                     * ndaysi

                  cpatch%qmean_par_level_diffd(t,ico)= cpatch%qmean_par_level_diffd(t,ico) &
                                                     + cpatch%fmean_par_level_diffd(ico)    &
                                                     * ndaysi

                  cpatch%qmean_par_level_diffu(t,ico)= cpatch%qmean_par_level_diffu(t,ico) &
                                                     + cpatch%fmean_par_level_diffu(ico)    &
                                                     * ndaysi

                  cpatch%qmean_light_level   (t,ico) = cpatch%qmean_light_level   (t,ico)  &
                                                     + cpatch%fmean_light_level     (ico)  &
                                                     * ndaysi

                  cpatch%qmean_light_level_beam(t,ico) =                                   &
                                                      cpatch%qmean_light_level_beam(t,ico) &
                                                    + cpatch%fmean_light_level_beam  (ico) &
                                                    * ndaysi
                  cpatch%qmean_light_level_diff(t,ico) =                                   &
                                                      cpatch%qmean_light_level_diff(t,ico) &
                                                    + cpatch%fmean_light_level_diff  (ico) &
                                                    * ndaysi
                  cpatch%qmean_par_l         (t,ico) = cpatch%qmean_par_l         (t,ico)  &
                                                     + cpatch%fmean_par_l           (ico)  &
                                                     * ndaysi
                  cpatch%qmean_par_l_beam    (t,ico) = cpatch%qmean_par_l_beam    (t,ico)  &
                                                     + cpatch%fmean_par_l_beam      (ico)  &
                                                     * ndaysi
                  cpatch%qmean_par_l_diff    (t,ico) = cpatch%qmean_par_l_diff    (t,ico)  &
                                                     + cpatch%fmean_par_l_diff      (ico)  &
                                                     * ndaysi
                  cpatch%qmean_rshort_l      (t,ico) = cpatch%qmean_rshort_l      (t,ico)  &
                                                     + cpatch%fmean_rshort_l        (ico)  &
                                                     * ndaysi
                  cpatch%qmean_rlong_l       (t,ico) = cpatch%qmean_rlong_l       (t,ico)  &
                                                     + cpatch%fmean_rlong_l         (ico)  &
                                                     * ndaysi
                  cpatch%qmean_sensible_lc   (t,ico) = cpatch%qmean_sensible_lc   (t,ico)  &
                                                     + cpatch%fmean_sensible_lc     (ico)  &
                                                     * ndaysi
                  cpatch%qmean_vapor_lc      (t,ico) = cpatch%qmean_vapor_lc      (t,ico)  &
                                                     + cpatch%fmean_vapor_lc        (ico)  &
                                                     * ndaysi
                  cpatch%qmean_transp        (t,ico) = cpatch%qmean_transp        (t,ico)  &
                                                     + cpatch%fmean_transp          (ico)  &
                                                     * ndaysi
                  cpatch%qmean_intercepted_al(t,ico) = cpatch%qmean_intercepted_al(t,ico)  &
                                                     + cpatch%fmean_intercepted_al  (ico)  &
                                                     * ndaysi
                  cpatch%qmean_wshed_lg      (t,ico) = cpatch%qmean_wshed_lg      (t,ico)  &
                                                     + cpatch%fmean_wshed_lg        (ico)  &
                                                     * ndaysi
                  cpatch%qmean_rshort_w      (t,ico) = cpatch%qmean_rshort_w      (t,ico)  &
                                                     + cpatch%fmean_rshort_w        (ico)  &
                                                     * ndaysi
                  cpatch%qmean_rlong_w       (t,ico) = cpatch%qmean_rlong_w       (t,ico)  &
                                                     + cpatch%fmean_rlong_w         (ico)  &
                                                     * ndaysi
                  cpatch%qmean_rad_profile (:,t,ico) = cpatch%qmean_rad_profile (:,t,ico)  &
                                                     + cpatch%fmean_rad_profile   (:,ico)  &
                                                     * ndaysi
                  cpatch%qmean_sensible_wc   (t,ico) = cpatch%qmean_sensible_wc   (t,ico)  &
                                                     + cpatch%fmean_sensible_wc     (ico)  &
                                                     * ndaysi
                  cpatch%qmean_vapor_wc      (t,ico) = cpatch%qmean_vapor_wc      (t,ico)  &
                                                     + cpatch%fmean_vapor_wc        (ico)  &
                                                     * ndaysi
                  cpatch%qmean_intercepted_aw(t,ico) = cpatch%qmean_intercepted_aw(t,ico)  &
                                                     + cpatch%fmean_intercepted_aw  (ico)  &
                                                     * ndaysi
                  cpatch%qmean_wshed_wg      (t,ico) = cpatch%qmean_wshed_wg      (t,ico)  &
                                                     + cpatch%fmean_wshed_wg        (ico)  &
                                                     * ndaysi
                  !------ Mean sum of squares. --------------------------------------------!
                  cpatch%qmsqu_gpp        (t,ico) = cpatch%qmsqu_gpp              (t,ico)  &
                                                  + isqu_ftz(cpatch%fmean_gpp       (ico)) &
                                                  * ndaysi
                  cpatch%qmsqu_npp        (t,ico) = cpatch%qmsqu_npp              (t,ico)  &
                                                  + isqu_ftz(cpatch%fmean_npp       (ico)) &
                                                  * ndaysi
                  cpatch%qmsqu_plresp     (t,ico) = cpatch%qmsqu_plresp           (t,ico)  &
                                                  + isqu_ftz(cpatch%fmean_plresp    (ico)) &
                                                  * ndaysi
                  cpatch%qmsqu_sensible_lc(t,ico) = cpatch%qmsqu_sensible_lc      (t,ico)  &
                                                + isqu_ftz(cpatch%fmean_sensible_lc (ico)) &
                                                * ndaysi
                  cpatch%qmsqu_vapor_lc   (t,ico) = cpatch%qmsqu_vapor_lc         (t,ico)  &
                                                  + isqu_ftz(cpatch%fmean_vapor_lc  (ico)) &
                                                  * ndaysi
                  cpatch%qmsqu_transp     (t,ico) = cpatch%qmsqu_transp           (t,ico)  &
                                                  + isqu_ftz(cpatch%fmean_transp    (ico)) &
                                                  * ndaysi
                  cpatch%qmsqu_sensible_wc(t,ico) = cpatch%qmsqu_sensible_wc      (t,ico)  &
                                                + isqu_ftz(cpatch%fmean_sensible_wc (ico)) &
                                                * ndaysi
                  cpatch%qmsqu_vapor_wc   (t,ico) = cpatch%qmsqu_vapor_wc         (t,ico)  &
                                                  + isqu_ftz(cpatch%fmean_vapor_wc  (ico)) &
                                                  * ndaysi
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!
      return
   end subroutine integrate_ed_qmean_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine normalises the daily mean variables of those variables that could  !
   ! not be integrated directly.  These are pretty much just temperature, density, and     !
   ! liquid water fraction.                                                                !
   !---------------------------------------------------------------------------------------!
   subroutine normalize_ed_qmean_vars(cgrid)
      use ed_state_vars        , only : edtype             & ! structure
                                      , polygontype        & ! structure
                                      , sitetype           & ! structure
                                      , patchtype          ! ! structure
      use ed_misc_coms         , only : ndcycle            ! ! intent(in)
      use therm_lib            , only : press2exner        & ! function
                                      , extheta2temp       & ! function
                                      , uextcm2tl          & ! subroutine
                                      , uint2tl            & ! subroutine
                                      , idealdenssh        ! ! function
      use soil_coms            , only : tiny_sfcwater_mass & ! intent(in)
                                      , isoilbc            & ! intent(in)
                                      , soil               & ! intent(in)
                                      , dslz               ! ! intent(in)
      use consts_coms          , only : t00                & ! intent(in)
                                      , wdns               ! ! intent(in)
      use grid_coms            , only : nzg                ! ! intent(in)
      implicit none
      !----- Argument ---------------------------------------------------------------------!
      type(edtype)                     , target  :: cgrid
      !----- Local variables --------------------------------------------------------------!
      type(polygontype)                , pointer :: cpoly
      type(sitetype)                   , pointer :: csite
      type(patchtype)                  , pointer :: cpatch
      real             , dimension(nzg)          :: cgrid_qmean_soil_hcap
      integer                                    :: ipy
      integer                                    :: isi
      integer                                    :: ipa
      integer                                    :: ico
      integer                                    :: nsoil
      integer                                    :: k
      integer                                    :: t
      real                                       :: can_exner
      real                                       :: atm_exner
      real                                       :: site_area_i
      real                                       :: poly_area_i
      real                                       :: site_wgt
      real                                       :: patch_wgt
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      Use the variables that have been already aggregated.                          !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)


         !----- Inverse of this polygon area (it should be always 1.) ---------------------!
         poly_area_i = 1./sum(cpoly%area)
         !---------------------------------------------------------------------------------!

         !----- Re-set some support variables. --------------------------------------------!
         cgrid_qmean_soil_hcap(:) = 0.0
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Site loop.                                                                 !
         !---------------------------------------------------------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            !----- Inverse of this site area (it should be always 1.) ---------------------!
            site_area_i = 1./sum(csite%area)
            !------------------------------------------------------------------------------!


            !----- Site weight. -----------------------------------------------------------!
            site_wgt = cpoly%area(isi) * poly_area_i
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Find the derived properties for the air above canopy.                   !
            !------------------------------------------------------------------------------!
            do t=1,ndcycle
               atm_exner                   = press2exner ( cpoly%qmean_atm_prss  (t,isi) )
               cpoly%qmean_atm_temp(t,isi) = extheta2temp( atm_exner                       &
                                                         , cpoly%qmean_atm_theta (t,isi) )
               cpoly%qmean_atm_rhos(t,isi) = idealdenssh ( cpoly%qmean_atm_prss  (t,isi)   &
                                                         , cpoly%qmean_atm_temp  (t,isi)   &
                                                         , cpoly%qmean_atm_shv   (t,isi) )
            end do
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Patch loop.                                                             !
            !------------------------------------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)


               !----- Site weight. --------------------------------------------------------!
               patch_wgt = csite%area(ipa) * site_area_i * site_wgt
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !     Soil matric potential, temperature, and liquid water.                 !
               !---------------------------------------------------------------------------!
               do k=1,nzg
                  nsoil = cpoly%ntext_soil(k,isi)

                  !----- Heat capacity stays outside the time loop. -----------------------!
                  cgrid_qmean_soil_hcap(k) = cgrid_qmean_soil_hcap(k)                      &
                                           + soil(nsoil)%slcpd * patch_wgt
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Find the soil properties, and integrate polygon-level soil matric  !
                  ! potential.                                                             !
                  !------------------------------------------------------------------------!
                  do t=1,ndcycle
                     call uextcm2tl( csite%qmean_soil_energy(k,t,ipa)                      &
                                   , csite%qmean_soil_water (k,t,ipa) * wdns               &
                                   , soil(nsoil)%slcpd                                     &
                                   , csite%qmean_soil_temp  (k,t,ipa)                      &
                                   , csite%qmean_soil_fliq  (k,t,ipa))
                  end do
                  !------------------------------------------------------------------------!
               end do
               !---------------------------------------------------------------------------!




               !---------------------------------------------------------------------------!
               !   If the patch had some temporary snow/pounding layer, convert the mean   !
               ! energy to J/kg, then find the mean temperature and liquid fraction.       !
               ! Otherwise, set them to either zero or default values.                     !
               !---------------------------------------------------------------------------!
               do t=1,ndcycle
                  if (csite%qmean_sfcw_mass(t,ipa) > tiny_sfcwater_mass) then
                     csite%qmean_sfcw_energy(t,ipa) = csite%qmean_sfcw_energy(t,ipa)       &
                                                    / csite%qmean_sfcw_mass  (t,ipa)
                     call uint2tl( csite%qmean_sfcw_energy(t,ipa)                          &
                                 , csite%qmean_sfcw_temp  (t,ipa)                          &
                                 , csite%qmean_sfcw_fliq  (t,ipa))
                  else
                     csite%qmean_sfcw_mass  (t,ipa)  = 0.
                     csite%qmean_sfcw_depth (t,ipa)  = 0.
                     csite%qmean_sfcw_energy(t,ipa)  = 0.
                     csite%qmean_sfcw_temp  (t,ipa)  = csite%qmean_soil_temp(nzg,t,ipa)
                     csite%qmean_sfcw_fliq  (t,ipa)  = csite%qmean_soil_fliq(nzg,t,ipa)
                  end if
               end do
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      Cohort loop.                                                         !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1,cpatch%ncohorts



                  !------------------------------------------------------------------------!
                  !     Find the vegetation temperature and liquid fraction.               !
                  !------------------------------------------------------------------------!
                  do t=1,ndcycle
                     !----- Leaf. ---------------------------------------------------------!
                     if (cpatch%qmean_leaf_hcap(t,ico) > 0.) then
                        call uextcm2tl( cpatch%qmean_leaf_energy(t,ico)                    &
                                      , cpatch%qmean_leaf_water (t,ico)                    &
                                      , cpatch%qmean_leaf_hcap  (t,ico)                    &
                                      , cpatch%qmean_leaf_temp  (t,ico)                    &
                                      , cpatch%qmean_leaf_fliq  (t,ico) )
                     else
                        cpatch%qmean_leaf_vpdef(t,ico) = csite%qmean_can_vpdef(t,ipa)
                        cpatch%qmean_leaf_temp (t,ico) = csite%qmean_can_temp (t,ipa)
                        if (csite%qmean_can_temp(t,ipa) > t00) then
                           cpatch%qmean_leaf_fliq(t,ico) = 1.0
                        elseif (csite%qmean_can_temp(t,ipa) == t00) then
                           cpatch%qmean_leaf_fliq(t,ico) = 0.5
                        else
                           cpatch%qmean_leaf_fliq(t,ico) = 0.0
                        end if
                     end if
                     !----- Wood. ---------------------------------------------------------!
                     if (cpatch%qmean_wood_hcap(t,ico) > 0.) then
                        call uextcm2tl( cpatch%qmean_wood_energy(t,ico)                    &
                                      , cpatch%qmean_wood_water (t,ico)                    &
                                      , cpatch%qmean_wood_hcap  (t,ico)                    &
                                      , cpatch%qmean_wood_temp  (t,ico)                    &
                                      , cpatch%qmean_wood_fliq  (t,ico) )
                     else
                        cpatch%qmean_wood_temp(t,ico) = csite%qmean_can_temp(t,ipa)
                        if (csite%qmean_can_temp(t,ipa) > t00) then
                           cpatch%qmean_wood_fliq(t,ico) = 1.0
                        elseif (csite%qmean_can_temp(t,ipa) == t00) then
                           cpatch%qmean_wood_fliq(t,ico) = 0.5
                        else
                           cpatch%qmean_wood_fliq(t,ico) = 0.0
                        end if
                        !------------------------------------------------------------------!
                     end if
                     !---------------------------------------------------------------------!
                  end do
                  !------------------------------------------------------------------------!
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!






         !---------------------------------------------------------------------------------!
         !      Find the derived properties for the air above canopy.                      !
         !---------------------------------------------------------------------------------!
         do t=1,ndcycle
            !------------------------------------------------------------------------------!
            !      Find the derived properties for the air above canopy.                   !
            !------------------------------------------------------------------------------!
            atm_exner                   = press2exner ( cgrid%qmean_atm_prss     (t,ipy) )
            cgrid%qmean_atm_temp(t,ipy) = extheta2temp( atm_exner                          &
                                                      , cgrid%qmean_atm_theta    (t,ipy) )
            cgrid%qmean_atm_rhos(t,ipy) = idealdenssh ( cgrid%qmean_atm_prss     (t,ipy)   &
                                                      , cgrid%qmean_atm_temp     (t,ipy)   &
                                                      , cgrid%qmean_atm_shv      (t,ipy) )
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !      Find the derived properties for the canopy air space.                   !
            !------------------------------------------------------------------------------!
            can_exner                   = press2exner ( cgrid%qmean_can_prss (t,ipy) )
            cgrid%qmean_can_temp(t,ipy) = extheta2temp( can_exner                          &
                                                      , cgrid%qmean_can_theta(t,ipy) )
            cgrid%qmean_can_rhos(t,ipy) = idealdenssh ( cgrid%qmean_can_prss (t,ipy)       &
                                                      , cgrid%qmean_can_temp (t,ipy)       &
                                                      , cgrid%qmean_can_shv  (t,ipy) )
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !   If the patch had some temporary snow/pounding layer, convert the mean      !
            ! energy to J/kg, then find the mean temperature and liquid fraction.  Other-  !
            ! wise, set them to either zero or default values.                             !
            !------------------------------------------------------------------------------!
            if (cgrid%qmean_sfcw_mass(t,ipy) > tiny_sfcwater_mass) then
               cgrid%qmean_sfcw_energy(t,ipy) = cgrid%qmean_sfcw_energy(t,ipy)             &
                                              / cgrid%qmean_sfcw_mass  (t,ipy)
               call uint2tl(cgrid%qmean_sfcw_energy(t,ipy),cgrid%qmean_sfcw_temp(t,ipy)    &
                           ,cgrid%qmean_sfcw_fliq(t,ipy))
            else
               cgrid%qmean_sfcw_mass  (t,ipy)  = 0.
               cgrid%qmean_sfcw_depth (t,ipy)  = 0.
               cgrid%qmean_sfcw_energy(t,ipy)  = 0.
               cgrid%qmean_sfcw_temp  (t,ipy)  = cgrid%qmean_soil_temp(nzg,t,ipy)
               cgrid%qmean_sfcw_fliq  (t,ipy)  = cgrid%qmean_soil_fliq(nzg,t,ipy)
            end if
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Find the temperature and the fraction of liquid water.                   !
            !------------------------------------------------------------------------------!
            do k=1,nzg
               call uextcm2tl( cgrid%qmean_soil_energy(k,t,ipy)                            &
                             , cgrid%qmean_soil_water (k,t,ipy) * wdns                     &
                             , cgrid_qmean_soil_hcap  (k)                                  &
                             , cgrid%qmean_soil_temp  (k,t,ipy)                            &
                             , cgrid%qmean_soil_fliq  (k,t,ipy) )
            end do
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Find the vegetation temperature and liquid fraction.                     !
            !------------------------------------------------------------------------------!
            !----- Leaf. ------------------------------------------------------------------!
            if (cgrid%qmean_leaf_hcap(t,ipy) > 0.) then
               call uextcm2tl( cgrid%qmean_leaf_energy(t,ipy)                              &
                             , cgrid%qmean_leaf_water (t,ipy)                              &
                             , cgrid%qmean_leaf_hcap  (t,ipy)                              &
                             , cgrid%qmean_leaf_temp  (t,ipy)                              &
                             , cgrid%qmean_leaf_fliq  (t,ipy) )
            else
               cgrid%qmean_leaf_temp (t,ipy) = cgrid%qmean_can_temp (t,ipy)
               if (cgrid%qmean_can_temp(t,ipy) > t00) then
                  cgrid%qmean_leaf_fliq(t,ipy) = 1.0
               elseif (cgrid%qmean_can_temp(t,ipy) == t00) then
                  cgrid%qmean_leaf_fliq(t,ipy) = 0.5
               else
                  cgrid%qmean_leaf_fliq(t,ipy) = 0.0
               end if
            end if
            !----- Wood. ------------------------------------------------------------------!
            if (cgrid%qmean_wood_hcap(t,ipy) > 0.) then
               call uextcm2tl( cgrid%qmean_wood_energy(t,ipy)                              &
                             , cgrid%qmean_wood_water (t,ipy)                              &
                             , cgrid%qmean_wood_hcap  (t,ipy)                              &
                             , cgrid%qmean_wood_temp  (t,ipy)                              &
                             , cgrid%qmean_wood_fliq  (t,ipy) )
            else
               cgrid%qmean_wood_temp(t,ipy) = cgrid%qmean_can_temp(t,ipy)
               if (cgrid%qmean_can_temp(t,ipy) > t00) then
                  cgrid%qmean_wood_fliq(t,ipy) = 1.0
               elseif (cgrid%qmean_can_temp(t,ipy) == t00) then
                  cgrid%qmean_wood_fliq(t,ipy) = 0.5
               else
                  cgrid%qmean_wood_fliq(t,ipy) = 0.0
               end if
            end if
         end do
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!
      return
   end subroutine normalize_ed_qmean_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine resets the mean diel once the "Q" file has been written.           !
   !---------------------------------------------------------------------------------------!
   subroutine zero_ed_qmean_vars(cgrid)
      use ed_state_vars, only : edtype        & ! structure
                              , polygontype   & ! structure
                              , sitetype      & ! structure
                              , patchtype     ! ! structure
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(edtype)     , target  :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer :: cpoly
      type(sitetype)   , pointer :: csite
      type(patchtype)  , pointer :: cpatch
      integer                    :: ipy
      integer                    :: isi
      integer                    :: ipa
      integer                    :: ico
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !       Loop over polygons.                                                          !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         cgrid%qmean_gpp                (:,ipy) = 0.0
         cgrid%qmean_npp                (:,ipy) = 0.0
         cgrid%qmean_leaf_resp          (:,ipy) = 0.0
         cgrid%qmean_root_resp          (:,ipy) = 0.0
         cgrid%qmean_growth_resp        (:,ipy) = 0.0
         cgrid%qmean_storage_resp       (:,ipy) = 0.0
         cgrid%qmean_vleaf_resp         (:,ipy) = 0.0
         cgrid%qmean_plresp             (:,ipy) = 0.0
         cgrid%qmean_leaf_energy        (:,ipy) = 0.0
         cgrid%qmean_leaf_water         (:,ipy) = 0.0
         cgrid%qmean_leaf_hcap          (:,ipy) = 0.0
         cgrid%qmean_leaf_vpdef         (:,ipy) = 0.0
         cgrid%qmean_leaf_temp          (:,ipy) = 0.0
         cgrid%qmean_leaf_fliq          (:,ipy) = 0.0
         cgrid%qmean_leaf_gsw           (:,ipy) = 0.0
         cgrid%qmean_leaf_gbw           (:,ipy) = 0.0
         cgrid%qmean_wood_energy        (:,ipy) = 0.0
         cgrid%qmean_wood_water         (:,ipy) = 0.0
         cgrid%qmean_wood_hcap          (:,ipy) = 0.0
         cgrid%qmean_wood_temp          (:,ipy) = 0.0
         cgrid%qmean_wood_fliq          (:,ipy) = 0.0
         cgrid%qmean_wood_gbw           (:,ipy) = 0.0
         cgrid%qmean_fs_open            (:,ipy) = 0.0
         cgrid%qmean_fsw                (:,ipy) = 0.0
         cgrid%qmean_fsn                (:,ipy) = 0.0
         cgrid%qmean_a_light            (:,ipy) = 0.0
         cgrid%qmean_a_rubp             (:,ipy) = 0.0
         cgrid%qmean_a_co2              (:,ipy) = 0.0
         cgrid%qmean_psi_open           (:,ipy) = 0.0
         cgrid%qmean_psi_closed         (:,ipy) = 0.0
         cgrid%qmean_water_supply       (:,ipy) = 0.0
         cgrid%qmean_par_l              (:,ipy) = 0.0
         cgrid%qmean_par_l_beam         (:,ipy) = 0.0
         cgrid%qmean_par_l_diff         (:,ipy) = 0.0
         cgrid%qmean_rshort_l           (:,ipy) = 0.0
         cgrid%qmean_rlong_l            (:,ipy) = 0.0
         cgrid%qmean_sensible_lc        (:,ipy) = 0.0
         cgrid%qmean_vapor_lc           (:,ipy) = 0.0
         cgrid%qmean_transp             (:,ipy) = 0.0
         cgrid%qmean_intercepted_al     (:,ipy) = 0.0
         cgrid%qmean_wshed_lg           (:,ipy) = 0.0
         cgrid%qmean_rshort_w           (:,ipy) = 0.0
         cgrid%qmean_rlong_w            (:,ipy) = 0.0
         cgrid%qmean_sensible_wc        (:,ipy) = 0.0
         cgrid%qmean_vapor_wc           (:,ipy) = 0.0
         cgrid%qmean_intercepted_aw     (:,ipy) = 0.0
         cgrid%qmean_wshed_wg           (:,ipy) = 0.0
         cgrid%qmean_rh                 (:,ipy) = 0.0
         cgrid%qmean_cwd_rh             (:,ipy) = 0.0
         cgrid%qmean_nep                (:,ipy) = 0.0
         cgrid%qmean_rk4step            (:,ipy) = 0.0
         cgrid%qmean_available_water    (:,ipy) = 0.0
         cgrid%qmean_can_theiv          (:,ipy) = 0.0
         cgrid%qmean_can_theta          (:,ipy) = 0.0
         cgrid%qmean_can_vpdef          (:,ipy) = 0.0
         cgrid%qmean_can_temp           (:,ipy) = 0.0
         cgrid%qmean_can_shv            (:,ipy) = 0.0
         cgrid%qmean_can_co2            (:,ipy) = 0.0
         cgrid%qmean_can_rhos           (:,ipy) = 0.0
         cgrid%qmean_can_prss           (:,ipy) = 0.0
         cgrid%qmean_gnd_temp           (:,ipy) = 0.0
         cgrid%qmean_gnd_shv            (:,ipy) = 0.0
         cgrid%qmean_can_ggnd           (:,ipy) = 0.0
         cgrid%qmean_sfcw_depth         (:,ipy) = 0.0
         cgrid%qmean_sfcw_energy        (:,ipy) = 0.0
         cgrid%qmean_sfcw_mass          (:,ipy) = 0.0
         cgrid%qmean_sfcw_temp          (:,ipy) = 0.0
         cgrid%qmean_sfcw_fliq          (:,ipy) = 0.0
         cgrid%qmean_soil_energy      (:,:,ipy) = 0.0
         cgrid%qmean_soil_mstpot      (:,:,ipy) = 0.0
         cgrid%qmean_soil_water       (:,:,ipy) = 0.0
         cgrid%qmean_soil_temp        (:,:,ipy) = 0.0
         cgrid%qmean_soil_fliq        (:,:,ipy) = 0.0
         cgrid%qmean_rshort_gnd         (:,ipy) = 0.0
         cgrid%qmean_par_gnd            (:,ipy) = 0.0
         cgrid%qmean_rlong_gnd          (:,ipy) = 0.0
         cgrid%qmean_rlongup            (:,ipy) = 0.0
         cgrid%qmean_parup              (:,ipy) = 0.0
         cgrid%qmean_nirup              (:,ipy) = 0.0
         cgrid%qmean_rshortup           (:,ipy) = 0.0
         cgrid%qmean_rnet               (:,ipy) = 0.0
         cgrid%qmean_albedo             (:,ipy) = 0.0
         cgrid%qmean_albedo_par         (:,ipy) = 0.0
         cgrid%qmean_albedo_nir         (:,ipy) = 0.0
         cgrid%qmean_rlong_albedo       (:,ipy) = 0.0
         cgrid%qmean_ustar              (:,ipy) = 0.0
         cgrid%qmean_tstar              (:,ipy) = 0.0
         cgrid%qmean_qstar              (:,ipy) = 0.0
         cgrid%qmean_cstar              (:,ipy) = 0.0
         cgrid%qmean_carbon_ac          (:,ipy) = 0.0
         cgrid%qmean_carbon_st          (:,ipy) = 0.0
         cgrid%qmean_vapor_gc           (:,ipy) = 0.0
         cgrid%qmean_vapor_ac           (:,ipy) = 0.0
         cgrid%qmean_smoist_gg        (:,:,ipy) = 0.0
         cgrid%qmean_throughfall        (:,ipy) = 0.0
         cgrid%qmean_transloss        (:,:,ipy) = 0.0
         cgrid%qmean_runoff             (:,ipy) = 0.0
         cgrid%qmean_drainage           (:,ipy) = 0.0
         cgrid%qmean_sensible_gc        (:,ipy) = 0.0
         cgrid%qmean_sensible_ac        (:,ipy) = 0.0
         cgrid%qmean_sensible_gg      (:,:,ipy) = 0.0
         cgrid%qmean_qthroughfall       (:,ipy) = 0.0
         cgrid%qmean_qrunoff            (:,ipy) = 0.0
         cgrid%qmean_qdrainage          (:,ipy) = 0.0
         cgrid%qmean_atm_theiv          (:,ipy) = 0.0
         cgrid%qmean_atm_theta          (:,ipy) = 0.0
         cgrid%qmean_atm_temp           (:,ipy) = 0.0
         cgrid%qmean_atm_vpdef          (:,ipy) = 0.0
         cgrid%qmean_atm_shv            (:,ipy) = 0.0
         cgrid%qmean_atm_rshort         (:,ipy) = 0.0
         cgrid%qmean_atm_rshort_diff    (:,ipy) = 0.0
         cgrid%qmean_atm_par            (:,ipy) = 0.0
         cgrid%qmean_atm_par_diff       (:,ipy) = 0.0
         cgrid%qmean_atm_rlong          (:,ipy) = 0.0
         cgrid%qmean_atm_vels           (:,ipy) = 0.0
         cgrid%qmean_atm_rhos           (:,ipy) = 0.0
         cgrid%qmean_atm_prss           (:,ipy) = 0.0
         cgrid%qmean_atm_co2            (:,ipy) = 0.0
         cgrid%qmean_pcpg               (:,ipy) = 0.0
         cgrid%qmean_qpcpg              (:,ipy) = 0.0
         cgrid%qmean_dpcpg              (:,ipy) = 0.0
         cgrid%qmsqu_gpp                (:,ipy) = 0.0
         cgrid%qmsqu_npp                (:,ipy) = 0.0
         cgrid%qmsqu_plresp             (:,ipy) = 0.0
         cgrid%qmsqu_sensible_lc        (:,ipy) = 0.0
         cgrid%qmsqu_vapor_lc           (:,ipy) = 0.0
         cgrid%qmsqu_transp             (:,ipy) = 0.0
         cgrid%qmsqu_sensible_wc        (:,ipy) = 0.0
         cgrid%qmsqu_vapor_wc           (:,ipy) = 0.0
         cgrid%qmsqu_rh                 (:,ipy) = 0.0
         cgrid%qmsqu_cwd_rh             (:,ipy) = 0.0
         cgrid%qmsqu_nep                (:,ipy) = 0.0
         cgrid%qmsqu_rlongup            (:,ipy) = 0.0
         cgrid%qmsqu_parup              (:,ipy) = 0.0
         cgrid%qmsqu_nirup              (:,ipy) = 0.0
         cgrid%qmsqu_rshortup           (:,ipy) = 0.0
         cgrid%qmsqu_rnet               (:,ipy) = 0.0
         cgrid%qmsqu_albedo             (:,ipy) = 0.0
         cgrid%qmsqu_ustar              (:,ipy) = 0.0
         cgrid%qmsqu_carbon_ac          (:,ipy) = 0.0
         cgrid%qmsqu_carbon_st          (:,ipy) = 0.0
         cgrid%qmsqu_vapor_gc           (:,ipy) = 0.0
         cgrid%qmsqu_vapor_ac           (:,ipy) = 0.0
         cgrid%qmsqu_sensible_gc        (:,ipy) = 0.0
         cgrid%qmsqu_sensible_ac        (:,ipy) = 0.0


         !---------------------------------------------------------------------------------!
         !       Loop over sites.                                                          !
         !---------------------------------------------------------------------------------!
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            cpoly%qmean_atm_theiv      (:,isi) = 0.0
            cpoly%qmean_atm_theta      (:,isi) = 0.0
            cpoly%qmean_atm_temp       (:,isi) = 0.0
            cpoly%qmean_atm_vpdef      (:,isi) = 0.0
            cpoly%qmean_atm_shv        (:,isi) = 0.0
            cpoly%qmean_atm_rshort     (:,isi) = 0.0
            cpoly%qmean_atm_rshort_diff(:,isi) = 0.0
            cpoly%qmean_atm_par        (:,isi) = 0.0
            cpoly%qmean_atm_par_diff   (:,isi) = 0.0
            cpoly%qmean_atm_rlong      (:,isi) = 0.0
            cpoly%qmean_atm_vels       (:,isi) = 0.0
            cpoly%qmean_atm_rhos       (:,isi) = 0.0
            cpoly%qmean_atm_prss       (:,isi) = 0.0
            cpoly%qmean_atm_co2        (:,isi) = 0.0
            cpoly%qmean_pcpg           (:,isi) = 0.0
            cpoly%qmean_qpcpg          (:,isi) = 0.0
            cpoly%qmean_dpcpg          (:,isi) = 0.0

            !------------------------------------------------------------------------------!
            !       Loop over sites.                                                       !
            !------------------------------------------------------------------------------!
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)
               csite%qmean_rh                     (:,ipa) = 0.0
               csite%qmean_cwd_rh                 (:,ipa) = 0.0
               csite%qmean_nep                    (:,ipa) = 0.0
               csite%qmean_rk4step                (:,ipa) = 0.0
               csite%qmean_available_water        (:,ipa) = 0.0
               csite%qmean_can_theiv              (:,ipa) = 0.0
               csite%qmean_can_theta              (:,ipa) = 0.0
               csite%qmean_can_vpdef              (:,ipa) = 0.0
               csite%qmean_can_temp               (:,ipa) = 0.0
               csite%qmean_can_shv                (:,ipa) = 0.0
               csite%qmean_can_co2                (:,ipa) = 0.0
               csite%qmean_can_rhos               (:,ipa) = 0.0
               csite%qmean_can_prss               (:,ipa) = 0.0
               csite%qmean_gnd_temp               (:,ipa) = 0.0
               csite%qmean_gnd_shv                (:,ipa) = 0.0
               csite%qmean_can_ggnd               (:,ipa) = 0.0
               csite%qmean_sfcw_depth             (:,ipa) = 0.0
               csite%qmean_sfcw_energy            (:,ipa) = 0.0
               csite%qmean_sfcw_mass              (:,ipa) = 0.0
               csite%qmean_sfcw_temp              (:,ipa) = 0.0
               csite%qmean_sfcw_fliq              (:,ipa) = 0.0
               csite%qmean_soil_energy          (:,:,ipa) = 0.0
               csite%qmean_soil_mstpot          (:,:,ipa) = 0.0
               csite%qmean_soil_water           (:,:,ipa) = 0.0
               csite%qmean_soil_temp            (:,:,ipa) = 0.0
               csite%qmean_soil_fliq            (:,:,ipa) = 0.0
               csite%qmean_rshort_gnd             (:,ipa) = 0.0
               csite%qmean_par_gnd                (:,ipa) = 0.0
               csite%qmean_rlong_gnd              (:,ipa) = 0.0
               csite%qmean_rlongup                (:,ipa) = 0.0
               csite%qmean_parup                  (:,ipa) = 0.0
               csite%qmean_nirup                  (:,ipa) = 0.0
               csite%qmean_rshortup               (:,ipa) = 0.0
               csite%qmean_rnet                   (:,ipa) = 0.0
               csite%qmean_albedo                 (:,ipa) = 0.0
               csite%qmean_albedo_par             (:,ipa) = 0.0
               csite%qmean_albedo_nir             (:,ipa) = 0.0
               csite%qmean_rlong_albedo           (:,ipa) = 0.0
               csite%qmean_ustar                  (:,ipa) = 0.0
               csite%qmean_tstar                  (:,ipa) = 0.0
               csite%qmean_qstar                  (:,ipa) = 0.0
               csite%qmean_cstar                  (:,ipa) = 0.0
               csite%qmean_carbon_ac              (:,ipa) = 0.0
               csite%qmean_carbon_st              (:,ipa) = 0.0
               csite%qmean_vapor_gc               (:,ipa) = 0.0
               csite%qmean_vapor_ac               (:,ipa) = 0.0
               csite%qmean_smoist_gg            (:,:,ipa) = 0.0
               csite%qmean_throughfall            (:,ipa) = 0.0
               csite%qmean_transloss            (:,:,ipa) = 0.0
               csite%qmean_runoff                 (:,ipa) = 0.0
               csite%qmean_drainage               (:,ipa) = 0.0
               csite%qmean_sensible_gc            (:,ipa) = 0.0
               csite%qmean_sensible_ac            (:,ipa) = 0.0
               csite%qmean_sensible_gg          (:,:,ipa) = 0.0
               csite%qmean_qthroughfall           (:,ipa) = 0.0
               csite%qmean_qrunoff                (:,ipa) = 0.0
               csite%qmean_qdrainage              (:,ipa) = 0.0
               csite%qmsqu_rh                     (:,ipa) = 0.0
               csite%qmsqu_cwd_rh                 (:,ipa) = 0.0
               csite%qmsqu_nep                    (:,ipa) = 0.0
               csite%qmsqu_rlongup                (:,ipa) = 0.0
               csite%qmsqu_parup                  (:,ipa) = 0.0
               csite%qmsqu_nirup                  (:,ipa) = 0.0
               csite%qmsqu_rshortup               (:,ipa) = 0.0
               csite%qmsqu_rnet                   (:,ipa) = 0.0
               csite%qmsqu_albedo                 (:,ipa) = 0.0
               csite%qmsqu_ustar                  (:,ipa) = 0.0
               csite%qmsqu_carbon_ac              (:,ipa) = 0.0
               csite%qmsqu_carbon_st              (:,ipa) = 0.0
               csite%qmsqu_vapor_gc               (:,ipa) = 0.0
               csite%qmsqu_vapor_ac               (:,ipa) = 0.0
               csite%qmsqu_sensible_gc            (:,ipa) = 0.0
               csite%qmsqu_sensible_ac            (:,ipa) = 0.0



               !---------------------------------------------------------------------------!
               !       Loop over cohorts.                                                  !
               !---------------------------------------------------------------------------!
               cohortloop: do ico=1, cpatch%ncohorts
                  cpatch%qmean_gpp                 (:,ico) = 0.0
                  cpatch%qmean_npp                 (:,ico) = 0.0
                  cpatch%qmean_leaf_resp           (:,ico) = 0.0
                  cpatch%qmean_root_resp           (:,ico) = 0.0
                  cpatch%qmean_growth_resp         (:,ico) = 0.0
                  cpatch%qmean_storage_resp        (:,ico) = 0.0
                  cpatch%qmean_vleaf_resp          (:,ico) = 0.0
                  cpatch%qmean_plresp              (:,ico) = 0.0
                  cpatch%qmean_leaf_energy         (:,ico) = 0.0
                  cpatch%qmean_leaf_water          (:,ico) = 0.0
                  cpatch%qmean_leaf_hcap           (:,ico) = 0.0
                  cpatch%qmean_leaf_vpdef          (:,ico) = 0.0
                  cpatch%qmean_leaf_temp           (:,ico) = 0.0
                  cpatch%qmean_leaf_fliq           (:,ico) = 0.0
                  cpatch%qmean_leaf_gsw            (:,ico) = 0.0
                  cpatch%qmean_leaf_gbw            (:,ico) = 0.0
                  cpatch%qmean_wood_energy         (:,ico) = 0.0
                  cpatch%qmean_wood_water          (:,ico) = 0.0
                  cpatch%qmean_wood_hcap           (:,ico) = 0.0
                  cpatch%qmean_wood_temp           (:,ico) = 0.0
                  cpatch%qmean_wood_fliq           (:,ico) = 0.0
                  cpatch%qmean_wood_gbw            (:,ico) = 0.0
                  cpatch%qmean_fs_open             (:,ico) = 0.0
                  cpatch%qmean_fsw                 (:,ico) = 0.0
                  cpatch%qmean_fsn                 (:,ico) = 0.0
                  cpatch%qmean_a_light             (:,ico) = 0.0
                  cpatch%qmean_a_rubp              (:,ico) = 0.0
                  cpatch%qmean_a_co2               (:,ico) = 0.0
                  cpatch%qmean_psi_open            (:,ico) = 0.0
                  cpatch%qmean_psi_closed          (:,ico) = 0.0
                  cpatch%qmean_water_supply        (:,ico) = 0.0

                  cpatch%qmean_par_level_beam      (:,ico) = 0.0
                  cpatch%qmean_par_level_diffu     (:,ico) = 0.0
                  cpatch%qmean_par_level_diffd     (:,ico) = 0.0

                  cpatch%qmean_light_level         (:,ico) = 0.0
                  cpatch%qmean_light_level_beam    (:,ico) = 0.0
                  cpatch%qmean_light_level_diff    (:,ico) = 0.0
                  cpatch%qmean_par_l               (:,ico) = 0.0
                  cpatch%qmean_par_l_beam          (:,ico) = 0.0
                  cpatch%qmean_par_l_diff          (:,ico) = 0.0
                  cpatch%qmean_rshort_l            (:,ico) = 0.0
                  cpatch%qmean_rlong_l             (:,ico) = 0.0
                  cpatch%qmean_sensible_lc         (:,ico) = 0.0
                  cpatch%qmean_vapor_lc            (:,ico) = 0.0
                  cpatch%qmean_transp              (:,ico) = 0.0
                  cpatch%qmean_intercepted_al      (:,ico) = 0.0
                  cpatch%qmean_wshed_lg            (:,ico) = 0.0
                  cpatch%qmean_rshort_w            (:,ico) = 0.0
                  cpatch%qmean_rlong_w             (:,ico) = 0.0
                  cpatch%qmean_rad_profile       (:,:,ico) = 0.0
                  cpatch%qmean_sensible_wc         (:,ico) = 0.0
                  cpatch%qmean_vapor_wc            (:,ico) = 0.0
                  cpatch%qmean_intercepted_aw      (:,ico) = 0.0
                  cpatch%qmean_wshed_wg            (:,ico) = 0.0
                  cpatch%qmsqu_gpp                 (:,ico) = 0.0
                  cpatch%qmsqu_npp                 (:,ico) = 0.0
                  cpatch%qmsqu_plresp              (:,ico) = 0.0
                  cpatch%qmsqu_sensible_lc         (:,ico) = 0.0
                  cpatch%qmsqu_vapor_lc            (:,ico) = 0.0
                  cpatch%qmsqu_transp              (:,ico) = 0.0
                  cpatch%qmsqu_sensible_wc         (:,ico) = 0.0
                  cpatch%qmsqu_vapor_wc            (:,ico) = 0.0
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!

      return
   end subroutine zero_ed_qmean_vars
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!










   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !                             |--------------------------------|                        !
   !                             |** YEARLY AVERAGE SUBROUTINES **|                        !
   !                             |--------------------------------|                        !
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This sub-routine updates the yearly variables.                                   !
   !---------------------------------------------------------------------------------------!
   subroutine update_ed_yearly_vars(cgrid)

      use ed_state_vars, only : edtype      & ! structure
                              , polygontype & ! structure
                              , sitetype    & ! structure
                              , patchtype   ! ! structure
      use ed_max_dims  , only : n_pft       & ! intent(in)
                              , n_dbh       ! ! intent(in)
      use consts_coms  , only : pi1         ! ! intent(in)
     
      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      type(edtype)     , target  :: cgrid
      !------ Local variables. ------------------------------------------------------------!
      type(polygontype), pointer :: cpoly
      type(sitetype)   , pointer :: csite
      type(patchtype)  , pointer :: cpatch
      integer                    :: ipy
      integer                    :: isi
      integer                    :: ipa
      integer                    :: ico
      real                       :: poly_area_i
      real                       :: site_area_i
      real                       :: site_wgt
      real                       :: patch_wgt
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      All above-ground biomass variables are in kgC/m2 or kgC/m2/yr; and all basal  !
      ! areas are in cm2/m2 or cm2/m2/yr.                                                  !
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !      Loop over polygons.                                                           !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)


         !----- Re-set the variables. -----------------------------------------------------!
         cgrid%total_basal_area        (ipy) = 0.0
         cgrid%total_basal_area_growth (ipy) = 0.0
         cgrid%total_basal_area_mort   (ipy) = 0.0
         cgrid%total_basal_area_recruit(ipy) = 0.0
         cgrid%total_agb               (ipy) = 0.0
         cgrid%total_agb_growth        (ipy) = 0.0
         cgrid%total_agb_mort          (ipy) = 0.0
         cgrid%total_agb_recruit       (ipy) = 0.0
         !---------------------------------------------------------------------------------!


         !----- Inverse of this polygon area (it should be always 1.) ---------------------!
         poly_area_i = 1./sum(cpoly%area)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !      Loop over polygons.                                                        !
         !---------------------------------------------------------------------------------!
         siteloop: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            !----- Inverse of this site area (it should be always 1.) ---------------------!
            site_area_i=1./sum(csite%area)
            !------------------------------------------------------------------------------!


            !----- Site weight. -----------------------------------------------------------!
            site_wgt = cpoly%area(isi) * poly_area_i
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !      Do growth, mortality, harvesting.                                       !
            !------------------------------------------------------------------------------!
            cgrid%total_agb              (ipy) = cgrid%total_agb              (ipy)        &
                                               + sum(cpoly%agb           (:,:,isi))        &
                                               * site_wgt
            cgrid%total_basal_area       (ipy) = cgrid%total_basal_area       (ipy)        &
                                               + sum(cpoly%basal_area    (:,:,isi))        &
                                               * site_wgt
            cgrid%total_agb_growth       (ipy) = cgrid%total_agb_growth       (ipy)        &
                                               + sum(cpoly%agb_growth    (:,:,isi))        &
                                               * site_wgt
            cgrid%total_agb_mort         (ipy) = cgrid%total_agb_mort         (ipy)        &
                                               + sum(cpoly%agb_mort(1:n_pft,1:n_dbh,isi))  &
                                               * site_wgt
            cgrid%total_basal_area_growth(ipy) = cgrid%total_basal_area_growth(ipy)        &
                                       + sum(cpoly%basal_area_growth(1:n_pft,2:n_dbh,isi)) &
                                       * site_wgt
            cgrid%total_basal_area_mort  (ipy) = cgrid%total_basal_area_mort  (ipy)        &
                                       + sum(cpoly%basal_area_mort  (1:n_pft,2:n_dbh,isi)) &
                                       * site_wgt
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !      Loop over patches and cohorts to get recruitment.                       !
            !------------------------------------------------------------------------------!
            patchloop: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)


               !----- Site weight. --------------------------------------------------------!
               patch_wgt = csite%area(ipa) * site_area_i * site_wgt
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Loop over patches and cohorts to get recruitment.                    !
               !---------------------------------------------------------------------------!
               cohortloop: do ico = 1,cpatch%ncohorts

                  if (cpatch%new_recruit_flag(ico) == 1) then
                     cgrid%total_agb_recruit       (ipy) =                                 &
                            cgrid%total_agb_recruit      (ipy)                             &
                          + cpatch%agb                   (ico)                             &
                          * cpatch%nplant                (ico)                             &
                          * patch_wgt
                     cgrid%total_basal_area_recruit(ipy) =                                 &
                            cgrid%total_basal_area_recruit(ipy)                            &
                          + cpatch%basarea(ico)                                            &
                          * cpatch%nplant (ico)                                            &
                          * patch_wgt
                     cpatch%new_recruit_flag       (ico) = 0
                  end if
                  cpatch%first_census (ico) = 1
               end do cohortloop
               !---------------------------------------------------------------------------!
            end do patchloop
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!
      end do polyloop
      !------------------------------------------------------------------------------------!
      return
   end subroutine update_ed_yearly_vars
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !      This sub-routine re-sets the yearly variables.                                   !
   !---------------------------------------------------------------------------------------!
   subroutine zero_ed_yearly_vars(cgrid)

      use ed_max_dims  , only : n_pft       & ! intent(in)
                              , n_dbh       ! ! intent(in)
      use ed_state_vars, only : edtype      & ! structure
                              , polygontype ! ! structure

      implicit none
      !------ Arguments. ------------------------------------------------------------------!
      type(edtype)     , target  :: cgrid
      !------ Local variables. ------------------------------------------------------------!
      integer                    :: ipy
      integer                    :: isi
      type(polygontype), pointer :: cpoly
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Loop over polygons.                                                            !
      !------------------------------------------------------------------------------------!
      do ipy = 1,cgrid%npolygons
         
         cpoly => cgrid%polygon(ipy)
         !---------------------------------------------------------------------------------!
         !     Loop over sites.                                                            !
         !---------------------------------------------------------------------------------!
         do isi = 1,cpoly%nsites
            cpoly%agb_growth        (:,:,isi) = 0.0
            cpoly%agb_mort          (:,:,isi) = 0.0
            cpoly%agb_cut           (:,:,isi) = 0.0
            cpoly%basal_area_growth (:,:,isi) = 0.0
            cpoly%basal_area_mort   (:,:,isi) = 0.0
            cpoly%basal_area_cut    (:,:,isi) = 0.0
         end do
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!

      return
   end subroutine zero_ed_yearly_vars
   !=======================================================================================!
   !=======================================================================================!




   !=======================================================================================!
   !=======================================================================================!
   !     This function finds the intergrator for sum of squares.  This checks that the     !
   ! squared term does not cause underflow, and if it does, it flushes results to zero.    !
   ! This is very similar to sngloff.                                                      !
   !---------------------------------------------------------------------------------------!
   real(kind=4) function isqu_ftz(x)
      implicit none

      !----- Arguments. -------------------------------------------------------------------!
      real(kind=4), intent(in) :: x
      !----- Local variables. -------------------------------------------------------------!
      real(kind=8)             :: x8
      real(kind=8)             :: xsqu8
      !----- Local parameters. ------------------------------------------------------------!
      real(kind=8), parameter  :: off8 = 1.d-30
      !------------------------------------------------------------------------------------!


      !----- Convert values to double precision. ------------------------------------------!
      x8    = dble(x)
      xsqu8 = x8 * x8
      !------------------------------------------------------------------------------------!

      !----- Check whether overflow is about to happen. -----------------------------------!
      if (abs(xsqu8) < off8) then
         isqu_ftz = 0.
      else
         isqu_ftz = sngl(xsqu8)
      end if
      !------------------------------------------------------------------------------------!

      return
   end function isqu_ftz
   !=======================================================================================!
   !=======================================================================================!
end module average_utils
!==========================================================================================!
!==========================================================================================!
