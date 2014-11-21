!==========================================================================================!
!==========================================================================================!
!     This is the main driver for file output in ED.                                       !
!------------------------------------------------------------------------------------------!
subroutine ed_output(analysis_time,new_day,dail_analy_time,mont_analy_time,dcyc_analy_time &
                    ,annual_time,history_time,dcycle_time,the_end)

   use ed_state_vars, only : edgrid_g                & ! structure
                           , filltab_alltypes        & ! subroutine
                           , filltables              ! ! intent(inout)
   use grid_coms    , only : ngrids                  & ! intent(in)
                           , nzg                     ! ! intent(in)
   use ed_node_coms , only : mynum                   & ! intent(in)
                           , nnodetot                ! ! intent(in)
   use ed_misc_coms , only : dtlsm                   & ! intent(in)
                           , current_time            & ! intent(in)
                           , isoutput                & ! intent(in)
                           , ifoutput                & ! intent(in)
                           , itoutput                & ! intent(in)
                           , writing_dail            & ! intent(in)
                           , writing_mont            & ! intent(in)
                           , writing_dcyc            & ! intent(in)
                           , iprintpolys             & ! intent(in)
                           , frqsum                  ! ! intent(in)
   use average_utils, only : aggregate_polygon_fmean & ! sub-routine
                           , normalize_ed_fmean_vars & ! sub-routine
                           , integrate_ed_dmean_vars & ! sub-routine
                           , integrate_ed_qmean_vars & ! sub-routine
                           , normalize_ed_dmean_vars & ! sub-routine
                           , integrate_ed_mmean_vars & ! sub-routine
                           , zero_ed_dmean_vars      & ! sub-routine
                           , normalize_ed_mmean_vars & ! sub-routine
                           , normalize_ed_qmean_vars & ! sub-routine
                           , zero_ed_mmean_vars      & ! sub-routine
                           , zero_ed_qmean_vars      ! ! sub-routine
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   logical, intent(in)  :: the_end
   logical, intent(in)  :: analysis_time
   logical, intent(in)  :: dail_analy_time
   logical, intent(in)  :: mont_analy_time
   logical, intent(in)  :: dcyc_analy_time
   logical, intent(in)  :: history_time
   logical, intent(in)  :: dcycle_time
   logical, intent(in)  :: new_day
   logical, intent(in)  :: annual_time
   !----- Local variables. ----------------------------------------------------------------!
   integer              :: ifm
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      If there is any IO, then we need to check whether the pointer tables are up to   !
   ! date or not.  They must be rehashed if there has been any change in the number of     !
   ! cohorts or patches (e.g, a cohort or a patch has been terminated or created).         !
   !---------------------------------------------------------------------------------------!
   if (analysis_time   .or. history_time    .or.                                           &
       dail_analy_time .or. mont_analy_time .or. dcyc_analy_time .or. annual_time ) then

      if (filltables) then
         
         !----- Re-hash the tables. -------------------------------------------------------!
         call filltab_alltypes
         !----- Reset the rehash flag. ----------------------------------------------------!
         filltables=.false.
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      If this is the time for an output, we shall call routines to prepare the vari-   !
   ! ables for output.                                                                     !
   !---------------------------------------------------------------------------------------!
   if ( analysis_time .or.   history_time .or. dcycle_time  .or.                           &
       (new_day       .and. (writing_dail .or. writing_mont .or. writing_dcyc))) then
      do ifm=1,ngrids
         call normalize_ed_fmean_vars    (edgrid_g(ifm))
         call aggregate_polygon_fmean    (edgrid_g(ifm))
      end do

      if (writing_dail .or. writing_mont .or. writing_dcyc) then
         do ifm=1,ngrids
            call integrate_ed_dmean_vars(edgrid_g(ifm))
            if (writing_dcyc) call integrate_ed_qmean_vars(edgrid_g(ifm))
         end do
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !----- Instantaneous analysis. ---------------------------------------------------------!
   if (analysis_time) then
      !----- Write out analysis fields - mostly polygon averages. -------------------------!
      if (ifoutput == 3) call h5_output('INST')
      if (itoutput == 3) call h5_output('OPTI')

      !----- If printpolys is on then print this info to the screen. ----------------------!
      if (iprintpolys == 1) then
         do ifm=1,ngrids
            call print_fields(ifm,edgrid_g(ifm))
         end do
      end if
   end if
   !---------------------------------------------------------------------------------------!



   !----- Daily analysis output and monthly integration. ----------------------------------!
   if (new_day .and. (writing_dail .or. writing_mont .or. writing_dcyc)) then

      do ifm=1,ngrids
         call normalize_ed_dmean_vars(edgrid_g(ifm))
         if (writing_mont .or. writing_dcyc) then
            call integrate_ed_mmean_vars(edgrid_g(ifm))
         end if
      end do

      if (dail_analy_time) call h5_output('DAIL')

      do ifm=1,ngrids
         call zero_ed_dmean_vars(edgrid_g(ifm))
      end do
   end if
   !---------------------------------------------------------------------------------------!




   !----- Monthly analysis and monthly mean diurnal cycle output. -------------------------!
   if (mont_analy_time .or. dcyc_analy_time) then
      do ifm=1,ngrids
         call normalize_ed_mmean_vars(edgrid_g(ifm))
         if (writing_dcyc) call normalize_ed_qmean_vars(edgrid_g(ifm))
      end do
      if (mont_analy_time) call h5_output('MONT')
      if (dcyc_analy_time) call h5_output('DCYC')
      do ifm=1,ngrids
         call zero_ed_mmean_vars(edgrid_g(ifm))
         if (writing_dcyc) call zero_ed_qmean_vars(edgrid_g(ifm))
      end do
   end if
   !---------------------------------------------------------------------------------------!



   !----- Yearly analysis output. ---------------------------------------------------------!
   if (annual_time) then
      call h5_output('YEAR')
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !      History files should only be output at a frequency which divides by frqanl, thus !
   ! the integrated fast-time variables are valid, but representative of the last frqanl   !
   ! period, not the last frqhist period.                                                  !
   !---------------------------------------------------------------------------------------!
   if (history_time .and. isoutput /= 0) then
      call h5_output('HIST')
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine ed_output
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     The following subroutine finds the polygon averages from site-, patch-, and cohort-  !
! -level properties that have fmean variables associated.                                  !
!------------------------------------------------------------------------------------------!
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
   !----- Arguments.      -----------------------------------------------------------------!
   type(edtype)         , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype)    , pointer :: cpoly
   type(sitetype)       , pointer :: csite
   type(patchtype)      , pointer :: cpatch
   real, dimension(nzg)           :: cgrid_fmean_soil_hcap
   integer                        :: ipy
   integer                        :: isi
   integer                        :: ipa
   integer                        :: ico
   integer                        :: p
   integer                        :: d
   integer                        :: k
   integer                        :: ksn
   integer                        :: lai_index
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
   real                           :: snowarea
   real                           :: dslzsum_i
   real                           :: rdepth
   real                           :: soil_mstpot
   real                           :: can_exner
   real                           :: atm_exner
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
      ! cgrid%blah(ipy) = 0. ! <<- This way works for all cases. 
      !------------------------------------------------------------------------------------!


      !----- Inverse of this polygon area (it should be always 1.) ------------------------!
      poly_area_i = 1./sum(cpoly%area)
      !------------------------------------------------------------------------------------!

      !----- Re-set some support variables. -----------------------------------------------!
      poly_lai                 = 0.0
      poly_wai                 = 0.0
      poly_nplant              = 0.0
      cgrid_fmean_soil_hcap(:) = 0.0
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

         !----- Inverse of the soil depth. ------------------------------------------------!
         dslzsum_i = 1./ sum(dslz(cpoly%lsl(isi):nzg))
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Loop over patches.                                                          !
         !---------------------------------------------------------------------------------!
         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)


            !----- Site weight. -----------------------------------------------------------!
            patch_wgt = csite%area(ipa) * site_area_i * site_wgt
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !     Loop over cohorts.                                                       !
            !------------------------------------------------------------------------------!
            cohortloop: do ico=1,cpatch%ncohorts
               !----- Aggregate AIs and density, they may be used to normalise averages. --!
               poly_nplant = poly_nplant + cpatch%nplant(ico) * patch_wgt
               poly_lai    = poly_lai    + cpatch%lai   (ico) * patch_wgt
               poly_wai    = poly_wai    + cpatch%wai   (ico) * patch_wgt
               !---------------------------------------------------------------------------!


               !----- Aggregate the other properties. -------------------------------------!
               cgrid%fmean_gpp           (ipy) = cgrid%fmean_gpp            (ipy)          &
                                               + cpatch%fmean_gpp           (ico)          &
                                               * cpatch%nplant              (ico)          &
                                               * patch_wgt
               cgrid%fmean_npp           (ipy) = cgrid%fmean_npp            (ipy)          &
                                               + cpatch%fmean_npp           (ico)          &
                                               * cpatch%nplant              (ico)          &
                                               * patch_wgt
               cgrid%fmean_leaf_resp     (ipy) = cgrid%fmean_leaf_resp      (ipy)          &
                                               + cpatch%fmean_leaf_resp     (ico)          &
                                               * cpatch%nplant              (ico)          &
                                               * patch_wgt
               cgrid%fmean_root_resp     (ipy) = cgrid%fmean_root_resp      (ipy)          &
                                               + cpatch%fmean_root_resp     (ico)          &
                                               * cpatch%nplant              (ico)          &
                                               * patch_wgt
               cgrid%fmean_growth_resp   (ipy) = cgrid%fmean_growth_resp    (ipy)          &
                                               + cpatch%fmean_growth_resp   (ico)          &
                                               * cpatch%nplant              (ico)          &
                                               * patch_wgt
               cgrid%fmean_storage_resp  (ipy) = cgrid%fmean_storage_resp   (ipy)          &
                                               + cpatch%fmean_storage_resp  (ico)          &
                                               * cpatch%nplant              (ico)          &
                                               * patch_wgt
               cgrid%fmean_vleaf_resp    (ipy) = cgrid%fmean_vleaf_resp     (ipy)          &
                                               + cpatch%fmean_vleaf_resp    (ico)          &
                                               * cpatch%nplant              (ico)          &
                                               * patch_wgt
               cgrid%fmean_plresp        (ipy) = cgrid%fmean_plresp         (ipy)          &
                                               + cpatch%fmean_plresp        (ico)          &
                                               * cpatch%nplant              (ico)          &
                                               * patch_wgt
               cgrid%fmean_leaf_energy   (ipy) = cgrid%fmean_leaf_energy    (ipy)          &
                                               + cpatch%fmean_leaf_energy   (ico)          &
                                               * patch_wgt
               cgrid%fmean_leaf_water    (ipy) = cgrid%fmean_leaf_water     (ipy)          &
                                               + cpatch%fmean_leaf_water    (ico)          &
                                               * patch_wgt
               cgrid%fmean_leaf_hcap     (ipy) = cgrid%fmean_leaf_hcap      (ipy)          &
                                               + cpatch%fmean_leaf_hcap     (ico)          &
                                               * patch_wgt
               cgrid%fmean_leaf_vpdef    (ipy) = cgrid%fmean_leaf_vpdef     (ipy)          &
                                               + cpatch%fmean_leaf_vpdef    (ico)          &
                                               * cpatch%lai                 (ico)          &
                                               * patch_wgt
               cgrid%fmean_leaf_gsw      (ipy) = cgrid%fmean_leaf_gsw       (ipy)          &
                                               + cpatch%fmean_leaf_gsw      (ico)          &
                                               * cpatch%lai                 (ico)          &
                                               * patch_wgt
               cgrid%fmean_leaf_gbw      (ipy) = cgrid%fmean_leaf_gbw       (ipy)          &
                                               + cpatch%fmean_leaf_gbw      (ico)          &
                                               * cpatch%lai                 (ico)          &
                                               * patch_wgt
               cgrid%fmean_wood_energy   (ipy) = cgrid%fmean_wood_energy    (ipy)          &
                                               + cpatch%fmean_wood_energy   (ico)          &
                                               * patch_wgt
               cgrid%fmean_wood_water    (ipy) = cgrid%fmean_wood_water     (ipy)          &
                                               + cpatch%fmean_wood_water    (ico)          &
                                               * patch_wgt
               cgrid%fmean_wood_hcap     (ipy) = cgrid%fmean_wood_hcap      (ipy)          &
                                               + cpatch%fmean_wood_hcap     (ico)          &
                                               * patch_wgt
               cgrid%fmean_wood_gbw      (ipy) = cgrid%fmean_wood_gbw       (ipy)          &
                                               + cpatch%fmean_wood_gbw      (ico)          &
                                               * cpatch%wai                 (ico)          &
                                               * patch_wgt
               cgrid%fmean_fs_open       (ipy) = cgrid%fmean_fs_open        (ipy)          &
                                               + cpatch%fmean_fs_open       (ico)          &
                                               * cpatch%lai                 (ico)          &
                                               * patch_wgt
               cgrid%fmean_fsw           (ipy) = cgrid%fmean_fsw            (ipy)          &
                                               + cpatch%fmean_fsw           (ico)          &
                                               * cpatch%lai                 (ico)          &
                                               * patch_wgt
               cgrid%fmean_fsn           (ipy) = cgrid%fmean_fsn            (ipy)          &
                                               + cpatch%fmean_fsn           (ico)          &
                                               * cpatch%lai                 (ico)          &
                                               * patch_wgt
               cgrid%fmean_psi_open      (ipy) = cgrid%fmean_psi_open       (ipy)          &
                                               + cpatch%fmean_psi_open      (ico)          &
                                               * cpatch%lai                 (ico)          &
                                               * patch_wgt
               cgrid%fmean_psi_closed    (ipy) = cgrid%fmean_psi_closed     (ipy)          &
                                               + cpatch%fmean_psi_closed    (ico)          &
                                               * cpatch%lai                 (ico)          &
                                               * patch_wgt
               cgrid%fmean_water_supply  (ipy) = cgrid%fmean_water_supply   (ipy)          &
                                               + cpatch%fmean_water_supply  (ico)          &
                                               * patch_wgt
               cgrid%fmean_par_l         (ipy) = cgrid%fmean_par_l          (ipy)          &
                                               + cpatch%fmean_par_l         (ico)          &
                                               * patch_wgt
               cgrid%fmean_par_l_beam    (ipy) = cgrid%fmean_par_l_beam     (ipy)          &
                                               + cpatch%fmean_par_l_beam    (ico)          &
                                               * patch_wgt
               cgrid%fmean_par_l_diff    (ipy) = cgrid%fmean_par_l_diff     (ipy)          &
                                               + cpatch%fmean_par_l_diff    (ico)          &
                                               * patch_wgt
               cgrid%fmean_rshort_l      (ipy) = cgrid%fmean_rshort_l       (ipy)          &
                                               + cpatch%fmean_rshort_l      (ico)          &
                                               * patch_wgt
               cgrid%fmean_rlong_l       (ipy) = cgrid%fmean_rlong_l        (ipy)          &
                                               + cpatch%fmean_rlong_l       (ico)          &
                                               * patch_wgt
               cgrid%fmean_sensible_lc   (ipy) = cgrid%fmean_sensible_lc    (ipy)          &
                                               + cpatch%fmean_sensible_lc   (ico)          &
                                               * patch_wgt
               cgrid%fmean_vapor_lc      (ipy) = cgrid%fmean_vapor_lc       (ipy)          &
                                               + cpatch%fmean_vapor_lc      (ico)          &
                                               * patch_wgt
               cgrid%fmean_transp        (ipy) = cgrid%fmean_transp         (ipy)          &
                                               + cpatch%fmean_transp        (ico)          &
                                               * patch_wgt
               cgrid%fmean_intercepted_al(ipy) = cgrid%fmean_intercepted_al (ipy)          &
                                               + cpatch%fmean_intercepted_al(ico)          &
                                               * patch_wgt
               cgrid%fmean_wshed_lg      (ipy) = cgrid%fmean_wshed_lg       (ipy)          &
                                               + cpatch%fmean_wshed_lg      (ico)          &
                                               * patch_wgt
               cgrid%fmean_rshort_w      (ipy) = cgrid%fmean_rshort_w       (ipy)          &
                                               + cpatch%fmean_rshort_w      (ico)          &
                                               * patch_wgt
               cgrid%fmean_rlong_w       (ipy) = cgrid%fmean_rlong_w        (ipy)          &
                                               + cpatch%fmean_rlong_w       (ico)          &
                                               * patch_wgt
               cgrid%fmean_sensible_wc   (ipy) = cgrid%fmean_sensible_wc    (ipy)          &
                                               + cpatch%fmean_sensible_wc   (ico)          &
                                               * patch_wgt
               cgrid%fmean_vapor_wc      (ipy) = cgrid%fmean_vapor_wc       (ipy)          &
                                               + cpatch%fmean_vapor_wc      (ico)          &
                                               * patch_wgt
               cgrid%fmean_intercepted_aw(ipy) = cgrid%fmean_intercepted_aw (ipy)          &
                                               + cpatch%fmean_intercepted_aw(ico)          &
                                               * patch_wgt
               cgrid%fmean_wshed_wg      (ipy) = cgrid%fmean_wshed_wg       (ipy)          &
                                               + cpatch%fmean_wshed_wg      (ico)          &
                                               * patch_wgt
            end do cohortloop
            !------------------------------------------------------------------------------!




            !------------------------------------------------------------------------------!
            !     Aggregate the patch-level variables.                                     !
            !------------------------------------------------------------------------------!
            cgrid%fmean_rh             (ipy) = cgrid%fmean_rh             (ipy)            &
                                             + csite%fmean_rh             (ipa)            &
                                             * patch_wgt
            cgrid%fmean_cwd_rh         (ipy) = cgrid%fmean_cwd_rh         (ipy)            &
                                             + csite%fmean_cwd_rh         (ipa)            &
                                             * patch_wgt
            cgrid%fmean_nep            (ipy) = cgrid%fmean_nep            (ipy)            &
                                             + csite%fmean_nep            (ipa)            &
                                             * patch_wgt
            cgrid%fmean_rk4step        (ipy) = cgrid%fmean_rk4step        (ipy)            &
                                             + csite%fmean_rk4step        (ipa)            &
                                             * patch_wgt
            cgrid%fmean_available_water(ipy) = cgrid%fmean_available_water(ipy)            &
                                             + csite%fmean_available_water(ipa)            &
                                             * patch_wgt
            cgrid%fmean_can_theiv      (ipy) = cgrid%fmean_can_theiv      (ipy)            &
                                             + csite%fmean_can_theiv      (ipa)            &
                                             * patch_wgt
            cgrid%fmean_can_theta      (ipy) = cgrid%fmean_can_theta      (ipy)            &
                                             + csite%fmean_can_theta      (ipa)            &
                                             * patch_wgt
            cgrid%fmean_can_vpdef      (ipy) = cgrid%fmean_can_vpdef      (ipy)            &
                                             + csite%fmean_can_vpdef      (ipa)            &
                                             * patch_wgt
            cgrid%fmean_can_shv        (ipy) = cgrid%fmean_can_shv        (ipy)            &
                                             + csite%fmean_can_shv        (ipa)            &
                                             * patch_wgt
            cgrid%fmean_can_co2        (ipy) = cgrid%fmean_can_co2        (ipy)            &
                                             + csite%fmean_can_co2        (ipa)            &
                                             * patch_wgt
            cgrid%fmean_can_prss       (ipy) = cgrid%fmean_can_prss       (ipy)            &
                                             + csite%fmean_can_prss       (ipa)            &
                                             * patch_wgt
            cgrid%fmean_gnd_temp       (ipy) = cgrid%fmean_gnd_temp       (ipy)            &
                                             + csite%fmean_gnd_temp       (ipa)            &
                                             * patch_wgt
            cgrid%fmean_gnd_shv        (ipy) = cgrid%fmean_gnd_shv        (ipy)            &
                                             + csite%fmean_gnd_shv        (ipa)            &
                                             * patch_wgt
            cgrid%fmean_can_ggnd       (ipy) = cgrid%fmean_can_ggnd       (ipy)            &
                                             + csite%fmean_can_ggnd       (ipa)            &
                                             * patch_wgt
            cgrid%fmean_sfcw_depth     (ipy) = cgrid%fmean_sfcw_depth     (ipy)            &
                                             + csite%fmean_sfcw_depth     (ipa)            &
                                             * patch_wgt
            !----- Temporarily convert pounding internal energy to J/m2. ------------------!
            cgrid%fmean_sfcw_energy    (ipy) = cgrid%fmean_sfcw_energy    (ipy)            &
                                             + csite%fmean_sfcw_energy    (ipa)            &
                                             * csite%fmean_sfcw_mass      (ipa)            &
                                             * patch_wgt
            !------------------------------------------------------------------------------!
            cgrid%fmean_sfcw_mass      (ipy) = cgrid%fmean_sfcw_mass      (ipy)            &
                                             + csite%fmean_sfcw_mass      (ipa)            &
                                             * patch_wgt
            cgrid%fmean_rshort_gnd     (ipy) = cgrid%fmean_rshort_gnd     (ipy)            &
                                             + csite%fmean_rshort_gnd     (ipa)            &
                                             * patch_wgt
            cgrid%fmean_par_gnd        (ipy) = cgrid%fmean_par_gnd        (ipy)            &
                                             + csite%fmean_par_gnd        (ipa)            &
                                             * patch_wgt
            cgrid%fmean_rlong_gnd      (ipy) = cgrid%fmean_rlong_gnd      (ipy)            &
                                             + csite%fmean_rlong_gnd      (ipa)            &
                                             * patch_wgt
            cgrid%fmean_rlongup        (ipy) = cgrid%fmean_rlongup        (ipy)            &
                                             + csite%fmean_rlongup        (ipa)            &
                                             * patch_wgt
            cgrid%fmean_parup          (ipy) = cgrid%fmean_parup          (ipy)            &
                                             + csite%fmean_parup          (ipa)            &
                                             * patch_wgt
            cgrid%fmean_nirup          (ipy) = cgrid%fmean_nirup          (ipy)            &
                                             + csite%fmean_nirup          (ipa)            &
                                             * patch_wgt
            cgrid%fmean_rshortup       (ipy) = cgrid%fmean_rshortup       (ipy)            &
                                             + csite%fmean_rshortup       (ipa)            &
                                             * patch_wgt
            cgrid%fmean_rnet           (ipy) = cgrid%fmean_rnet           (ipy)            &
                                             + csite%fmean_rnet           (ipa)            &
                                             * patch_wgt
            cgrid%fmean_albedo         (ipy) = cgrid%fmean_albedo         (ipy)            &
                                             + csite%fmean_albedo         (ipa)            &
                                             * patch_wgt
            cgrid%fmean_albedo_par     (ipy) = cgrid%fmean_albedo_par     (ipy)            &
                                             + csite%fmean_albedo_par     (ipa)            &
                                             * patch_wgt
            cgrid%fmean_albedo_nir     (ipy) = cgrid%fmean_albedo_nir     (ipy)            &
                                             + csite%fmean_albedo_nir     (ipa)            &
                                             * patch_wgt
            cgrid%fmean_rlong_albedo   (ipy) = cgrid%fmean_rlong_albedo   (ipy)            &
                                             + csite%fmean_rlong_albedo   (ipa)            &
                                             * patch_wgt
            cgrid%fmean_ustar          (ipy) = cgrid%fmean_ustar          (ipy)            &
                                             + csite%fmean_ustar          (ipa)            &
                                             * patch_wgt
            cgrid%fmean_tstar          (ipy) = cgrid%fmean_tstar          (ipy)            &
                                             + csite%fmean_tstar          (ipa)            &
                                             * patch_wgt
            cgrid%fmean_qstar          (ipy) = cgrid%fmean_qstar          (ipy)            &
                                             + csite%fmean_qstar          (ipa)            &
                                             * patch_wgt
            cgrid%fmean_cstar          (ipy) = cgrid%fmean_cstar          (ipy)            &
                                             + csite%fmean_cstar          (ipa)            &
                                             * patch_wgt
            cgrid%fmean_carbon_ac      (ipy) = cgrid%fmean_carbon_ac      (ipy)            &
                                             + csite%fmean_carbon_ac      (ipa)            &
                                             * patch_wgt
            cgrid%fmean_carbon_st      (ipy) = cgrid%fmean_carbon_st      (ipy)            &
                                             + csite%fmean_carbon_st      (ipa)            &
                                             * patch_wgt
            cgrid%fmean_vapor_gc       (ipy) = cgrid%fmean_vapor_gc       (ipy)            &
                                             + csite%fmean_vapor_gc       (ipa)            &
                                             * patch_wgt
            cgrid%fmean_vapor_ac       (ipy) = cgrid%fmean_vapor_ac       (ipy)            &
                                             + csite%fmean_vapor_ac       (ipa)            &
                                             * patch_wgt
            cgrid%fmean_throughfall    (ipy) = cgrid%fmean_throughfall    (ipy)            &
                                             + csite%fmean_throughfall    (ipa)            &
                                             * patch_wgt
            cgrid%fmean_runoff         (ipy) = cgrid%fmean_runoff         (ipy)            &
                                             + csite%fmean_runoff         (ipa)            &
                                             * patch_wgt
            cgrid%fmean_drainage       (ipy) = cgrid%fmean_drainage       (ipy)            &
                                             + csite%fmean_drainage       (ipa)            &
                                             * patch_wgt
            cgrid%fmean_sensible_gc    (ipy) = cgrid%fmean_sensible_gc    (ipy)            &
                                             + csite%fmean_sensible_gc    (ipa)            &
                                             * patch_wgt
            cgrid%fmean_sensible_ac    (ipy) = cgrid%fmean_sensible_ac    (ipy)            &
                                             + csite%fmean_sensible_ac    (ipa)            &
                                             * patch_wgt
            cgrid%fmean_qthroughfall   (ipy) = cgrid%fmean_qthroughfall   (ipy)            &
                                             + csite%fmean_qthroughfall   (ipa)            &
                                             * patch_wgt
            cgrid%fmean_qrunoff        (ipy) = cgrid%fmean_qrunoff        (ipy)            &
                                             + csite%fmean_qrunoff        (ipa)            &
                                             * patch_wgt
            cgrid%fmean_qdrainage      (ipy) = cgrid%fmean_qdrainage      (ipy)            &
                                             + csite%fmean_qdrainage      (ipa)            &
                                             * patch_wgt

            !----- Soil (extensive) properties. -------------------------------------------!
            do k=1,nzg
               nsoil = cpoly%ntext_soil(k,isi)
               cgrid%fmean_soil_energy (k,ipy) = cgrid%fmean_soil_energy (k,ipy)           &
                                               + csite%fmean_soil_energy (k,ipa)           &
                                               * patch_wgt
               cgrid%fmean_soil_mstpot (k,ipy) = cgrid%fmean_soil_mstpot (k,ipy)           &
                                               + csite%fmean_soil_mstpot (k,ipa)           &
                                               * patch_wgt
               cgrid%fmean_soil_water  (k,ipy) = cgrid%fmean_soil_water  (k,ipy)           &
                                               + csite%fmean_soil_water  (k,ipa)           &
                                               * patch_wgt
               cgrid%fmean_smoist_gg   (k,ipy) = cgrid%fmean_smoist_gg   (k,ipy)           &
                                               + csite%fmean_smoist_gg   (k,ipa)           &
                                               * patch_wgt
               cgrid%fmean_transloss   (k,ipy) = cgrid%fmean_transloss   (k,ipy)           &
                                               + csite%fmean_transloss   (k,ipa)           &
                                               * patch_wgt
               cgrid%fmean_sensible_gg (k,ipy) = cgrid%fmean_sensible_gg (k,ipy)           &
                                               + csite%fmean_sensible_gg (k,ipa)           &
                                               * patch_wgt
               cgrid_fmean_soil_hcap   (k)     = cgrid_fmean_soil_hcap   (k)               &
                                               + soil(nsoil)%slcpd                         &
                                               * patch_wgt

               !----- Find the mean soil wetness. -----------------------------------------!
               cgrid%fmean_soil_wetness  (ipy) = cgrid%fmean_soil_wetness    (ipy)         &
                                               + ( ( csite%fmean_soil_water(k,ipa)         &
                                                   - soil(nsoil)%soilwp) )                 &
                                               / (soil(nsoil)%slmsts - soil(nsoil)%soilwp) &
                                               * dslz(k) * dslzsum_i * patch_wgt
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!
         end do patchloop
         !---------------------------------------------------------------------------------!




         !---------------------------------------------------------------------------------!
         !     Aggregate the patch-level variables.                                        !
         !---------------------------------------------------------------------------------!
         cgrid%fmean_atm_theiv      (ipy) = cgrid%fmean_atm_theiv      (ipy)               &
                                          + cpoly%fmean_atm_theiv      (isi)               &
                                          * site_wgt
         cgrid%fmean_atm_theta      (ipy) = cgrid%fmean_atm_theta      (ipy)               &
                                          + cpoly%fmean_atm_theta      (isi)               &
                                          * site_wgt
         cgrid%fmean_atm_temp       (ipy) = cgrid%fmean_atm_temp       (ipy)               &
                                          + cpoly%fmean_atm_temp       (isi)               &
                                          * site_wgt
         cgrid%fmean_atm_vpdef      (ipy) = cgrid%fmean_atm_vpdef      (ipy)               &
                                          + cpoly%fmean_atm_vpdef      (isi)               &
                                          * site_wgt
         cgrid%fmean_atm_shv        (ipy) = cgrid%fmean_atm_shv        (ipy)               &
                                          + cpoly%fmean_atm_shv        (isi)               &
                                          * site_wgt
         cgrid%fmean_atm_rshort     (ipy) = cgrid%fmean_atm_rshort     (ipy)               &
                                          + cpoly%fmean_atm_rshort     (isi)               &
                                          * site_wgt
         cgrid%fmean_atm_rshort_diff(ipy) = cgrid%fmean_atm_rshort_diff(ipy)               &
                                          + cpoly%fmean_atm_rshort_diff(isi)               &
                                          * site_wgt
         cgrid%fmean_atm_par        (ipy) = cgrid%fmean_atm_par        (ipy)               &
                                          + cpoly%fmean_atm_par        (isi)               &
                                          * site_wgt
         cgrid%fmean_atm_par_diff   (ipy) = cgrid%fmean_atm_par_diff   (ipy)               &
                                          + cpoly%fmean_atm_par_diff   (isi)               &
                                          * site_wgt
         cgrid%fmean_atm_rlong      (ipy) = cgrid%fmean_atm_rlong      (ipy)               &
                                          + cpoly%fmean_atm_rlong      (isi)               &
                                          * site_wgt
         cgrid%fmean_atm_vels       (ipy) = cgrid%fmean_atm_vels       (ipy)               &
                                          + cpoly%fmean_atm_vels       (isi)               &
                                          * site_wgt
         cgrid%fmean_atm_rhos       (ipy) = cgrid%fmean_atm_rhos       (ipy)               &
                                          + cpoly%fmean_atm_rhos       (isi)               &
                                          * site_wgt
         cgrid%fmean_atm_prss       (ipy) = cgrid%fmean_atm_prss       (ipy)               &
                                          + cpoly%fmean_atm_prss       (isi)               &
                                          * site_wgt
         cgrid%fmean_atm_co2        (ipy) = cgrid%fmean_atm_co2        (ipy)               &
                                          + cpoly%fmean_atm_co2        (isi)               &
                                          * site_wgt
         cgrid%fmean_pcpg           (ipy) = cgrid%fmean_pcpg           (ipy)               &
                                          + cpoly%fmean_pcpg           (isi)               &
                                          * site_wgt
         cgrid%fmean_qpcpg          (ipy) = cgrid%fmean_qpcpg          (ipy)               &
                                          + cpoly%fmean_qpcpg          (isi)               &
                                          * site_wgt
         cgrid%fmean_dpcpg          (ipy) = cgrid%fmean_dpcpg          (ipy)               &
                                          + cpoly%fmean_dpcpg          (isi)               &
                                          * site_wgt
      end do siteloop
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Find the derived properties for the air above canopy.                         !
      !------------------------------------------------------------------------------------!
      atm_exner                 = press2exner (cgrid%fmean_atm_prss(ipy))
      cgrid%fmean_atm_temp(ipy) = extheta2temp(atm_exner,cgrid%fmean_atm_theta(ipy))
      cgrid%fmean_atm_rhos(ipy) = idealdenssh ( cgrid%fmean_atm_prss  (ipy)                &
                                              , cgrid%fmean_atm_temp  (ipy)                &
                                              , cgrid%fmean_atm_shv   (ipy) )
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Find the derived properties for the canopy air space.                         !
      !------------------------------------------------------------------------------------!
      can_exner                 = press2exner (cgrid%fmean_can_prss(ipy))
      cgrid%fmean_can_temp(ipy) = extheta2temp(can_exner,cgrid%fmean_can_theta(ipy))
      cgrid%fmean_can_rhos(ipy) = idealdenssh ( cgrid%fmean_can_prss  (ipy)                &
                                              , cgrid%fmean_can_temp  (ipy)                &
                                              , cgrid%fmean_can_shv   (ipy) )
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !   If the patch had some temporary snow/pounding layer, convert the mean energy to  !
      ! J/kg, then find the mean temperature and liquid fraction.  Otherwise, set them to  !
      ! either zero or default values.                                                     !
      !------------------------------------------------------------------------------------!
      if (cgrid%fmean_sfcw_mass(ipy) > tiny_sfcwater_mass) then
         cgrid%fmean_sfcw_energy(ipy) = cgrid%fmean_sfcw_energy(ipy)                       &
                                      / cgrid%fmean_sfcw_mass(ipy)
         call uint2tl(cgrid%fmean_sfcw_energy(ipy),cgrid%fmean_sfcw_temp(ipy)              &
                     ,cgrid%fmean_sfcw_fliq(ipy))
      else
         cgrid%fmean_sfcw_mass  (ipy)  = 0.
         cgrid%fmean_sfcw_depth (ipy)  = 0.
         cgrid%fmean_sfcw_energy(ipy)  = 0.
         cgrid%fmean_sfcw_temp  (ipy)  = cgrid%fmean_soil_temp(nzg,ipy)
         cgrid%fmean_sfcw_fliq  (ipy)  = cgrid%fmean_soil_fliq(nzg,ipy)
      end if
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Find the temperature and the fraction of liquid water.                         !
      !------------------------------------------------------------------------------------!
      do k=1,nzg
         call uextcm2tl(cgrid%fmean_soil_energy(k,ipy),cgrid%fmean_soil_water(k,ipy)*wdns  &
                       ,cgrid_fmean_soil_hcap  (k)    ,cgrid%fmean_soil_temp (k,ipy)       &
                       ,cgrid%fmean_soil_fliq  (k,ipy))
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the vegetation temperature and liquid fraction.                           !
      !------------------------------------------------------------------------------------!
      !----- Leaf. ------------------------------------------------------------------------!
      if (cgrid%fmean_leaf_hcap(ipy) > 0.) then
         call uextcm2tl( cgrid%fmean_leaf_energy(ipy), cgrid%fmean_leaf_water (ipy)        &
                       , cgrid%fmean_leaf_hcap  (ipy), cgrid%fmean_leaf_temp  (ipy)        &
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
      !----- Wood. ------------------------------------------------------------------------!
      if (cgrid%fmean_wood_hcap(ipy) > 0.) then
         call uextcm2tl( cgrid%fmean_wood_energy(ipy)                                     &
                       , cgrid%fmean_wood_water (ipy)                                     &
                       , cgrid%fmean_wood_hcap  (ipy)                                     &
                       , cgrid%fmean_wood_temp  (ipy)                                     &
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
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !    Normalise the "intensive" properties.  The weight was either the LAI, WAI, or   !
      ! plant density.  In case none of the cohorts qualified to contribute, then we       !
      ! assign either the canopy air space property, or a default number.                  !
      !------------------------------------------------------------------------------------!
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
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !    Compute the average amongst all surfaces (soil, temporary surface water, and    !
      ! vegetation, the last two only if they really exist).  All energy terms are         !
      ! converted to J/m2, all water terms to kg/m2, and the heat capacities of everything !
      ! that is not water is in J/m2/K.                                                    !
      !------------------------------------------------------------------------------------!
      skin_energy = cgrid%fmean_leaf_energy(ipy) + cgrid%fmean_wood_energy(ipy)            &
                  + cgrid%fmean_sfcw_energy(ipy) * cgrid%fmean_sfcw_mass  (ipy)            &
                  + cgrid%fmean_soil_energy(nzg,ipy) * dslz(nzg)
      skin_water  = cgrid%fmean_leaf_water(ipy)  + cgrid%fmean_wood_water(ipy)             &
                  + cgrid%fmean_sfcw_mass (ipy)                                            &
                  + cgrid%fmean_soil_water(nzg,ipy) * dslz(nzg) * wdns
      skin_hcap   = cgrid%fmean_leaf_hcap(ipy) + cgrid%fmean_wood_hcap(ipy)                &
                  + cgrid_fmean_soil_hcap(nzg) * dslz(nzg)
      call uextcm2tl(skin_energy,skin_water,skin_hcap,cgrid%fmean_skin_temp(ipy),skin_fliq)
      !------------------------------------------------------------------------------------!
   end do polyloop
   !---------------------------------------------------------------------------------------!
   return
end subroutine aggregate_polygon_fmean
!==========================================================================================!
!==========================================================================================!
