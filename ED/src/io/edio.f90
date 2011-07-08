!==========================================================================================!
!==========================================================================================!
!     This is the main driver for file output in ED.                                       !
!------------------------------------------------------------------------------------------!
subroutine ed_output(analysis_time,new_day,dail_analy_time,mont_analy_time,dcyc_analy_time &
                    ,annual_time,writing_dail,writing_mont,writing_dcyc,history_time       &
                    ,dcycle_time,the_end)

   use ed_state_vars, only : edgrid_g          & ! structure
                           , filltab_alltypes  & ! subroutine
                           , filltables        ! ! intent(inout)
   use grid_coms    , only : ngrids            & ! intent(in)
                           , nzg               ! ! intent(in)
   use ed_node_coms , only : mynum             & ! intent(in)
                           , nnodetot          ! ! intent(in)
   use ed_misc_coms , only : dtlsm             & ! intent(in)
                           , current_time      & ! intent(in)
                           , idoutput          & ! intent(in)
                           , imoutput          & ! intent(in)
                           , iqoutput          & ! intent(in)
                           , iyoutput          & ! intent(in)
                           , isoutput          & ! intent(in)
                           , ifoutput          & ! intent(in)
                           , itoutput          & ! intent(in)
                           , iprintpolys       & ! intent(in)
                           , frqsum            ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   logical, intent(in)  :: the_end
   logical, intent(in)  :: analysis_time
   logical, intent(in)  :: dail_analy_time
   logical, intent(in)  :: mont_analy_time
   logical, intent(in)  :: dcyc_analy_time
   logical, intent(in)  :: writing_dail
   logical, intent(in)  :: writing_mont
   logical, intent(in)  :: writing_dcyc
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
         call normalize_averaged_vars(edgrid_g(ifm),frqsum,dtlsm)
      end do

      !----- Perform averaging and data preparation. --------------------------------------!
      call spatial_averages
      
      if (writing_dail .or. writing_mont .or. writing_dcyc) then
         do ifm=1,ngrids
            call integrate_ed_daily_output_flux(edgrid_g(ifm))
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
      call avg_ed_daily_output_pool()

      do ifm=1,ngrids
         call normalize_ed_daily_output_vars(edgrid_g(ifm))
         if (writing_mont .or. writing_dcyc) then
            call integrate_ed_monthly_output_vars(edgrid_g(ifm))
         end if
      end do

      if (dail_analy_time) call h5_output('DAIL')

      do ifm=1,ngrids
         call zero_ed_daily_output_vars(edgrid_g(ifm))
      end do
   end if
   !---------------------------------------------------------------------------------------!




   !----- Monthly analysis and monthly mean diurnal cycle output. -------------------------!
   if (mont_analy_time .or. dcyc_analy_time) then
      do ifm=1,ngrids
         call normalize_ed_monthly_output_vars(edgrid_g(ifm))
      end do
      if (mont_analy_time) call h5_output('MONT')
      if (dcyc_analy_time) call h5_output('DCYC')
      do ifm=1,ngrids
         call zero_ed_monthly_output_vars(edgrid_g(ifm))
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
!      This subroutine calculates the polygon average of the carbon and nitrogen POOLS for !
! outputting at a DAILY timestep.                                                          !
!------------------------------------------------------------------------------------------!
subroutine avg_ed_daily_output_pool()
   use ed_state_vars, only : edtype       & ! structure
                           , polygontype  & ! structure
                           , sitetype     & ! structure
                           , patchtype    & ! structure
                           , edgrid_g     ! ! structure
   use grid_coms    , only : ngrids       & ! intent(in)
                           , nzg          & ! intent(in)
                           , nzs          ! ! intent(in)
   use pft_coms     , only : c2n_leaf     & ! intent(in)
                           , c2n_stem     & ! intent(in)
                           , c2n_storage  ! ! intent(in)
   implicit none

   !----- Local variables. ----------------------------------------------------------------!
   type(edtype)     , pointer :: cgrid
   type(polygontype), pointer :: cpoly
   type(sitetype)   , pointer :: csite
   type(patchtype)  , pointer :: cpatch
   integer                    :: igr
   integer                    :: ipy
   integer                    :: isi
   integer                    :: ipa
   integer                    :: ico
   integer                    :: ipft
   real                       :: area_si
   real                       :: area_pa
   !---------------------------------------------------------------------------------------!

   gridloop: do igr=1,ngrids
      cgrid => edgrid_g(igr)

      !------------------------------------------------------------------------------------!
      !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !
      !------------------------------------------------------------------------------------!
      !     Please, don't initialise polygon-level (cgrid) variables outside polyloop.     !
      ! This works in off-line runs, but it causes memory leaks (and crashes) in the       !
      ! coupled runs over the ocean, where cgrid%npolygons can be 0 if one of the sub-     !
      ! domains falls entirely over the ocean.  Thanks!                                    !
      !------------------------------------------------------------------------------------!
      ! cgrid%blah = 0. !<<--- This is a bad way of doing, look inside the loop for the
      !                 !      safe way of initialising the variable.
      !------------------------------------------------------------------------------------!
      polygonloop: do ipy=1,cgrid%npolygons
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
         !----- Zero variables. -----------------------------------------------------------!
         cgrid%Cleaf (ipy)  = 0.0
         cgrid%Croot (ipy)  = 0.0
         cgrid%Cstore(ipy)  = 0.0
         cgrid%Ccwd  (ipy)  = -9999.  !! haven't figured out simple way to get this yet
         cgrid%Nleaf (ipy)  = 0.0
         cgrid%Ndead (ipy)  = 0.0
         cgrid%Nroot (ipy)  = 0.0
         cgrid%Nstore(ipy)  = 0.0
         cgrid%Ncwd  (ipy)  = -9999.  !! haven't figured out simple way to get this yet
                                      !! b.c. CWD is an implict part of stsc

         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            area_si = cpoly%area(isi)
            patchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)
               area_pa = area_si * csite%area(ipa)
               !---------------------------------------------------------------------------!
               !     Here we must include a loop through all cohorts, because this may be  !
               ! an empty patch and vector operations cannot be done if the patchtype      !
               ! structure is not allocated.  This actually happens in both off-line and   !
               ! coupled runs, especially over deserts...                                  !
               !---------------------------------------------------------------------------!
               cohortloop: do ico = 1,cpatch%ncohorts
                  ipft = cpatch%pft(ico)

                  cgrid%Cleaf(ipy)   = cgrid%Cleaf(ipy)                                    &
                                     + cpatch%bleaf(ico) * cpatch%nplant(ico) * area_pa
                  cgrid%Cstore(ipy)  = cgrid%Cstore(ipy)                                   &
                                     + cpatch%bstorage(ico) * cpatch%nplant(ico) * area_pa
                  cgrid%Croot(ipy)   = cgrid%Croot(ipy)                                    &
                                     + cpatch%broot(ico) * cpatch%nplant(ico) * area_pa

                  cgrid%Nleaf(ipy)   = cgrid%Nleaf(ipy)                                    &
                                     + cpatch%bleaf(ico) * cpatch%nplant(ico)              &
                                     / c2n_leaf(ipft) * area_pa
                  cgrid%Nstore(ipy)  = cgrid%Nstore(ipy)                                   &
                                     * cpatch%bstorage(ico) * cpatch%nplant(ico)           &
                                     / c2n_storage * area_pa ! C:N not pft specific
                  !----- It appears we assume leaf and root have same C:N. ----------------!
                  cgrid%Nroot(ipy)   = cgrid%Nroot(ipy)                                    &
                                     + cpatch%broot(ico) * cpatch%nplant(ico)              &
                                     / c2n_leaf(ipft) * area_pa 
                  cgrid%Ndead(ipy)   = cgrid%Ndead(ipy)                                    &
                                     + cpatch%bdead(ico) * cpatch%nplant(ico)              &
                                     / c2n_stem(ipft) * area_pa
               end do cohortloop
            end do patchloop
         end do siteloop
      end do polygonloop
   end do gridloop
   return
end subroutine avg_ed_daily_output_pool
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     The following subroutine performs several spatial averaging functions and temporal   !
! integrations.  Specifically, it area averages patch level quantities to the site, and    !
! site level quantities to the polygon.                                                    !
!------------------------------------------------------------------------------------------!
subroutine spatial_averages

   use ed_state_vars         , only : edtype             & ! structure
                                    , polygontype        & ! structure
                                    , sitetype           & ! structure
                                    , patchtype          & ! structure
                                    , edgrid_g           ! ! structure
   use grid_coms             , only : ngrids             & ! intent(in)
                                    , nzg                & ! intent(in)
                                    , nzs                ! ! intent(in)
   use consts_coms           , only : alvl               & ! intent(in)
                                    , cpi                & ! intent(in)
                                    , wdns               & ! intent(in)
                                    , p00i               & ! intent(in)
                                    , t00                & ! intent(in)
                                    , rocp               & ! intent(in)
                                    , umol_2_kgC         & ! intent(in)
                                    , day_sec            ! ! intent(in)
   use ed_misc_coms          , only : frqsum             ! ! intent(in)
   use therm_lib             , only : qwtk               & ! subroutine
                                    , qtk                & ! subroutine
                                    , idealdenssh        ! ! function
   use soil_coms             , only : tiny_sfcwater_mass & ! intent(in)
                                    , isoilbc            & ! intent(in)
                                    , soil               & ! intent(in)
                                    , dslz               ! ! intent(in)
   use c34constants          , only : n_stoma_atts       ! ! intent(in)
   use ed_max_dims           , only : n_pft              ! ! intent(in)
   implicit none
   !----- Local variables -----------------------------------------------------------------!
   type(edtype)         , pointer :: cgrid
   type(polygontype)    , pointer :: cpoly
   type(sitetype)       , pointer :: csite
   type(patchtype)      , pointer :: cpatch
   real, dimension(3)             :: area_sum
   real, dimension(nzg)           :: site_avg_soil_hcap
   real, dimension(nzg)           :: poly_avg_soil_hcap
   integer                        :: igr,ipy,isi,ipa,ico,ipft,iatt
   integer                        :: k,ksn
   integer                        :: lai_index
   integer                        :: nsoil
   real                           :: lai_patch
   real                           :: laiarea_site
   real                           :: laiarea_poly
   real                           :: site_area_i
   real                           :: poly_area_i
   real                           :: frqsumi
   real                           :: skin_energy
   real                           :: skin_water
   real                           :: skin_hcap
   real                           :: skin_fliq
   real                           :: snowarea
   real                           :: dslzsum_i
   real                           :: rdepth
   real                           :: soil_mstpot
   !---------------------------------------------------------------------------------------!

   !----- Time scale for output.  We will use the inverse more often. ---------------------!
   frqsumi = 1.0 / frqsum

   gridloop: do igr=1,ngrids
      cgrid => edgrid_g(igr)

      !------------------------------------------------------------------------------------!
      !   WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! !
      !------------------------------------------------------------------------------------!
      !     Please, don't initialise polygon-level (cgrid) variables outside polyloop.     !
      ! This works in off-line runs, but it causes memory leaks (and crashes) in the       !
      ! coupled runs over the ocean, where cgrid%npolygons can be 0 if one of the sub-     !
      ! domains falls entirely over the ocean.  Thanks!                                    !
      !------------------------------------------------------------------------------------!
      ! cgrid%blah = 0. !<<--- This is a bad way of doing, look inside the loop for the
      !                 !      safe way of initialising the variable.
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         !----- Initialise some integrated variables --------------------------------------!
         area_sum           = 0.0
         laiarea_poly       = 0.0
         poly_avg_soil_hcap = 0.0

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
         cgrid%avg_lai_ebalvars(:,:,ipy) = 0.0
         cgrid%avg_balive          (ipy) = 0.0
         cgrid%avg_bleaf           (ipy) = 0.0
         cgrid%avg_broot           (ipy) = 0.0
         cgrid%avg_bsapwood        (ipy) = 0.0
         cgrid%avg_bdead           (ipy) = 0.0
         cgrid%avg_bstorage        (ipy) = 0.0
         cgrid%avg_bseeds          (ipy) = 0.0
         cgrid%lai                 (ipy) = 0.0
         cgrid%avg_gpp             (ipy) = 0.0
         cgrid%avg_nppleaf         (ipy) = 0.0
         cgrid%avg_nppfroot        (ipy) = 0.0
         cgrid%avg_nppsapwood      (ipy) = 0.0
         cgrid%avg_nppcroot        (ipy) = 0.0
         cgrid%avg_nppseeds        (ipy) = 0.0
         cgrid%avg_nppwood         (ipy) = 0.0
         cgrid%avg_nppdaily        (ipy) = 0.0
         cgrid%avg_leaf_resp       (ipy) = 0.0
         cgrid%avg_root_resp       (ipy) = 0.0
         cgrid%avg_growth_resp     (ipy) = 0.0
         cgrid%avg_storage_resp    (ipy) = 0.0
         cgrid%avg_vleaf_resp      (ipy) = 0.0
         cgrid%avg_plant_resp      (ipy) = 0.0
         cgrid%avg_htroph_resp     (ipy) = 0.0
         cgrid%avg_leaf_drop       (ipy) = 0.0
         cgrid%avg_leaf_maintenance(ipy) = 0.0
         cgrid%avg_root_maintenance(ipy) = 0.0
         cgrid%avg_available_water (ipy) = 0.0
         cgrid%max_leaf_temp       (ipy) = -huge(1.)
         cgrid%min_leaf_temp       (ipy) =  huge(1.)
         cgrid%max_wood_temp       (ipy) = -huge(1.)
         cgrid%min_wood_temp       (ipy) =  huge(1.)

         !----- Inverse of this polygon area (it should be always 1.) ---------------------!
         poly_area_i = 1./sum(cpoly%area)
         siteloop: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)
            
            if (csite%npatches == 0) then
               call fatal_error('No patches in this site, impossible!!!'                   &
                               &,'spatial_averages','edio.f90')
            end if

            !----- Inverse of this site area (it should be always 1.) ---------------------!
            site_area_i=1./sum(csite%area)

            !----- LAI --------------------------------------------------------------------!
            cpoly%lai(isi)  = sum(csite%lai  * csite%area ) * site_area_i
            cpoly%wpa(isi)  = sum(csite%wpa  * csite%area ) * site_area_i
            cpoly%wai(isi)  = sum(csite%wai  * csite%area ) * site_area_i


            !----- Average fast time flux dynamics over sites. ----------------------------!
            cpoly%avg_rshort_gnd(isi)= sum(csite%avg_rshort_gnd* csite%area ) * site_area_i
            cpoly%avg_rlong_gnd(isi) = sum(csite%avg_rlong_gnd * csite%area ) * site_area_i
            cpoly%avg_rlongup(isi)   = sum(csite%avg_rlongup   * csite%area ) * site_area_i
            cpoly%avg_carbon_ac(isi) = sum(csite%avg_carbon_ac * csite%area ) * site_area_i
            cpoly%avg_vapor_lc(isi)  = sum(csite%avg_vapor_lc  * csite%area ) * site_area_i
            cpoly%avg_vapor_wc(isi)  = sum(csite%avg_vapor_wc  * csite%area ) * site_area_i
            cpoly%avg_dew_cg(isi)    = sum(csite%avg_dew_cg    * csite%area ) * site_area_i
            cpoly%avg_vapor_gc(isi)  = sum(csite%avg_vapor_gc  * csite%area ) * site_area_i
            cpoly%avg_wshed_vg(isi)  = sum(csite%avg_wshed_vg  * csite%area ) * site_area_i
            cpoly%avg_vapor_ac(isi)  = sum(csite%avg_vapor_ac  * csite%area ) * site_area_i
            cpoly%avg_transp(isi)    = sum(csite%avg_transp    * csite%area ) * site_area_i
            cpoly%avg_evap(isi)      = sum(csite%avg_evap      * csite%area ) * site_area_i
            cpoly%aux(isi)           = sum(csite%aux           * csite%area ) * site_area_i
            cpoly%avg_drainage(isi)  = sum(csite%avg_drainage  * csite%area ) * site_area_i
            cpoly%avg_runoff(isi)    = sum(csite%avg_runoff    * csite%area ) * site_area_i
            cpoly%avg_drainage_heat(isi) = sum(csite%avg_drainage_heat * csite%area )      &
                                         * site_area_i
            cpoly%avg_runoff_heat(isi)   = sum(csite%avg_runoff_heat   * csite%area )      &
                                         * site_area_i
            cpoly%avg_intercepted(isi)   = sum(csite%avg_intercepted    * csite%area )     &
                                         * site_area_i
            cpoly%avg_qintercepted(isi)  = sum(csite%avg_qintercepted   * csite%area )     &
                                         * site_area_i
            cpoly%avg_throughfall(isi)   = sum(csite%avg_throughfall    * csite%area )     &
                                         * site_area_i
            cpoly%avg_qthroughfall(isi)  = sum(csite%avg_qthroughfall   * csite%area )     &
                                         * site_area_i
            cpoly%avg_sensible_lc(isi)   = sum(csite%avg_sensible_lc    * csite%area )     &
                                         * site_area_i
            cpoly%avg_sensible_wc(isi)   = sum(csite%avg_sensible_wc    * csite%area )     &
                                         * site_area_i
            cpoly%avg_qwshed_vg(isi)     = sum(csite%avg_qwshed_vg      * csite%area )     &
                                         * site_area_i
            cpoly%avg_sensible_gc(isi)   = sum(csite%avg_sensible_gc    * csite%area )     &
                                         * site_area_i
            cpoly%avg_sensible_ac(isi)   = sum(csite%avg_sensible_ac    * csite%area )     &
                                         * site_area_i

            !----- Average albedo values. -------------------------------------------------!
            cpoly%avg_albedo         (isi) = sum(csite%avg_albedo         * csite%area)    &
                                           * site_area_i
            cpoly%avg_albedo_beam    (isi) = sum(csite%avg_albedo_beam    * csite%area)    &
                                           * site_area_i
            cpoly%avg_albedo_diffuse (isi) = sum(csite%avg_albedo_diffuse * csite%area)    &
                                           * site_area_i
            cpoly%avg_rlong_albedo   (isi) = sum(csite%avg_rlong_albedo   * csite%area)    &
                                           * site_area_i

            !----- Extra variables for NACP intercomparision (MCD) ------------------------!
            cpoly%avg_fsc(isi)          = sum(csite%fast_soil_C        * csite%area )      &
                                        * site_area_i
            cpoly%avg_ssc(isi)          = sum(csite%slow_soil_C        * csite%area )      &
                                        * site_area_i
            cpoly%avg_stsc(isi)         = sum(csite%structural_soil_C  * csite%area )      &
                                        * site_area_i
            cpoly%avg_fsn(isi)          = sum(csite%fast_soil_N        * csite%area )      &
                                        * site_area_i
            cpoly%avg_msn(isi)          = sum(csite%mineralized_soil_N * csite%area )      &
                                        * site_area_i

            !----- Available water. -------------------------------------------------------!
            cpoly%avg_available_water(isi) = sum(csite%avg_available_water * csite%area)   &
                                           * site_area_i
            !------------------------------------------------------------------------------!
            !------------------------------------------------------------------------------!
            !     Find average soil properties.  The average soil temperature and liquid   !
            ! fraction is not directly computed, since the total mass and therefore the    !
            ! total heat capacity (dry + water/ice) is different from each patch.          !
            !------------------------------------------------------------------------------!
            cpoly%avg_soil_wetness(isi)     = 0.0
            cpoly%avg_sensible_gg(:,isi) = matmul(csite%avg_sensible_gg,csite%area)        &
                                         * site_area_i
            cpoly%avg_smoist_gg(:,isi)   = matmul(csite%avg_smoist_gg  ,csite%area)        &
                                         * site_area_i
            cpoly%avg_transloss(:,isi)   = matmul(csite%avg_transloss  ,csite%area)        &
                                         * site_area_i
            cpoly%aux_s(:,isi)           = matmul(csite%aux_s          ,csite%area)        &
                                         * site_area_i
            cpoly%avg_soil_energy(:,isi) = matmul(csite%soil_energy    ,csite%area)        &
                                         * site_area_i
            cpoly%avg_soil_water(:,isi)  = matmul(csite%soil_water     ,csite%area)        &
                                         * site_area_i

            !------------------------------------------------------------------------------!
            !     Initialise soil moisture potential.                                      !
            !------------------------------------------------------------------------------!
            cpoly%avg_soil_mstpot(:,isi) = 0.0

            do k=cpoly%lsl(isi),nzg
               !---------------------------------------------------------------------------!
               !    Find the mean heat capacity. This shouldn't matter in the current     !
               ! version as all patches within the same site have the same soil texture.   !
               ! Since heat capacity is given in J/m3/K, it's safe to average it even if   !
               ! soil texture changed.                                                     !
               !---------------------------------------------------------------------------!
               site_avg_soil_hcap(k) = 0.
               dslzsum_i = 1./ sum(dslz(cpoly%lsl(isi):nzg))

               do ipa=1,csite%npatches
                  nsoil = cpoly%ntext_soil(k,isi)
                  site_avg_soil_hcap(k) = site_avg_soil_hcap(k)                            &
                                        + soil(cpoly%ntext_soil(k,isi))%slcpd              &
                                        * csite%area(ipa) * site_area_i
                  !----- Integrate soil wetness. ------------------------------------------!
                  cpoly%avg_soil_wetness(isi) = cpoly%avg_soil_wetness(isi)                &
                       + ((csite%soil_water(k,ipa) - soil(nsoil)%soilwp))                  &
                       / (soil(nsoil)%slmsts - soil(nsoil)%soilwp)                         &
                       * dslz(k) * dslzsum_i * csite%area(ipa) * site_area_i

                  !------ Integrate the average soil moisture potential. ------------------!
                  soil_mstpot = soil(nsoil)%slpots                                         &
                              * (soil(nsoil)%slmsts / csite%soil_water(k,ipa))             &
                              ** soil(nsoil)%slbs
                  cpoly%avg_soil_mstpot(k,isi) = cpoly%avg_soil_mstpot(k,isi)              &
                                               + soil_mstpot * csite%area(ipa)             &
                                               * site_area_i
               end do
               !----- Also integrate the polygon-level average. ---------------------------!
               poly_avg_soil_hcap(k) = poly_avg_soil_hcap(k)                               &
                                     + site_avg_soil_hcap(k) * cpoly%area(isi)*poly_area_i

               !----- Finding the average temperature and liquid fraction. ----------------!
               call qwtk(cpoly%avg_soil_energy(k,isi),cpoly%avg_soil_water(k,isi)*wdns     &
                        ,site_avg_soil_hcap(k),cpoly%avg_soil_temp(k,isi)                  &
                        ,cpoly%avg_soil_fracliq(k,isi))
            end do
            
            !------------------------------------------------------------------------------! 
            !     For layers beneath the lowest soil level, assign a default soil          !
            ! potential and soil moisture consistent with the boundary condition.          !
            !------------------------------------------------------------------------------!
            select case (isoilbc)
            case (0,1)
               !----- Copy the potential from bottommost layer ----------------------------!
               do k=1,cpoly%lsl(isi)-1
                  cpoly%avg_soil_mstpot(k,isi) = cpoly%avg_soil_mstpot(cpoly%lsl(isi),isi)
               end do

            case (2)
               !----- Super-drainage, assume dry air soil beneath. ------------------------!
               do k=1,cpoly%lsl(isi)-1
                  nsoil = cpoly%ntext_soil(k,isi)
                  if (nsoil /= 13) then
                     cpoly%avg_soil_mstpot(k,isi) = soil(nsoil)%slpots                     &
                                                  * ( soil(nsoil)%slmsts                   &
                                                    / soil(nsoil)%soilcp )                 &
                                                  ** soil(nsoil)%slbs
                  end if
               end do
            case (3)
               !----- Aquifer, assume saturated soil moisture. ----------------------------!
               do k=1,cpoly%lsl(isi)-1
                  nsoil = cpoly%ntext_soil(k,isi)
                  cpoly%avg_soil_mstpot(k,isi) = soil(nsoil)%slpots
               end do
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Temporary water/snow layer must be averaged differently, since each patch !
            ! may or may not have the layer, and the number of layers may vary from patch  !
            ! to patch.  Here we average the integrated energy, mass and depth, and        !
            ! compute the average temperature and liquid fraction from the average.        !
            !------------------------------------------------------------------------------!
            cpoly%avg_sfcw_depth(isi)   = 0.0
            cpoly%avg_sfcw_energy(isi)  = 0.0
            cpoly%avg_sfcw_mass(isi)    = 0.0

            do ipa=1,csite%npatches
               ksn = csite%nlev_sfcwater(ipa)
               !----- Check whether mass is present.  If so, add this patch. --------------!
               if (ksn > 0) then
                  do k=1,ksn
                     cpoly%avg_sfcw_depth(isi)  = cpoly%avg_sfcw_depth(isi)                &
                                                + csite%sfcwater_depth(k,ipa)              &
                                                * csite%area(ipa) * site_area_i
                     cpoly%avg_sfcw_mass(isi)   = cpoly%avg_sfcw_mass(isi)                 &
                                                + csite%sfcwater_mass(k,ipa)               &
                                                * csite%area(ipa) * site_area_i
                     !---------------------------------------------------------------------!
                     !     Internal energy.  For averaging we use the extensive one (J/m2) !
                     ! we will switch back to J/kg after the mean mass is found.           !
                     !---------------------------------------------------------------------!
                     cpoly%avg_sfcw_energy(isi) = cpoly%avg_sfcw_energy(isi)               &
                                                + csite%sfcwater_energy(k,ipa)             &
                                                * csite%sfcwater_mass(k,ipa)               &
                                                * csite%area(ipa) * site_area_i
                  end do
               end if
            end do

            !------------------------------------------------------------------------------!
            !   If the site had some layer, convert the mean energy to J/kg, then find     !
            ! the mean temperature and liquid fraction.  Otherwise, make them with zero/   !
            ! default values.                                                              !
            !------------------------------------------------------------------------------!
            if (cpoly%avg_sfcw_mass(isi) > tiny_sfcwater_mass) then
               cpoly%avg_sfcw_energy(isi) = cpoly%avg_sfcw_energy(isi)                     &
                                          / cpoly%avg_sfcw_mass(isi)
               call qtk(cpoly%avg_sfcw_energy(isi),cpoly%avg_sfcw_tempk(isi)               &
                       ,cpoly%avg_sfcw_fracliq(isi))
            else
               cpoly%avg_sfcw_mass(isi)    = 0.
               cpoly%avg_sfcw_depth(isi)   = 0.
               cpoly%avg_sfcw_energy(isi)  = 0.
               cpoly%avg_sfcw_tempk(isi)   = cpoly%avg_soil_temp(nzg,isi)
               cpoly%avg_sfcw_fracliq(isi) = cpoly%avg_soil_fracliq(nzg,isi)
            end if
            !------------------------------------------------------------------------------!

            !----- Average over patches. --------------------------------------------------!
            longpatchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               !----- Zero the rootfraction diagnostic. -----------------------------------!
               csite%rootdense(:,ipa) = 0.

               !---------------------------------------------------------------------------!
               !     Adding cohort "extensive" variables. Those that are not must be       !
               ! scaled by nplant. Just make sure that we have at least one cohort.        !
               !---------------------------------------------------------------------------!
               if (cpatch%ncohorts > 0) then
                  lai_patch = sum(cpatch%lai, cpatch%leaf_resolvable)

                  !------------------------------------------------------------------------!
                  !     Leaf water and energy properties.                                  !
                  !------------------------------------------------------------------------!
                  csite%avg_leaf_energy(ipa)  = sum(cpatch%leaf_energy)
                  csite%avg_leaf_water (ipa)  = sum(cpatch%leaf_water)
                  csite%avg_leaf_hcap  (ipa)  = sum(cpatch%leaf_hcap)
                  !----- Check whether there is any heat storage. -------------------------!
                  if (csite%avg_leaf_hcap(ipa) > 0.) then
                     !----- Yes, use the default thermodynamics. --------------------------!
                     call qwtk(csite%avg_leaf_energy(ipa),csite%avg_leaf_water(ipa)        &
                              ,csite%avg_leaf_hcap(ipa),csite%avg_leaf_temp(ipa)           &
                              ,csite%avg_leaf_fliq(ipa))
                  else
                     !----- No, copy the canopy air properties. ---------------------------!
                     csite%avg_leaf_temp(ipa) = csite%can_temp(ipa)
                     if (csite%can_temp(ipa) > t00) then
                        csite%avg_leaf_fliq(ipa) = 1.0
                     elseif (csite%can_temp(ipa) == t00) then
                        csite%avg_leaf_fliq(ipa) = 0.5 
                     else
                        csite%avg_leaf_fliq(ipa) = 0.0 
                     end if
                  end if
                  !------------------------------------------------------------------------!
                  !     Stem water and energy properties.
                  !------------------------------------------------------------------------!
                  csite%avg_wood_energy(ipa)  = sum(cpatch%wood_energy)
                  csite%avg_wood_water (ipa)  = sum(cpatch%wood_water)
                  csite%avg_wood_hcap  (ipa)  = sum(cpatch%wood_hcap)
                  !----- Check whether there is any heat storage. -------------------------!
                  if (csite%avg_wood_hcap(ipa) > 0.) then
                     !----- Yes, use the default thermodynamics. --------------------------!
                     call qwtk(csite%avg_wood_energy(ipa),csite%avg_wood_water(ipa)        &
                              ,csite%avg_wood_hcap(ipa),csite%avg_wood_temp(ipa)           &
                              ,csite%avg_wood_fliq(ipa))
                  else
                     !----- No, copy the canopy air properties. ---------------------------!
                     csite%avg_wood_temp(ipa) = csite%can_temp(ipa)
                     if (csite%can_temp(ipa) > t00) then
                        csite%avg_wood_fliq(ipa) = 1.0
                     elseif (csite%can_temp(ipa) == t00) then
                        csite%avg_wood_fliq(ipa) = 0.5 
                     else
                        csite%avg_wood_fliq(ipa) = 0.0 
                     end if
                  end if
                  !------------------------------------------------------------------------!

                  cgrid%avg_gpp(ipy)          = cgrid%avg_gpp(ipy)                         &
                                              + sum(cpatch%mean_gpp)                       &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i

                  cgrid%avg_nppleaf(ipy)      = cgrid%avg_nppleaf(ipy)                     &
                                              + sum(cpatch%today_nppleaf                   &
                                              * cpatch%nplant)                             &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i
                                              
                  cgrid%avg_nppfroot(ipy)     = cgrid%avg_nppfroot(ipy)                    &
                                              + sum(cpatch%today_nppfroot                  &
                                              * cpatch%nplant)                             &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i
                                              
                  cgrid%avg_nppsapwood(ipy)   = cgrid%avg_nppsapwood(ipy)                  &
                                              + sum(cpatch%today_nppsapwood                &
                                              * cpatch%nplant)                             &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i
                                              
                  cgrid%avg_nppcroot(ipy)     = cgrid%avg_nppcroot(ipy)                    &
                                              + sum(cpatch%today_nppcroot                  &
                                              * cpatch%nplant)                             &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i
                                              
                  cgrid%avg_nppseeds(ipy)     = cgrid%avg_nppseeds(ipy)                    &
                                              + sum(cpatch%today_nppseeds                  &
                                              * cpatch%nplant)                             &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i
                                              
                  cgrid%avg_nppwood(ipy)      = cgrid%avg_nppwood(ipy)                     &
                                              + sum(cpatch%today_nppwood                   &
                                              * cpatch%nplant)                             &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i
                                              
                  cgrid%avg_nppdaily(ipy)     = cgrid%avg_nppdaily(ipy)                    &
                                              + sum(cpatch%today_nppdaily                  &
                                              * cpatch%nplant)                             &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i
                                              
                  cgrid%avg_leaf_resp(ipy)    = cgrid%avg_leaf_resp(ipy)                   &
                                              + sum(cpatch%mean_leaf_resp)                 &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i

                  cgrid%avg_root_resp(ipy)    = cgrid%avg_root_resp(ipy)                   &
                                              + sum(cpatch%mean_root_resp)                 &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i

                  cgrid%avg_growth_resp(ipy)  = cgrid%avg_growth_resp(ipy)                 &
                                              + sum(cpatch%mean_growth_resp)               &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i

                  cgrid%avg_storage_resp(ipy) = cgrid%avg_storage_resp(ipy)                &
                                              + sum(cpatch%mean_storage_resp)              &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i

                  cgrid%avg_vleaf_resp(ipy)   = cgrid%avg_vleaf_resp(ipy)                  &
                                              + sum(cpatch%mean_vleaf_resp)                &
                                              * csite%area(ipa)*cpoly%area(isi)            &
                                              * site_area_i * poly_area_i
                  !------------------------------------------------------------------------!

                  cgrid%avg_balive(ipy)     = cgrid%avg_balive(ipy)                        &
                                            + sum(cpatch%balive*cpatch%nplant)             &
                                            * csite%area(ipa)*cpoly%area(isi)              &
                                            * site_area_i * poly_area_i

                  cgrid%avg_bleaf(ipy)      = cgrid%avg_bleaf(ipy)                         &
                                            + sum(cpatch%bleaf*cpatch%nplant)              &
                                            * csite%area(ipa)*cpoly%area(isi)              &
                                            * site_area_i * poly_area_i

                  cgrid%avg_broot(ipy)      = cgrid%avg_broot(ipy)                         &
                                            + sum(cpatch%broot*cpatch%nplant)              &
                                            * csite%area(ipa)*cpoly%area(isi)              &
                                            * site_area_i * poly_area_i

                  cgrid%avg_bsapwood(ipy)   = cgrid%avg_bsapwood(ipy)                      &
                                            + sum(cpatch%bsapwood*cpatch%nplant)           &
                                            * csite%area(ipa)*cpoly%area(isi)              &
                                            * site_area_i * poly_area_i

                  cgrid%avg_bdead(ipy)      = cgrid%avg_bdead(ipy)                         &
                                            + sum(cpatch%bdead*cpatch%nplant)              &
                                            * csite%area(ipa)*cpoly%area(isi)              &
                                            * site_area_i * poly_area_i

                  cgrid%avg_bstorage(ipy)   = cgrid%avg_bstorage(ipy)                      &
                                            + sum(cpatch%bstorage*cpatch%nplant)           &
                                            * csite%area(ipa)*cpoly%area(isi)              &
                                            * site_area_i * poly_area_i

                  cgrid%avg_bseeds(ipy)     = cgrid%avg_bseeds(ipy)                        &
                                            + sum(cpatch%bseeds*cpatch%nplant)             &
                                            * csite%area(ipa)*cpoly%area(isi)              &
                                            * site_area_i * poly_area_i

                  !------------------------------------------------------------------------!
                  !      Leaf drop due to phenology and maintenance costs... I left it     !
                  ! in kgC/m2/yr...                                                        !
                  !------------------------------------------------------------------------!
                  cgrid%avg_leaf_drop(ipy)         = cgrid%avg_leaf_drop(ipy)              &
                                                   + sum( cpatch%leaf_drop                 &
                                                        * cpatch%nplant)                   &
                                                   * csite%area(ipa)*cpoly%area(isi)       &
                                                   * site_area_i * poly_area_i
                  cgrid%avg_leaf_maintenance(ipy)  = cgrid%avg_leaf_maintenance(ipy)       &
                                                   + sum( cpatch%leaf_maintenance          &
                                                        * cpatch%nplant)                   &
                                                   * csite%area(ipa)*cpoly%area(isi)       &
                                                   * site_area_i * poly_area_i
                  cgrid%avg_root_maintenance(ipy)  = cgrid%avg_root_maintenance(ipy)       &
                                                   + sum( cpatch%root_maintenance          &
                                                        * cpatch%nplant)                   &
                                                   * csite%area(ipa)*cpoly%area(isi)       &
                                                   * site_area_i * poly_area_i

                  !----- Check the extremes and update if necessary. ----------------------!
                  if (maxval(cpatch%leaf_temp) > cgrid%max_leaf_temp(ipy)) then
                     cgrid%max_leaf_temp(ipy) = maxval(cpatch%leaf_temp)
                  end if
                  if (minval(cpatch%leaf_temp) < cgrid%min_leaf_temp(ipy)) then
                     cgrid%min_leaf_temp(ipy) = minval(cpatch%leaf_temp)
                  end if
               else
                  lai_patch                   = 0.
                  csite%avg_leaf_energy(ipa)  = 0.
                  csite%avg_leaf_water(ipa)   = 0.
                  csite%avg_leaf_hcap(ipa)    = 0.
                  csite%avg_leaf_temp(ipa)    = csite%can_temp(ipa)
                  csite%avg_wood_energy(ipa)  = 0.
                  csite%avg_wood_water(ipa)   = 0.
                  csite%avg_wood_hcap(ipa)    = 0.
                  csite%avg_wood_temp(ipa)    = csite%can_temp(ipa)
                  if (csite%can_temp(ipa) > t00) then
                     csite%avg_leaf_fliq(ipa) = 1.
                     csite%avg_wood_fliq(ipa) = 1.
                  elseif (csite%can_temp(ipa) == t00) then
                     csite%avg_leaf_fliq(ipa) = 0.5
                     csite%avg_wood_fliq(ipa) = 0.5
                  else
                     csite%avg_leaf_fliq(ipa) = 0.0
                     csite%avg_wood_fliq(ipa) = 0.0
                  end if
               end if

               if (lai_patch > 0.) then
                  csite%laiarea(ipa) = csite%area(ipa)
               else
                  csite%laiarea(ipa) = 0.
               end if

               cgrid%avg_htroph_resp(ipy) = cgrid%avg_htroph_resp(ipy)                     &
                                          + csite%mean_rh(ipa)                             &
                                          * csite%area(ipa)*cpoly%area(isi)                &
                                          * site_area_i * poly_area_i

               !----- Not sure what these variables do. -----------------------------------!
               lai_index = min(3,max(1, floor(csite%lai(ipa)/2.0) + 1)  )
               area_sum(lai_index) = area_sum(lai_index) + csite%area(ipa)

               !----- Net radiation. ------------------------------------------------------!
               cgrid%avg_lai_ebalvars(lai_index,1,ipy) =                                   &
                      cgrid%avg_lai_ebalvars(lai_index,1,ipy)                              &
                    + (csite%avg_rshort_gnd(ipa) + csite%avg_rlong_gnd(ipa))               &
                    * csite%area(ipa)*cpoly%area(isi)                                      &
                    * site_area_i * poly_area_i

               !----- Latent heat flux. ---------------------------------------------------!
               cgrid%avg_lai_ebalvars(lai_index,2,ipy) =                                   &
                      cgrid%avg_lai_ebalvars(lai_index,2,ipy)                              &
                    + csite%avg_vapor_ac(ipa)                                              &
                    * csite%area(ipa)*cpoly%area(isi)                                      &
                    * site_area_i * poly_area_i

               !----- Sensible heat flux. -------------------------------------------------!
               cgrid%avg_lai_ebalvars(lai_index,3,ipy) =                                   &
                      cgrid%avg_lai_ebalvars(lai_index,3,ipy)                              &
                    + csite%avg_sensible_ac(ipa)                                           &
                    * csite%area(ipa)*cpoly%area(isi)                                      &
                    * site_area_i * poly_area_i

               !----- Canopy temperature --------------------------------------------------!
               cgrid%avg_lai_ebalvars(lai_index,4,ipy) =                                   &
                      cgrid%avg_lai_ebalvars(lai_index,4,ipy)                              &
                    + csite%can_temp(ipa)                                                  &
                    * csite%area(ipa)*cpoly%area(isi)                                      &
                    * site_area_i * poly_area_i


               !------ Updating maximum and minimum soil temperature. ---------------------!
               do k=1,nzg
                  if (csite%soil_tempk(k,ipa) > cgrid%max_soil_temp(ipy)) then
                     cgrid%max_soil_temp(ipy) = csite%soil_tempk(k,ipa)
                  end if
                 
                  if (csite%soil_tempk(k,ipa) < cgrid%min_soil_temp(ipy)) then
                     cgrid%min_soil_temp(ipy) = csite%soil_tempk(k,ipa)
                  end if
               end do


               cohortloop: do ico=1,cpatch%ncohorts
                  cpatch%old_stoma_vector(1,ico) = real(cpatch%old_stoma_data(ico)%recalc)
                  cpatch%old_stoma_vector(2,ico) = cpatch%old_stoma_data(ico)%T_L
                  cpatch%old_stoma_vector(3,ico) = cpatch%old_stoma_data(ico)%e_A
                  cpatch%old_stoma_vector(4,ico) = cpatch%old_stoma_data(ico)%PAR
                  cpatch%old_stoma_vector(5,ico) = cpatch%old_stoma_data(ico)%rb_factor
                  cpatch%old_stoma_vector(6,ico) = cpatch%old_stoma_data(ico)%prss
                  cpatch%old_stoma_vector(7,ico) = cpatch%old_stoma_data(ico)%phenology_factor
                  cpatch%old_stoma_vector(8,ico) = cpatch%old_stoma_data(ico)%gsw_open
                  cpatch%old_stoma_vector(9,ico) = real(cpatch%old_stoma_data(ico)%ilimit)
                  
                  cpatch%old_stoma_vector(10,ico) = cpatch%old_stoma_data(ico)%T_L_residual
                  cpatch%old_stoma_vector(11,ico) = cpatch%old_stoma_data(ico)%e_a_residual
                  cpatch%old_stoma_vector(12,ico) = cpatch%old_stoma_data(ico)%par_residual
                  cpatch%old_stoma_vector(13,ico) = cpatch%old_stoma_data(ico)%rb_residual
                  cpatch%old_stoma_vector(14,ico) = cpatch%old_stoma_data(ico)%prss_residual
                  cpatch%old_stoma_vector(15,ico) = cpatch%old_stoma_data(ico)%leaf_residual
                  cpatch%old_stoma_vector(16,ico) = cpatch%old_stoma_data(ico)%gsw_residual


                  !------------------------------------------------------------------------!
                  !    Rooting fraction: step 1, find root biomass per cubic meter         !
                  !    broot*nplant/rooting_depth   [kg/plant]*[plant/m2]/[m]              !
                  !------------------------------------------------------------------------!
                  rdepth = sum(dslz(cpatch%krdepth(ico):nzg))
                  do k=cpatch%krdepth(ico),nzg
                     csite%rootdense(k,ipa) = csite%rootdense(k,ipa)                       &
                                            + cpatch%broot(ico)*cpatch%nplant(ico)/rdepth
                  end do
                  !------------------------------------------------------------------------!


               end do cohortloop
                  
               pftloop: do ipft = 1,n_pft
                  csite%old_stoma_vector_max(1,ipft,ipa) = real(csite%old_stoma_data_max(ipft,ipa)%recalc)
                  csite%old_stoma_vector_max(2,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%T_L
                  csite%old_stoma_vector_max(3,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%e_A
                  csite%old_stoma_vector_max(4,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%PAR
                  csite%old_stoma_vector_max(5,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%rb_factor
                  csite%old_stoma_vector_max(6,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%prss
                  csite%old_stoma_vector_max(7,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%phenology_factor
                  csite%old_stoma_vector_max(8,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%gsw_open
                  csite%old_stoma_vector_max(9,ipft,ipa) = real(csite%old_stoma_data_max(ipft,ipa)%ilimit)
                  
                  csite%old_stoma_vector_max(10,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%T_L_residual
                  csite%old_stoma_vector_max(11,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%e_a_residual
                  csite%old_stoma_vector_max(12,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%par_residual
                  csite%old_stoma_vector_max(13,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%rb_residual
                  csite%old_stoma_vector_max(14,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%prss_residual
                  csite%old_stoma_vector_max(15,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%leaf_residual
                  csite%old_stoma_vector_max(16,ipft,ipa) = csite%old_stoma_data_max(ipft,ipa)%gsw_residual
               end do pftloop
            
            end do longpatchloop

            laiarea_site  = sum(csite%laiarea)
            laiarea_poly  = laiarea_poly + laiarea_site
            if (laiarea_site == 0.0) then
               csite%laiarea = 0.0
            else
               csite%laiarea = csite%laiarea / laiarea_site
            end if


            ! Take an area weighted average of the root density to get site level fraction
            !------------------------------------------------------------------------------!
            cpoly%avg_soil_rootfrac(:,isi) = matmul(csite%rootdense,csite%area)        &
                 * site_area_i

            
            !------------------------------------------------------------------------------!
            !    Site average of canopy thermodynamic state.  We average the variables     !
            ! that are insensitive to changes in pressure (potential temperature, specific !
            ! humidity, CO2 mixing ratio and pressure itself) by averaging over sites.     !
            ! The other variables are found using the prescribed diagnostic equations.     !
            !------------------------------------------------------------------------------!
            cpoly%avg_can_depth   (isi) = sum(csite%can_depth  * csite%area) * site_area_i
            cpoly%avg_can_theta   (isi) = sum(csite%can_theta  * csite%area) * site_area_i
            cpoly%avg_can_theiv   (isi) = sum(csite%can_theiv  * csite%area) * site_area_i
            cpoly%avg_can_shv     (isi) = sum(csite%can_shv    * csite%area) * site_area_i
            cpoly%avg_can_co2     (isi) = sum(csite%can_co2    * csite%area) * site_area_i
            cpoly%avg_can_prss    (isi) = sum(csite%can_prss   * csite%area) * site_area_i
            cpoly%avg_can_temp    (isi) = cpoly%avg_can_theta(isi)                         &
                                        * (p00i * cpoly%avg_can_prss(isi)) ** rocp
            cpoly%avg_can_rhos    (isi) = idealdenssh(cpoly%avg_can_prss(isi)              &
                                                     ,cpoly%avg_can_temp(isi)              &
                                                     ,cpoly%avg_can_shv (isi) )
            !------------------------------------------------------------------------------!
            !   Site average of leaf and stem properties.  Again, we average "extensive"   !
            ! properties and find the average temperature based on the average leaf and    !
            ! stem internal energy mass, and heat capacity.                                !
            !------------------------------------------------------------------------------!
            cpoly%avg_leaf_energy(isi) = sum(csite%avg_leaf_energy * csite%area)           &
                                       * site_area_i
            cpoly%avg_leaf_water (isi) = sum(csite%avg_leaf_water  * csite%area)           &
                                       * site_area_i
            cpoly%avg_leaf_hcap  (isi) = sum(csite%avg_leaf_hcap   * csite%area)           &
                                       * site_area_i
            cpoly%avg_wood_energy(isi) = sum(csite%avg_wood_energy * csite%area)           &
                                       * site_area_i
            cpoly%avg_wood_water (isi) = sum(csite%avg_wood_water  * csite%area)           &
                                       * site_area_i
            cpoly%avg_wood_hcap  (isi) = sum(csite%avg_wood_hcap   * csite%area)           &
                                       * site_area_i
            !------------------------------------------------------------------------------!
            !    Unless this is a bare site or it is absolute leafless, there should be    !
            ! some heat capacity, then compute the average leaf temperature.  Otherwise,   !
            ! assign mean canopy temperature.                                              !
            !------------------------------------------------------------------------------!
            if (cpoly%avg_leaf_hcap(isi) > 0.) then
               call qwtk(cpoly%avg_leaf_energy(isi),cpoly%avg_leaf_water(isi)                &
                        ,cpoly%avg_leaf_hcap(isi),cpoly%avg_leaf_temp(isi)                   &
                        ,cpoly%avg_leaf_fliq(isi))
            else
               cpoly%avg_leaf_temp(isi) = cpoly%avg_can_temp(isi)
               if (cpoly%avg_can_temp(isi) > t00) then
                  cpoly%avg_leaf_fliq(isi) = 1.0
               elseif (cpoly%avg_can_temp(isi) == t00) then
                  cpoly%avg_leaf_fliq(isi) = 0.5
               else
                  cpoly%avg_leaf_fliq(isi) = 0.0
               end if
            end if
            !------------------------------------------------------------------------------!
            !     Unless this is a bare site, or the user doesn't want to solve branches,  !
            ! there should be some heat capacity, in which case the average stem temper-   !
            ! ature.  Otherwise, assign mean canopy temperature.                           !
            !------------------------------------------------------------------------------!
            if (cpoly%avg_wood_hcap(isi) > 0.) then
               call qwtk(cpoly%avg_wood_energy(isi),cpoly%avg_wood_water(isi)                &
                        ,cpoly%avg_wood_hcap(isi),cpoly%avg_wood_temp(isi)                   &
                        ,cpoly%avg_wood_fliq(isi))
            else
               cpoly%avg_wood_temp(isi) = cpoly%avg_can_temp(isi)
               if (cpoly%avg_can_temp(isi) > t00) then
                  cpoly%avg_wood_fliq(isi) = 1.0
               elseif (cpoly%avg_can_temp(isi) == t00) then
                  cpoly%avg_wood_fliq(isi) = 0.5
               else
                  cpoly%avg_wood_fliq(isi) = 0.0
               end if
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !    Compute the average amongst all surfaces (soil, temporary surface water,  !
            ! and vegetation, the last two only if they really exist).  All energy terms   !
            ! are converted to J/m2, all water terms to kg/m2, and the heat capacities of  !
            ! everything that is not water is in J/m2/K.                                   !
            !------------------------------------------------------------------------------!
            skin_energy = cpoly%avg_leaf_energy(isi)                                       &
                        + cpoly%avg_wood_energy(isi)                                       &
                        + cpoly%avg_sfcw_energy(isi) * cpoly%avg_sfcw_mass(isi)            &
                        + cpoly%avg_soil_energy(nzg,isi) * dslz(nzg)
            skin_water  = cpoly%avg_leaf_water(isi)                                        &
                        + cpoly%avg_leaf_water(isi)                                        &
                        + cpoly%avg_sfcw_mass(isi)                                         &
                        + cpoly%avg_soil_water(nzg,isi) * dslz(nzg) * wdns
            skin_hcap   = cpoly%avg_leaf_hcap(isi)                                         &
                        + cpoly%avg_wood_hcap(isi)                                         &
                        + site_avg_soil_hcap(nzg) * dslz(nzg)
            call qwtk(skin_energy,skin_water,skin_hcap,cpoly%avg_skin_temp(isi),skin_fliq)
            !------------------------------------------------------------------------------!
         end do siteloop
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Find plant respiration using the other averaged quantities.                  !
         !---------------------------------------------------------------------------------!
         cgrid%avg_plant_resp(ipy)  = cgrid%avg_leaf_resp     (ipy)                        &
                                    + cgrid%avg_root_resp     (ipy)                        &
                                    + cgrid%avg_growth_resp   (ipy)                        &
                                    + cgrid%avg_storage_resp  (ipy)                        &
                                    + cgrid%avg_vleaf_resp    (ipy)
         !---------------------------------------------------------------------------------!



         !----- Normalize the lai specific quantities -------------------------------------!
         if (area_sum(1)>0.) then
            cgrid%avg_lai_ebalvars(1,1,ipy) = cgrid%avg_lai_ebalvars(1,1,ipy)/area_sum(1)
            cgrid%avg_lai_ebalvars(1,2,ipy) = cgrid%avg_lai_ebalvars(1,2,ipy)/area_sum(1)
            cgrid%avg_lai_ebalvars(1,3,ipy) = cgrid%avg_lai_ebalvars(1,3,ipy)/area_sum(1)
            cgrid%avg_lai_ebalvars(1,4,ipy) = cgrid%avg_lai_ebalvars(1,4,ipy)/area_sum(1)
         else
            cgrid%avg_lai_ebalvars(1,:,ipy) = -9999.0
         end if
         if (area_sum(2)>0.) then
            cgrid%avg_lai_ebalvars(2,1,ipy) = cgrid%avg_lai_ebalvars(2,1,ipy)/area_sum(2)
            cgrid%avg_lai_ebalvars(2,2,ipy) = cgrid%avg_lai_ebalvars(2,2,ipy)/area_sum(2)
            cgrid%avg_lai_ebalvars(2,3,ipy) = cgrid%avg_lai_ebalvars(2,3,ipy)/area_sum(2)
            cgrid%avg_lai_ebalvars(2,4,ipy) = cgrid%avg_lai_ebalvars(2,4,ipy)/area_sum(2)
         else
            cgrid%avg_lai_ebalvars(2,:,ipy) = -9999.0
         end if
         if (area_sum(3)>0.) then
            cgrid%avg_lai_ebalvars(3,1,ipy) = cgrid%avg_lai_ebalvars(3,1,ipy)/area_sum(3)
            cgrid%avg_lai_ebalvars(3,2,ipy) = cgrid%avg_lai_ebalvars(3,2,ipy)/area_sum(3)
            cgrid%avg_lai_ebalvars(3,3,ipy) = cgrid%avg_lai_ebalvars(3,3,ipy)/area_sum(3)
            cgrid%avg_lai_ebalvars(3,4,ipy) = cgrid%avg_lai_ebalvars(3,4,ipy)/area_sum(3)
         else
            cgrid%avg_lai_ebalvars(3,:,ipy) = -9999.0
         end if

         !----- Finding the polygon mean LAI ----------------------------------------------!
         cgrid%lai(ipy)  = sum(cpoly%lai  * cpoly%area ) * poly_area_i
         cgrid%wpa(ipy)  = sum(cpoly%wpa  * cpoly%area ) * poly_area_i
         cgrid%wai(ipy)  = sum(cpoly%wai  * cpoly%area ) * poly_area_i
         !----- Average fast time flux dynamics over polygons. ----------------------------!
         cgrid%avg_rshort_gnd(ipy)   = sum(cpoly%avg_rshort_gnd   *cpoly%area)*poly_area_i
         cgrid%avg_rlong_gnd(ipy)    = sum(cpoly%avg_rlong_gnd    *cpoly%area)*poly_area_i
         cgrid%avg_rlongup  (ipy)    = sum(cpoly%avg_rlongup      *cpoly%area)*poly_area_i
         cgrid%avg_carbon_ac(ipy)    = sum(cpoly%avg_carbon_ac    *cpoly%area)*poly_area_i
         cgrid%avg_vapor_lc(ipy)     = sum(cpoly%avg_vapor_lc     *cpoly%area)*poly_area_i
         cgrid%avg_vapor_wc(ipy)     = sum(cpoly%avg_vapor_wc     *cpoly%area)*poly_area_i
         cgrid%avg_dew_cg(ipy)       = sum(cpoly%avg_dew_cg       *cpoly%area)*poly_area_i
         cgrid%avg_vapor_gc(ipy)     = sum(cpoly%avg_vapor_gc     *cpoly%area)*poly_area_i
         cgrid%avg_wshed_vg(ipy)     = sum(cpoly%avg_wshed_vg     *cpoly%area)*poly_area_i
         cgrid%avg_intercepted(ipy)  = sum(cpoly%avg_intercepted  *cpoly%area)*poly_area_i
         cgrid%avg_throughfall(ipy)  = sum(cpoly%avg_throughfall  *cpoly%area)*poly_area_i
         cgrid%avg_fsc(ipy)          = sum(cpoly%avg_fsc          *cpoly%area)*poly_area_i
         cgrid%avg_stsc(ipy)         = sum(cpoly%avg_stsc         *cpoly%area)*poly_area_i
         cgrid%avg_ssc(ipy)          = sum(cpoly%avg_ssc          *cpoly%area)*poly_area_i
         cgrid%avg_fsn(ipy)          = sum(cpoly%avg_fsn          *cpoly%area)*poly_area_i
         cgrid%avg_msn(ipy)          = sum(cpoly%avg_msn          *cpoly%area)*poly_area_i
         cgrid%avg_runoff_heat(ipy)  = sum(cpoly%avg_runoff_heat  *cpoly%area)*poly_area_i
         cgrid%avg_runoff(ipy)       = sum(cpoly%avg_runoff       *cpoly%area)*poly_area_i
         cgrid%avg_drainage(ipy)     = sum(cpoly%avg_drainage     *cpoly%area)*poly_area_i
         cgrid%avg_drainage_heat(ipy)= sum(cpoly%avg_drainage_heat*cpoly%area)*poly_area_i
         cgrid%avg_vapor_ac(ipy)     = sum(cpoly%avg_vapor_ac     *cpoly%area)*poly_area_i
         cgrid%avg_transp(ipy)       = sum(cpoly%avg_transp       *cpoly%area)*poly_area_i
         cgrid%avg_evap(ipy)         = sum(cpoly%avg_evap         *cpoly%area)*poly_area_i
         cgrid%aux(ipy)              = sum(cpoly%aux              *cpoly%area)*poly_area_i
         cgrid%avg_sensible_lc(ipy)  = sum(cpoly%avg_sensible_lc  *cpoly%area)*poly_area_i
         cgrid%avg_sensible_wc(ipy)  = sum(cpoly%avg_sensible_wc  *cpoly%area)*poly_area_i
         cgrid%avg_qwshed_vg(ipy)    = sum(cpoly%avg_qwshed_vg    *cpoly%area)*poly_area_i
         cgrid%avg_qintercepted(ipy) = sum(cpoly%avg_qintercepted *cpoly%area)*poly_area_i
         cgrid%avg_qthroughfall(ipy) = sum(cpoly%avg_qthroughfall *cpoly%area)*poly_area_i
         cgrid%avg_sensible_gc(ipy)  = sum(cpoly%avg_sensible_gc  *cpoly%area)*poly_area_i
         cgrid%avg_sensible_ac(ipy)  = sum(cpoly%avg_sensible_ac  *cpoly%area)*poly_area_i

         !----- Average albedo values. ----------------------------------------------------!
         cgrid%avg_albedo         (ipy) = sum(cpoly%avg_albedo         * cpoly%area)       &
                                        * poly_area_i
         cgrid%avg_albedo_beam    (ipy) = sum(cpoly%avg_albedo_beam    * cpoly%area)       &
                                        * poly_area_i
         cgrid%avg_albedo_diffuse (ipy) = sum(cpoly%avg_albedo_diffuse * cpoly%area)       &
                                        * poly_area_i
         cgrid%avg_rlong_albedo   (ipy) = sum(cpoly%avg_rlong_albedo   * cpoly%area)       &
                                        * poly_area_i
         !---------------------------------------------------------------------------------!

         cgrid%avg_available_water(ipy) = sum(cpoly%avg_available_water  *cpoly%area)      &
                                        * poly_area_i

         !---------------------------------------------------------------------------------!
         !    Polygon average of canopy thermodynamic state.  We average the variables     !
         ! that are insensitive to changes in pressure (potential temperature, specific    !
         ! humidity, CO2 mixing ratio and pressure itself) by averaging over sites.  The   !
         ! other variables are found using the prescribed diagnostic equations.            !
         !---------------------------------------------------------------------------------!
         cgrid%avg_can_depth   (ipy) = sum(cpoly%avg_can_depth * cpoly%area) * poly_area_i
         cgrid%avg_can_theiv   (ipy) = sum(cpoly%avg_can_theiv * cpoly%area) * poly_area_i
         cgrid%avg_can_theta   (ipy) = sum(cpoly%avg_can_theta * cpoly%area) * poly_area_i
         cgrid%avg_can_shv     (ipy) = sum(cpoly%avg_can_shv   * cpoly%area) * poly_area_i
         cgrid%avg_can_co2     (ipy) = sum(cpoly%avg_can_co2   * cpoly%area) * poly_area_i
         cgrid%avg_can_prss    (ipy) = sum(cpoly%avg_can_prss  * cpoly%area) * poly_area_i
         cgrid%avg_can_temp    (ipy) = cgrid%avg_can_theta(ipy)                            &
                                     * (p00i * cgrid%avg_can_prss(ipy)) ** rocp
         cgrid%avg_can_rhos    (ipy) = idealdenssh(cgrid%avg_can_prss(ipy)                 &
                                                  ,cgrid%avg_can_temp(ipy)                 &
                                                  ,cgrid%avg_can_shv (ipy) )
         !---------------------------------------------------------------------------------!
         !    Similar to the site level, average mass, heat capacity and energy then find  !
         ! the average temperature and liquid water fraction.                              !
         !---------------------------------------------------------------------------------!
         cgrid%avg_sensible_gg  (:,ipy) = matmul(cpoly%avg_sensible_gg  , cpoly%area)      &
                                        * poly_area_i
         cgrid%avg_smoist_gg    (:,ipy) = matmul(cpoly%avg_smoist_gg    , cpoly%area)      &
                                        * poly_area_i
         cgrid%avg_transloss    (:,ipy) = matmul(cpoly%avg_transloss    , cpoly%area)      &
                                        * poly_area_i
         cgrid%aux_s            (:,ipy) = matmul(cpoly%aux_s            , cpoly%area)      &
                                        * poly_area_i
         cgrid%avg_soil_energy  (:,ipy) = matmul(cpoly%avg_soil_energy  , cpoly%area)      &
                                        * poly_area_i
         cgrid%avg_soil_water   (:,ipy) = matmul(cpoly%avg_soil_water   , cpoly%area)      &
                                        * poly_area_i
         cgrid%avg_soil_mstpot  (:,ipy) = matmul(cpoly%avg_soil_mstpot  , cpoly%area)      &
                                        * poly_area_i
         cgrid%avg_soil_rootfrac(:,ipy) = matmul(cpoly%avg_soil_rootfrac, cpoly%area)      &
                                        * poly_area_i


         do k=cgrid%lsl(ipy),nzg
            !------------------------------------------------------------------------------!
            !     Finding the average temperature and liquid fraction.  The polygon-level  !
            ! mean heat capacity was already found during the site loop.                   !
            !------------------------------------------------------------------------------!
            call qwtk(cgrid%avg_soil_energy(k,ipy),cgrid%avg_soil_water(k,ipy)*wdns        &
                     ,poly_avg_soil_hcap(k),cgrid%avg_soil_temp(k,ipy)                     &
                     ,cgrid%avg_soil_fracliq(k,ipy))
         end do
         cgrid%avg_soil_wetness(ipy) = sum(cpoly%avg_soil_wetness * cpoly%area)            &
                                     * poly_area_i

         !---------------------------------------------------------------------------------!
         !    Also using the same idea as the site-level: average energy, mass, and depth, !
         ! but find temperature and liquid fraction with the averaged values. Again, use   !
         ! the energy in J/m2 instead of J/kg to do the polygon averaging.                 !
         !---------------------------------------------------------------------------------!
         cgrid%avg_sfcw_mass(ipy)   = sum(cpoly%avg_sfcw_mass  * cpoly%area) * poly_area_i
         cgrid%avg_sfcw_depth(ipy)  = sum(cpoly%avg_sfcw_depth * cpoly%area) * poly_area_i
         cgrid%avg_sfcw_energy(ipy) = sum( cpoly%avg_sfcw_energy                           &
                                         * cpoly%avg_sfcw_mass * cpoly%area) * poly_area_i

         !----- Scale energy and find temp and fracliq if there is enogh mass -------------!
         if (cgrid%avg_sfcw_mass(ipy) > tiny_sfcwater_mass) then
            cgrid%avg_sfcw_energy(ipy) = cgrid%avg_sfcw_energy(ipy)                        &
                                       / cgrid%avg_sfcw_mass(ipy)
            call qtk(cgrid%avg_sfcw_energy(ipy),cgrid%avg_sfcw_tempk(ipy)                  &
                    ,cgrid%avg_sfcw_fracliq(ipy))
         else
            cgrid%avg_sfcw_mass(ipy)    = 0.
            cgrid%avg_sfcw_depth(ipy)   = 0.
            cgrid%avg_sfcw_energy(ipy)  = 0.
            cgrid%avg_sfcw_tempk(ipy)   = cgrid%avg_soil_temp(nzg,ipy)
            cgrid%avg_sfcw_fracliq(ipy) = cgrid%avg_soil_fracliq(nzg,ipy)
         end if
         !---------------------------------------------------------------------------------!
         !    Similar to site level, compute mean leaf and stem internal energy, water     !
         ! mass, and heat capacity.  If there is some heat capacity, then find the mean    !
         ! temperature, otherwise just assume to be the canopy temperature.                !
         !---------------------------------------------------------------------------------!
         cgrid%avg_leaf_energy(ipy) = sum(cpoly%avg_leaf_energy * cpoly%area) * poly_area_i
         cgrid%avg_leaf_water(ipy)  = sum(cpoly%avg_leaf_water  * cpoly%area) * poly_area_i
         cgrid%avg_leaf_hcap(ipy)   = sum(cpoly%avg_leaf_hcap   * cpoly%area) * poly_area_i
         cgrid%avg_wood_energy(ipy) = sum(cpoly%avg_wood_energy * cpoly%area) * poly_area_i
         cgrid%avg_wood_water(ipy)  = sum(cpoly%avg_wood_water  * cpoly%area) * poly_area_i
         cgrid%avg_wood_hcap(ipy)   = sum(cpoly%avg_wood_hcap   * cpoly%area) * poly_area_i
         if (cgrid%avg_leaf_hcap(ipy) > 0.) then
            call qwtk(cgrid%avg_leaf_energy(ipy),cgrid%avg_leaf_water(ipy)                   &
                     ,cgrid%avg_leaf_hcap(ipy),cgrid%avg_leaf_temp(ipy)                      &
                     ,cgrid%avg_leaf_fliq(ipy))
         else
            cgrid%avg_leaf_temp(ipy) = cgrid%avg_can_temp(ipy)
            if (cgrid%avg_can_temp(ipy) > 0.0) then
               cgrid%avg_leaf_fliq(ipy) = 1.0
            elseif (cgrid%avg_can_temp(ipy) == 0.0) then
               cgrid%avg_leaf_fliq(ipy) = 0.5
            else
               cgrid%avg_leaf_fliq(ipy) = 0.0
            end if
         end if
         if (cgrid%avg_wood_hcap(ipy) > 0.) then
            call qwtk(cgrid%avg_wood_energy(ipy),cgrid%avg_wood_water(ipy)                   &
                     ,cgrid%avg_wood_hcap(ipy),cgrid%avg_wood_temp(ipy)                      &
                     ,cgrid%avg_wood_fliq(ipy))
         else
            cgrid%avg_wood_temp(ipy) = cgrid%avg_can_temp(ipy)
            if (cgrid%avg_can_temp(ipy) > 0.0) then
               cgrid%avg_wood_fliq(ipy) = 1.0
            elseif (cgrid%avg_can_temp(ipy) == 0.0) then
               cgrid%avg_wood_fliq(ipy) = 0.5
            else
               cgrid%avg_wood_fliq(ipy) = 0.0
            end if
         end if
         !---------------------------------------------------------------------------------!





         !---------------------------------------------------------------------------------!
         !    Compute the average amongst all surfaces (soil, temporary surface water, and !
         ! vegetation, the last two only if they really exist).  All energy terms are      !
         ! converted to J/m2, all water terms to kg/m2, and the heat capacities of every-  !
         ! thing that is not water is in J/m2/K.                                           !
         !---------------------------------------------------------------------------------!
         skin_energy = cgrid%avg_leaf_energy(ipy)                                          &
                     + cgrid%avg_wood_energy(ipy)                                          &
                     + cgrid%avg_sfcw_energy(ipy) * cgrid%avg_sfcw_mass(ipy)               &
                     + cgrid%avg_soil_energy(nzg,ipy) * dslz(nzg)
         skin_water  = cgrid%avg_leaf_water(ipy)                                           &
                     + cgrid%avg_wood_water(ipy)                                           &
                     + cgrid%avg_sfcw_mass(ipy)                                            &
                     + cgrid%avg_soil_water(nzg,ipy) * dslz(nzg) * wdns
         skin_hcap   = cgrid%avg_leaf_hcap(ipy)                                            &
                     + cgrid%avg_wood_hcap(ipy)                                            &
                     + poly_avg_soil_hcap(nzg) * dslz(nzg)
         call qwtk(skin_energy,skin_water,skin_hcap,cgrid%avg_skin_temp(ipy),skin_fliq)
         !---------------------------------------------------------------------------------!

      end do polyloop
   end do gridloop

   return
end subroutine spatial_averages
!==========================================================================================!
!==========================================================================================!





