module radiate_driver
   contains

   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will control the two-stream radiation scheme.  This is called     !
   ! every step, but not every sub-step.                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_radiation(cgrid)
      use ed_misc_coms          , only : current_time          & ! intent(in)
                                       , radfrq                & ! intent(in)
                                       , dtlsm                 ! ! intent(in)
      use ed_state_vars         , only : edtype                & ! structure
                                       , polygontype           & ! structure
                                       , sitetype              & ! structure
                                       , patchtype             ! ! structure
      use canopy_radiation_coms , only : par_beam_norm         & ! intent(in)
                                       , par_diff_norm         & ! intent(in)
                                       , nir_beam_norm         & ! intent(in)
                                       , nir_diff_norm         & ! intent(in)
                                       , cosz_min              & ! intent(in)
                                       , rshort_twilight_min   ! ! intent(in)
      use consts_coms           , only : pio180                ! ! intent(in)
      use grid_coms             , only : nzg                   & ! intent(in)
                                       , nzs                   ! ! intent(in)
      use radiate_utils         , only : angle_of_incid        & ! sub-routine
                                       , ed_zen                & ! function
                                       , update_rad_avg        ! ! function
      !$ use omp_lib
      implicit none
      !----- Argument. --------------------------------------------------------------------!
      type(edtype)     , target   :: cgrid
      !----- Local variables. -------------------------------------------------------------!
      type(polygontype), pointer  :: cpoly
      type(sitetype)   , pointer  :: csite
      type(patchtype)  , pointer  :: cpatch
      real                        :: rshort_tot
      integer                     :: maxcohort
      integer                     :: ipy
      integer                     :: isi
      integer                     :: ipa
      logical                     :: daytime
      logical                     :: twilight
      real                        :: hrangl
      real                        :: sloperad
      real                        :: aspectrad
      real                        :: sum_norm
      integer                     :: ibuff
      !------------------------------------------------------------------------------------!


      !----- Check whether it is time to update radiative fluxes and heating rates --------!
      if (mod(current_time%time + .001,radfrq) < dtlsm) then

         !----- Loop over polygons and sites. ---------------------------------------------!

         polyloop: do ipy = 1,cgrid%npolygons

            !----- Find the solar zenith angle [cosz] --------------------------------------!
            cgrid%cosz(ipy) = ed_zen(cgrid%lon(ipy),cgrid%lat(ipy),current_time)
            !------------------------------------------------------------------------------!

            cpoly => cgrid%polygon(ipy)

            siteloop: do isi = 1,cpoly%nsites

               csite => cpoly%site(isi)

               !---------------------------------------------------------------------------!
               !     Update angle of incidence.                                            !
               !---------------------------------------------------------------------------!
               hrangl    = 15. * pio180                                                    &
                         * (mod(current_time%time + cgrid%lon(ipy) / 15. + 24., 24.) - 12.)
               sloperad  = cpoly%slope(isi)  * pio180
               aspectrad = cpoly%aspect(isi) * pio180
               call angle_of_incid(cpoly%cosaoi(isi),cgrid%cosz(ipy),hrangl                &
                                  ,sloperad,aspectrad)
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    Find the two logicals, that will tell which part of the day we are at  !
               ! least.                                                                    !
               !---------------------------------------------------------------------------!
               daytime  = cpoly%cosaoi(isi) > cosz_min .and.                               &
                          cpoly%met(isi)%rshort > rshort_twilight_min
               twilight = cpoly%met(isi)%rshort > rshort_twilight_min
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !      Update the daylight length and nighttime flag.                       !
               !---------------------------------------------------------------------------!
               cpoly%nighttime(isi)              = .not. twilight
               if (twilight) cpoly%daylight(isi) = cpoly%daylight(isi) + radfrq
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !    Loop over subgrid-scale patches.  These routines can be done as        !
               ! arrays.                                                                   !
               !---------------------------------------------------------------------------!
               maxcohort = 1
               do ipa = 1,csite%npatches
                  cpatch=>csite%patch(ipa)
                  if ( cpatch%ncohorts>maxcohort ) maxcohort = cpatch%ncohorts
               end do
               !---------------------------------------------------------------------------!


               !---------------------------------------------------------------------------!
               !    MLO. Changed this part of the code so it also works with the hori-     !
               ! zontal shading.  Essentially I took the loop through patches out of       !
               ! sfcrad_ed and scale_ed_radiation, so each patch could use different       !
               ! normalisation factors.                                                    !
               !---------------------------------------------------------------------------!
               !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(                                  &
               !$OMP  ibuff,par_beam_norm,par_diff_norm,nir_beam_norm,nir_diff_norm        &
               !$OMP ,rshort_tot,sum_norm)
               patchloop: do ipa=1,csite%npatches

                  ibuff = 1
                  !$ ibuff = OMP_get_thread_num()+1

                  !------------------------------------------------------------------------!
                  !     Copy radiation components to the local variables; they will be     !
                  ! normalised in the 'if' block that follows this block.  In case IHRZRAD !
                  ! is 1 or 2, the correction term 'fbeam' changes the amount of incoming  !
                  ! radiation to account for lateral shading or lateral illumination.      !
                  ! Otherwise, fbeam is always 1.                                          !
                  !------------------------------------------------------------------------!
                  par_beam_norm = dble(cpoly%met(isi)%par_beam*csite%fbeam(ipa))
                  par_diff_norm = dble(cpoly%met(isi)%par_diffuse              )
                  nir_beam_norm = dble(cpoly%met(isi)%nir_beam*csite%fbeam(ipa))
                  nir_diff_norm = dble(cpoly%met(isi)%nir_diffuse              )
                  rshort_tot    = cpoly%met(isi)%par_beam * csite%fbeam(ipa)               &
                                + cpoly%met(isi)%par_diffuse                               &
                                + cpoly%met(isi)%nir_beam * csite%fbeam(ipa)               &
                                + cpoly%met(isi)%nir_diffuse
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !      In case the angle of incidence is too high (i.e., its cosine is   !
                  ! too close to zero), eliminate all direct radiation.   This is differ-  !
                  ! ent from the cosine of the zenith angle because mountains can hide the !
                  ! sun even when it is still above the horizon.                           !
                  !------------------------------------------------------------------------!
                  if (daytime) then
                     !---------------------------------------------------------------------!
                     !    Daytime.  Normalise the four components by total radiation.  In  !
                     ! case IHRZRAD = 0 or 4, incoming radiation is simply the site level. !
                     ! Otherwise, incoming direct radiation may be scaled to account for   !
                     ! lateral shading or additional lateral illumination based on crown   !
                     ! closure.  In both cases the canopy radiation model uses normalised  !
                     ! profiles to find the vertical structure.                            !
                     !---------------------------------------------------------------------!
                     par_beam_norm = max( 1.d-5, par_beam_norm / dble(rshort_tot) )
                     par_diff_norm = max( 1.d-5, par_diff_norm / dble(rshort_tot) )
                     nir_beam_norm = max( 1.d-5, nir_beam_norm / dble(rshort_tot) )
                     nir_diff_norm = max( 1.d-5, nir_diff_norm / dble(rshort_tot) )
                     !---------------------------------------------------------------------!
                  elseif (twilight) then
                     !---------------------------------------------------------------------!
                     !     Twilight, if for some reason there is any direct radiation,     !
                     ! copy it to diffuse light.                                           !
                     !---------------------------------------------------------------------!
                     par_diff_norm = max( 1.d-5                                            &
                                        , (par_beam_norm+par_diff_norm) / dble(rshort_tot))
                     nir_diff_norm = max( 1.d-5                                            &
                                        , (nir_beam_norm+nir_diff_norm) / dble(rshort_tot))
                     par_beam_norm = 1.d-5
                     nir_beam_norm = 1.d-5
                     !---------------------------------------------------------------------!
                  else 
                     !---------------------------------------------------------------------!
                     !     Night-time, nothing will happen, fill split equally to the 4    !
                     ! components.                                                         !
                     !---------------------------------------------------------------------!
                     rshort_tot    = 0.0
                     par_beam_norm = 2.5d-1
                     par_diff_norm = 2.5d-1
                     nir_beam_norm = 2.5d-1
                     nir_diff_norm = 2.5d-1
                  end if
                  !------------------------------------------------------------------------!



                  !------------------------------------------------------------------------!
                  !     Because we must tweak the radiation so none of the terms are zero, !
                  ! we must correct the normalised radiation variables so they add up to   !
                  ! one.                                                                   !
                  !------------------------------------------------------------------------!
                  sum_norm      = par_beam_norm + par_diff_norm                            &
                                + nir_beam_norm + nir_diff_norm
                  par_beam_norm = par_beam_norm / sum_norm
                  par_diff_norm = par_diff_norm / sum_norm
                  nir_beam_norm = nir_beam_norm / sum_norm
                  nir_diff_norm = nir_diff_norm / sum_norm
                  !------------------------------------------------------------------------!


                  !----- Get normalised radiative transfer information. -------------------!
                  call sfcrad_ed(cpoly%cosaoi(isi),csite,ipa,ibuff,nzg,nzs                 &
                                ,cpoly%ntext_soil(:,isi),cpoly%ncol_soil(isi)              &
                                ,cpoly%met(isi)%rlong,twilight)
                  !------------------------------------------------------------------------!



                  !----- Scale radiation back to the actual values. -----------------------!
                  call scale_ed_radiation(rshort_tot,cpoly%met(isi)%rlong                  &
                                         ,cpoly%nighttime(isi),csite,ipa)
                  !------------------------------------------------------------------------!
               end do patchloop
               !$OMP END PARALLEL DO
               !---------------------------------------------------------------------------!
            end do siteloop
            !------------------------------------------------------------------------------!
         end do polyloop
         !---------------------------------------------------------------------------------!

         !----- Update the average radiation for phenology. -------------------------------!
         call update_rad_avg(cgrid)
         !---------------------------------------------------------------------------------!

      end if

      return
   end subroutine canopy_radiation
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine will drive the distribution of radiation among crowns, snow       !
   ! layers, and soil.                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine sfcrad_ed(cosaoi,csite,ipa,ibuff,mzg,mzs,ntext_soil,ncol_soil,rlong,twilight)

      use ed_state_vars        , only : sitetype             & ! structure
                                      , patchtype            ! ! structure
      use canopy_layer_coms    , only : crown_mod            & ! intent(in)
                                      , tai_lyr_max          ! ! intent(in)
      use canopy_radiation_coms, only : icanrad              & ! intent(in)
                                      , clumping_factor      & ! intent(in)
                                      , par_beam_norm        & ! intent(in)
                                      , par_diff_norm        & ! intent(in)
                                      , nir_beam_norm        & ! intent(in)
                                      , nir_diff_norm        & ! intent(in)
                                      , leaf_scatter_vis     & ! intent(in)
                                      , wood_scatter_vis     & ! intent(in)
                                      , leaf_scatter_nir     & ! intent(in)
                                      , wood_scatter_nir     & ! intent(in)
                                      , leaf_emiss_tir       & ! intent(in)
                                      , wood_emiss_tir       & ! intent(in)
                                      , snow_albedo_vis      & ! intent(in)
                                      , snow_albedo_nir      & ! intent(in)
                                      , snow_emiss_tir       & ! intent(in)
                                      , radscr               ! ! intent(inout)
      use soil_coms            , only : soil                 & ! intent(in)
                                      , soilcol              ! ! intent(in)
      use consts_coms          , only : stefan               & ! intent(in)
                                      , lnexp_max            ! ! intent(in)
      use ed_max_dims          , only : n_pft                & ! intent(in)
                                      , n_radprof            ! ! intent(in)
      use allometry            , only : h2crownbh            ! ! intent(in)
      use ed_misc_coms         , only : ibigleaf             ! ! intent(in)
      use old_twostream_rad    , only : old_lw_two_stream    & ! sub-routine
                                      , old_sw_two_stream    ! ! sub-routine
      use multiple_scatter     , only : lw_multiple_scatter  & ! sub-routine
                                      , sw_multiple_scatter  ! ! sub-routine
      use twostream_rad        , only : lw_two_stream        & ! sub-routine
                                      , sw_two_stream        ! ! sub-routine
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                  , target      :: csite
      integer                         , intent(in)  :: ipa
      integer                         , intent(in)  :: mzg
      integer                         , intent(in)  :: mzs
      integer         , dimension(mzg), intent(in)  :: ntext_soil
      integer                         , intent(in)  :: ncol_soil
      real                            , intent(in)  :: rlong
      real                            , intent(in)  :: cosaoi
      logical                         , intent(in)  :: twilight
      integer                         , intent(in)  :: ibuff
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer                     :: cpatch
      integer                                       :: il
      integer                                       :: ico
      integer                                       :: ipft
      integer                                       :: cohort_count
      integer                                       :: nsoil
      integer                                       :: colour
      integer                                       :: k
      integer                                       :: ksn
      integer                                       :: tuco_leaf
      real                                          :: fcpct
      real                                          :: albedo_soil_par
      real                                          :: albedo_soil_nir
      real                                          :: albedo_damp_par
      real                                          :: albedo_damp_nir
      real                                          :: albedo_sfcw_par
      real                                          :: albedo_sfcw_nir
      real                                          :: rad_sfcw_par
      real                                          :: rad_sfcw_nir
      real                                          :: sfcw_odepth
      real                                          :: fractrans_par
      real                                          :: fractrans_nir
      real            , dimension(mzs)              :: abs_sfcw_par
      real            , dimension(mzs)              :: abs_sfcw_nir
      real                                          :: abs_ground_par
      real                                          :: abs_ground_nir
      real                                          :: albedo_ground_par
      real                                          :: albedo_ground_nir
      real                                          :: downward_par_below_beam
      real                                          :: downward_nir_below_beam
      real                                          :: downward_par_below_diffuse
      real                                          :: upward_par_above_diffuse
      real                                          :: downward_nir_below_diffuse
      real                                          :: upward_nir_above_diffuse
      real                                          :: T_surface
      real                                          :: emissivity
      real                                          :: downward_lw_below
      real                                          :: upward_lw_below
      real                                          :: upward_lw_above
      real                                          :: downward_rshort_below_beam
      real                                          :: downward_rshort_below_diffuse
      real                                          :: upward_rshort_above_diffuse
      real                                          :: surface_netabs_longwave
      real                                          :: nir_v_beam
      real                                          :: nir_v_diffuse
      real                                          :: wleaf_vis
      real                                          :: wleaf_nir
      real                                          :: wleaf_tir
      real                                          :: wwood_vis
      real                                          :: wwood_nir
      real                                          :: wwood_tir
      real                                          :: bl_lai_each
      real                                          :: bl_wai_each
     !----- External function. -----------------------------------------------------------!
      real            , external                    :: sngloff
      !----- Local constants. -------------------------------------------------------------!
      real(kind=8)    , parameter                   :: tiny_offset = 1.d-20
      !------------------------------------------------------------------------------------!


    
      !----- Link current cohort ----------------------------------------------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Initialise the variables.                                                      !
      !------------------------------------------------------------------------------------!
      radscr(ibuff)%pft_array(:)                 = -1
      radscr(ibuff)%leaf_temp_array(:)           = 0.d0
      radscr(ibuff)%wood_temp_array(:)           = 0.d0
      radscr(ibuff)%lai_array(:)                 = 0.d0
      radscr(ibuff)%wai_array(:)                 = 0.d0
      radscr(ibuff)%CA_array(:)                  = 0.d0
      radscr(ibuff)%htop_array(:)                = 0.d0
      radscr(ibuff)%hbot_array(:)                = 0.d0
      radscr(ibuff)%par_level_beam(:)            = 0.d0
      radscr(ibuff)%par_level_diffd(:)           = 0.d0
      radscr(ibuff)%par_level_diffu(:)           = 0.d0
      radscr(ibuff)%light_level_array(:)         = 0.d0
      radscr(ibuff)%light_beam_level_array(:)    = 0.d0
      radscr(ibuff)%light_diff_level_array(:)    = 0.d0
      radscr(ibuff)%par_v_beam_array(:)          = 0.
      radscr(ibuff)%rshort_v_beam_array(:)       = 0.
      radscr(ibuff)%par_v_diffuse_array(:)       = 0.
      radscr(ibuff)%rshort_v_diffuse_array(:)    = 0.
      radscr(ibuff)%lw_v_array(:)                = 0.
      radscr(ibuff)%radprof_array(:,:)           = 0.
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Cohort_count is the number of cohorts that affect the radiation balance (i.e.  !
      ! those which are flagged as resolvable.                                             !
      !------------------------------------------------------------------------------------!
      cohort_count = 0
      tuco_leaf    = 0
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Initialise par_l_ as zero just in case it is night time or if there is no      !
      ! resolvable cohort.                                                                 !
      !------------------------------------------------------------------------------------!
      csite%par_l_beam_max   (ipa) = 0.0
      csite%par_l_diffuse_max(ipa) = 0.0
      csite%par_l_max        (ipa) = 0.0
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Transfer information from linked lists to arrays.  Here we must check          !
      ! whether is running as a true size-and-age structure model, or as big leaf.         !
      !------------------------------------------------------------------------------------!
      select case (ibigleaf)
      case (0)

         !---------------------------------------------------------------------------------!
         !     Size- and age-structure.  Each layer in the radiation will correspond to    !
         ! one cohort.  Unusually, here we go from shortest to tallest, as required by the !
         ! radiation schemes.                                                              !
         !---------------------------------------------------------------------------------!
         do ico = cpatch%ncohorts,1,-1

            !----- Initialize values. -----------------------------------------------------!
            cpatch%par_l             (ico)    = 0.0
            cpatch%par_l_beam        (ico)    = 0.0
            cpatch%par_l_diffuse     (ico)    = 0.0

            cpatch%rshort_l          (ico)    = 0.0
            cpatch%rshort_l_beam     (ico)    = 0.0
            cpatch%rshort_l_diffuse  (ico)    = 0.0

            cpatch%rlong_l           (ico)    = 0.0

            cpatch%rshort_w          (ico)    = 0.0
            cpatch%rshort_w_beam     (ico)    = 0.0
            cpatch%rshort_w_diffuse  (ico)    = 0.0

            cpatch%rlong_w           (ico)    = 0.0

            cpatch%light_level       (ico)    = 0.0
            cpatch%light_level_beam  (ico)    = 0.0
            cpatch%light_level_diff  (ico)    = 0.0

            cpatch%rad_profile     (:,ico)    = 0.0

            if (cpatch%leaf_resolvable(ico) .or. cpatch%wood_resolvable(ico)) then

               cohort_count                          = cohort_count + 1
               radscr(ibuff)%pft_array(cohort_count) = cpatch%pft(ico)

               !---------------------------------------------------------------------------!
               !     Here we only tell the true LAI if the leaf is resolvable, and the     !
               ! true WAI if the wood is resolvable.  Also, for photosynthesis, we must    !
               ! keep track of the tallest cohort that has leaves (we track the array      !
               ! counters because we will extract the information directly from the        !
               ! arrays.                                                                   !
               !---------------------------------------------------------------------------!
               if (cpatch%leaf_resolvable(ico)) then
                  tuco_leaf                             = cohort_count
                  radscr(ibuff)%lai_array(cohort_count) = dble(cpatch%lai(ico))
               else
                  radscr(ibuff)%lai_array(cohort_count) = 0.d0
               end if
               if (cpatch%wood_resolvable(ico)) then
                  radscr(ibuff)%wai_array(cohort_count) = dble(cpatch%wai(ico))
               else
                  radscr(ibuff)%wai_array(cohort_count) = 0.d0
               end if
               !---------------------------------------------------------------------------!

               ! TEST REMOVING THESE ZERO CALLS (REDUNDANT)

               radscr(ibuff)%leaf_temp_array(cohort_count) = dble(cpatch%leaf_temp(ico))
               radscr(ibuff)%wood_temp_array(cohort_count) = dble(cpatch%wood_temp(ico))
               radscr(ibuff)%rshort_v_beam_array   (  cohort_count) = 0.0
               radscr(ibuff)%par_v_beam_array      (  cohort_count) = 0.0
               radscr(ibuff)%rshort_v_diffuse_array(  cohort_count) = 0.0
               radscr(ibuff)%par_v_diffuse_array   (  cohort_count) = 0.0
               radscr(ibuff)%lw_v_array            (  cohort_count) = 0.0
               radscr(ibuff)%radprof_array         (:,cohort_count) = 0.0
               radscr(ibuff)%par_level_beam        (  cohort_count) = 0.d0
               radscr(ibuff)%par_level_diffd       (  cohort_count) = 0.d0
               radscr(ibuff)%par_level_diffu       (  cohort_count) = 0.d0
               radscr(ibuff)%light_level_array     (  cohort_count) = 0.d0
               radscr(ibuff)%light_beam_level_array(  cohort_count) = 0.d0
               radscr(ibuff)%light_diff_level_array(  cohort_count) = 0.d0
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !      Decide whether to assume infinite crown, or the crown area allometry !
               ! method as in Dietze and Clark (2008).                                     !
               !---------------------------------------------------------------------------!
               select case (crown_mod)
               case (0)
                  radscr(ibuff)%CA_array(cohort_count) = 1.d0
               case (1)
                  radscr(ibuff)%CA_array(cohort_count) = dble(cpatch%crown_area(ico))
               end select
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!
      case (1)
         !---------------------------------------------------------------------------------!
         !     Big-leaf solver.  We have either 0 (desert) or 1 cohort.  In case of zero,  !
         ! we bypass the entire cohort loop, otherwise we check how many layers we will    !
         ! create.                                                                         !
         !---------------------------------------------------------------------------------!
         if (cpatch%ncohorts > 0) then
            !----- Initialize values. -----------------------------------------------------!
            cpatch%par_l(1)                 = 0.0
            cpatch%par_l_beam(1)            = 0.0
            cpatch%par_l_diffuse(1)         = 0.0
            cpatch%rshort_l(1)              = 0.0
            cpatch%rshort_l_beam(1)         = 0.0
            cpatch%rshort_l_diffuse(1)      = 0.0
            cpatch%rlong_l(1)               = 0.0
            cpatch%rshort_w(1)              = 0.0
            cpatch%rshort_w_beam(1)         = 0.0
            cpatch%rshort_w_diffuse(1)      = 0.0
            cpatch%rlong_w(1)               = 0.0
            cpatch%rad_profile     (:,1)    = 0.0
            cpatch%light_level     (1)      = 0.0
            cpatch%light_level_beam(1)      = 0.0
            cpatch%light_level_diff(1)      = 0.0

            !------------------------------------------------------------------------------!
            !     Check whether the cohort is resolvable or not.                           !
            !------------------------------------------------------------------------------!
           if (cpatch%leaf_resolvable(1) .or. cpatch%wood_resolvable(1)) then

               !---------------------------------------------------------------------------!
               !     Here we only tell the true LAI if the leaf is resolvable, and the     !
               ! true WAI if the wood is resolvable.  Also, for photosynthesis, we must    !
               ! keep track of the tallest cohort that has leaves (we track the array      !
               ! counters because we will extract the information directly from the        !
               ! arrays.                                                                   !
               !---------------------------------------------------------------------------!
               if (cpatch%leaf_resolvable(1) .and. cpatch%wood_resolvable(1)) then
                  cohort_count = ceiling( (cpatch%lai(1) + cpatch%wai(1)) / tai_lyr_max )
                  bl_lai_each  = cpatch%lai(1) / real(cohort_count)
                  bl_wai_each  = cpatch%wai(1) / real(cohort_count) 
                  tuco_leaf    = cohort_count
               elseif (cpatch%leaf_resolvable(1)) then
                  cohort_count = ceiling( cpatch%lai(1) / tai_lyr_max )
                  bl_lai_each  = cpatch%lai(1) / real(cohort_count)
                  bl_wai_each  = 0.0
                  tuco_leaf    = cohort_count
               elseif (cpatch%wood_resolvable(1)) then
                  cohort_count = ceiling( cpatch%wai(1) / tai_lyr_max )
                  bl_lai_each  = 0.0
                  bl_wai_each  = cpatch%wai(1) / real(cohort_count)
               end if
               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !     Loop over all layers, and assign equal amounts of LAI and WAI such    !
               ! that they add back to the total amount.  We impose crown model off for    !
               ! big leaf, so CA_array must be always set to 1.                            !
               !---------------------------------------------------------------------------!
               do ico = 1, cohort_count
                  radscr(ibuff)%pft_array                (ico) = cpatch%pft(1)
                  radscr(ibuff)%lai_array                (ico) = dble(bl_lai_each)
                  radscr(ibuff)%wai_array                (ico) = dble(bl_wai_each)
                  radscr(ibuff)%CA_array                 (ico) = 1.d0
                  radscr(ibuff)%leaf_temp_array          (ico) = dble(cpatch%leaf_temp(1))
                  radscr(ibuff)%wood_temp_array          (ico) = dble(cpatch%wood_temp(1))
                  radscr(ibuff)%rshort_v_beam_array      (ico) = 0.0
                  radscr(ibuff)%par_v_beam_array         (ico) = 0.0
                  radscr(ibuff)%rshort_v_diffuse_array   (ico) = 0.0
                  radscr(ibuff)%lw_v_array             (  ico) = 0.0
                  radscr(ibuff)%radprof_array          (:,ico) = 0.0
                  radscr(ibuff)%par_v_diffuse_array      (ico) = 0.0
                  radscr(ibuff)%par_level_beam           (ico) = 0.d0
                  radscr(ibuff)%par_level_diffu          (ico) = 0.d0
                  radscr(ibuff)%par_level_diffd          (ico) = 0.d0
                  radscr(ibuff)%light_level_array        (ico) = 0.d0
                  radscr(ibuff)%light_beam_level_array   (ico) = 0.d0
                  radscr(ibuff)%light_diff_level_array   (ico) = 0.d0
               end do
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      csite%rshort_s_diffuse(:,ipa) = 0.0
      csite%rshort_s_beam   (:,ipa) = 0.0
      csite%par_s_diffuse   (:,ipa) = 0.0
      csite%par_s_beam      (:,ipa) = 0.0
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the ground albedo as a function of soil water relative moisture of the    !
      ! top layer.                                                                         !
      !------------------------------------------------------------------------------------!
      nsoil  = ntext_soil(mzg)
      colour = ncol_soil
      select case (nsoil)
      case (13)
         !----- Bedrock, use constants soil value for granite. ----------------------------!
         albedo_soil_par = soilcol(colour)%alb_vis_dry
         albedo_soil_nir = soilcol(colour)%alb_nir_dry
         !----- Damp soil, for temporary surface water albedo. ----------------------------!
         albedo_damp_par = albedo_soil_par
         albedo_damp_nir = albedo_damp_nir
      case (12)
         !----- Peat, follow McCumber and Pielke (1981). ----------------------------------!
         fcpct = csite%soil_water(mzg,ipa) / soil(nsoil)%slmsts
         albedo_soil_par = max (0.07, 0.14 * (1.0 - fcpct))
         albedo_soil_nir = albedo_soil_par
         !----- Damp soil, for temporary surface water albedo. ----------------------------!
         albedo_damp_par = 0.14
         albedo_damp_nir = 0.14
      case default
         select case (colour)
         case (21)
            !------------------------------------------------------------------------------!
            !     ED-2.1 soil colour.  Also, we use the ED-2.1 default method to determine !
            ! the albedo.                                                                  !
            !------------------------------------------------------------------------------!
            fcpct           = csite%soil_water(mzg,ipa) / soil(nsoil)%slmsts
            albedo_soil_par = max(0.14,0.31-0.34*fcpct)
            albedo_soil_nir = albedo_soil_par
            !----- Damp soil, for temporary surface water albedo. -------------------------!
            albedo_damp_par = 0.14
            albedo_damp_nir = 0.14
         case default
            !------------------------------------------------------------------------------!
            !      Other soils, we use the soil numbers from CLM-4.  The colour class must !
            ! be given at RAMSIN.  At this point the value is the same for all points, but !
            ! in the future we may read their files if the results are promising.          !
            !------------------------------------------------------------------------------!
            fcpct           = max(0., 0.11 - 0.40 * csite%soil_water(mzg,ipa))
            albedo_soil_par = min(soilcol(colour)%alb_vis_dry                              &
                                 ,soilcol(colour)%alb_vis_wet  + fcpct)
            albedo_soil_nir = min(soilcol(colour)%alb_nir_dry                              &
                                 ,soilcol(colour)%alb_nir_wet  + fcpct)
            !----- Damp soil, for temporary surface water albedo. -------------------------!
            fcpct           = max(0., 0.11 - 0.40 * soil(nsoil)%slmsts)
            albedo_damp_par = min(soilcol(colour)%alb_vis_dry                              &
                                 ,soilcol(colour)%alb_vis_wet  + fcpct)
            albedo_damp_nir = min(soilcol(colour)%alb_nir_dry                              &
                                 ,soilcol(colour)%alb_nir_wet  + fcpct)
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Decide what is our surface temperature.  When the soil is exposed, then that   !
      ! is the surface temperature.  Otherwise, we pick the temporary surface water or     !
      ! snow layer.                                                                        !
      !------------------------------------------------------------------------------------!
      rad_sfcw_par      = 1.0
      rad_sfcw_nir      = 1.0
      albedo_ground_par = 1.0
      albedo_ground_nir = 1.0
      abs_sfcw_par      = 0.0
      abs_sfcw_nir      = 0.0
      ksn               = csite%nlev_sfcwater(ipa)
      if (ksn == 0) then
         emissivity = soilcol(colour)%emiss_tir
         T_surface  = csite%soil_tempk(mzg,ipa)
      else
         !---------------------------------------------------------------------------------!
         !      Sfcwater albedo ALS ranges from wet-soil value for all-liquid to typical   !
         ! snow albedo for ice.  In the future, we should consider a more realistic snow   !
         ! albedo model that takes snow age and snow melt into account.                    !
         !                                                                                 !
         !  Potential starting points:                                                     !
         !                                                                                 !
         !  Roesch, A., et al., 2002: Comparison of spectral surface albedos and their     !
         !      impact on the general circulation model simulated surface climate.         !
         !      J. Geophys. Res.-Atmosph., 107(D14), 4221, 10.1029/2001JD000809.           !
         !                                                                                 !
         !  Oleson, K.W., et al., 2010: Technical description of version 4.0 of the        !
         !      Community Land Model (CLM). NCAR Technical Note NCAR/TN-478+STR.           !
         !                                                                                 !
         !---------------------------------------------------------------------------------!
         albedo_sfcw_par = albedo_damp_par + csite%sfcwater_fracliq(ksn,ipa)               &
                                           * ( snow_albedo_vis - albedo_damp_par )
         albedo_sfcw_nir = albedo_damp_nir + csite%sfcwater_fracliq(ksn,ipa)               &
                                           * ( snow_albedo_nir - albedo_damp_nir )
         !---------------------------------------------------------------------------------!


         !----- Fraction shortwave absorbed into sfcwater + soil. -------------------------!
         rad_sfcw_par = 1.0 - albedo_sfcw_par
         rad_sfcw_nir = 1.0 - albedo_sfcw_nir
         do k = ksn,1,-1
            !------------------------------------------------------------------------------!
            !      Fractrans is fraction of shortwave entering each sfcwater layer that    !
            ! gets transmitted through that layer.                                         !
            !------------------------------------------------------------------------------!
            sfcw_odepth   = min( lnexp_max,   20.0 * csite%sfcwater_depth(k,ipa)           &
                                                   / csite%snowfac         (ipa) )
            fractrans_par = exp( - sfcw_odepth )
            fractrans_nir = fractrans_par
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      abs_sfcw_???(k) is fraction of total incident shortwave (at top of top  !
            ! sfcwater layer) that is absorbed in each sfcwater layer.                     !
            !------------------------------------------------------------------------------!
            abs_sfcw_par(k) = rad_sfcw_par * (1.0 - fractrans_par) * csite%snowfac(ipa)
            abs_sfcw_nir(k) = rad_sfcw_nir * (1.0 - fractrans_nir) * csite%snowfac(ipa)
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      Rad is fraction of total incident shortwave (at top of top sfcwater     !
            ! layer) that remains at bottom of current sfcwater layer.                     !
            !------------------------------------------------------------------------------!
            rad_sfcw_par = rad_sfcw_par * fractrans_par
            rad_sfcw_nir = rad_sfcw_nir * fractrans_nir
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      Albedo_ground will ultimately be the albedo of the soil+sfcwater.  So   !
            ! subtract out whatever is being absorbed by sfcwater.                         !
            !------------------------------------------------------------------------------!
            albedo_ground_par = albedo_ground_par - abs_sfcw_par(k)
            albedo_ground_nir = albedo_ground_nir - abs_sfcw_nir(k)
            !------------------------------------------------------------------------------!
         end do
         !---------------------------------------------------------------------------------!


         !----- Long wave parameter if sfcwater exists. -----------------------------------!
         emissivity = snow_emiss_tir            *        csite%snowfac(ipa)                &
                    + soilcol(colour)%emiss_tir * ( 1. - csite%snowfac(ipa) )
         T_surface  = sqrt(sqrt( ( csite%sfcwater_tempk (ksn,ipa)** 4                      &
                                 * snow_emiss_tir            *        csite%snowfac(ipa)   &
                                 + csite%soil_tempk     (mzg,ipa)** 4                      &
                                 * soilcol(colour)%emiss_tir * ( 1. - csite%snowfac(ipa))) &
                               / emissivity ) )
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     This is the fraction of below-canopy radiation that is absorbed by the ground. !
      !------------------------------------------------------------------------------------!
      abs_ground_par = (1.0 - albedo_soil_par)                                             &
                     * (1.0 - csite%snowfac(ipa) + csite%snowfac(ipa) * rad_sfcw_par )
      abs_ground_nir = (1.0 - albedo_soil_nir)                                             &
                     * (1.0 - csite%snowfac(ipa) + csite%snowfac(ipa) * rad_sfcw_nir )
      !------------------------------------------------------------------------------------!

      !----- Subtract off ground absorption to obtain the soil+sfcwater albedo. -----------!
      albedo_ground_par = albedo_ground_par - abs_ground_par
      albedo_ground_nir = albedo_ground_nir - abs_ground_nir
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Decide whether to call the radiation or not.  If there is no cohort, we can    !
      ! bypass it entirely.                                                                !
      !------------------------------------------------------------------------------------!
      if (cohort_count > 0) then

         !---------------------------------------------------------------------------------!
         !     First, let's solve the long-wave.  Here we must check which model to use,   !
         ! two-stream or multiple scattering.                                              !
         !---------------------------------------------------------------------------------!
         select case (icanrad)
         case (0) 
            !------------------------------------------------------------------------------!
            !    Original two-stream model.                                                !
            !------------------------------------------------------------------------------!
            call old_lw_two_stream(emissivity,T_surface,rlong,cohort_count                 &
                                  ,radscr(ibuff)%pft_array(1:cohort_count)                 &
                                  ,radscr(ibuff)%LAI_array(1:cohort_count)                 &
                                  ,radscr(ibuff)%WAI_array(1:cohort_count)                 &
                                  ,radscr(ibuff)%leaf_temp_array(1:cohort_count)           &
                                  ,radscr(ibuff)%wood_temp_array(1:cohort_count)           &
                                  ,radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count) &
                                  ,radscr(ibuff)%lw_v_array(1:cohort_count)                &
                                  ,downward_lw_below,upward_lw_below,upward_lw_above)
            !------------------------------------------------------------------------------!

         case (1)
            !------------------------------------------------------------------------------!
            !      Multiple-scatter model.  Here there is one important difference: we do  !
            ! NOT scale longwave radiation, and we do NOT split longwave radiation coming  !
            ! from the sky and coming from the surface, as they interact in the middle     !
            ! layers.  We save all radiation in the "incid" arrays, and scale them after   !
            ! the solution, so the code called after this step can be used the same way    !
            ! for both radiations.                                                         !
            !------------------------------------------------------------------------------!
            call lw_multiple_scatter(emissivity,T_surface,rlong,cohort_count               &
                            ,radscr(ibuff)%pft_array(1:cohort_count)                       &
                            ,radscr(ibuff)%LAI_array(1:cohort_count)                       &
                            ,radscr(ibuff)%WAI_array(1:cohort_count)                       &
                            ,radscr(ibuff)%CA_array(1:cohort_count)                        &
                            ,radscr(ibuff)%leaf_temp_array(1:cohort_count)                 &
                            ,radscr(ibuff)%wood_temp_array(1:cohort_count)                 &
                            ,radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count)       &
                            ,radscr(ibuff)%lw_v_array(1:cohort_count)                      &
                            ,downward_lw_below,upward_lw_below,upward_lw_above)
            !------------------------------------------------------------------------------!





         case (2) 
            !------------------------------------------------------------------------------!
            !    Updated two-stream model.                                                 !
            !------------------------------------------------------------------------------!
            call lw_two_stream(emissivity,T_surface,rlong,cohort_count                     &
                              ,radscr(ibuff)%pft_array(1:cohort_count)                     &
                              ,radscr(ibuff)%LAI_array(1:cohort_count)                     &
                              ,radscr(ibuff)%WAI_array(1:cohort_count)                     &
                              ,radscr(ibuff)%CA_array(1:cohort_count)                      &
                              ,radscr(ibuff)%leaf_temp_array(1:cohort_count)               &
                              ,radscr(ibuff)%wood_temp_array(1:cohort_count)               &
                              ,radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count)     &
                              ,radscr(ibuff)%lw_v_array(1:cohort_count)                    &
                              ,downward_lw_below,upward_lw_below,upward_lw_above)
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!


         !----- Upwelling long wave radiation at the top of the canopy. -------------------!
         csite%rlongup(ipa)      = upward_lw_above
         csite%rlong_albedo(ipa) = upward_lw_above / rlong
         !---------------------------------------------------------------------------------!

         !----- Net long wave absorption by either soil or sfcwater. ----------------------!
         surface_netabs_longwave = downward_lw_below - upward_lw_below
         !---------------------------------------------------------------------------------!

         !---------------------------------------------------------------------------------!
         !     Compute short wave only if it is daytime or at least twilight.              !
         !---------------------------------------------------------------------------------!
         if (twilight) then
            !------------------------------------------------------------------------------!
            !     Solve the short-wave.  Decide which canopy radiation method we are going !
            ! to use.  Unlike the long-wave, here we scale radiation the same way in both  !
            ! cases, because PAR and NIR are truly independent bands.                      !
            !------------------------------------------------------------------------------!
            select case (icanrad)
            case (0)
               !---------------------------------------------------------------------------!
               !    Two-stream model.                                                      !
               !---------------------------------------------------------------------------!
               call old_sw_two_stream(albedo_ground_par,albedo_ground_nir,cosaoi           &
                                  ,cohort_count                                            &
                                  ,radscr(ibuff)%pft_array(1:cohort_count)                 &
                                  ,radscr(ibuff)%LAI_array(1:cohort_count)                 &
                                  ,radscr(ibuff)%WAI_array(1:cohort_count)                 &
                                  ,radscr(ibuff)%CA_array(1:cohort_count)                  &
                                  ,radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count) &
                                  ,radscr(ibuff)%par_v_beam_array(1:cohort_count)          &
                                  ,radscr(ibuff)%par_v_diffuse_array(1:cohort_count)       &
                                  ,radscr(ibuff)%rshort_v_beam_array(1:cohort_count)       &
                                  ,radscr(ibuff)%rshort_v_diffuse_array(1:cohort_count)    &
                                  ,downward_par_below_beam                                 &
                                  ,downward_par_below_diffuse                              &
                                  ,upward_par_above_diffuse                                &
                                  ,downward_nir_below_beam                                 &
                                  ,downward_nir_below_diffuse                              &
                                  ,upward_nir_above_diffuse                                &
                                  ,radscr(ibuff)%par_level_beam(1:cohort_count)            &
                                  ,radscr(ibuff)%par_level_diffd(1:cohort_count)           &
                                  ,radscr(ibuff)%par_level_diffu(1:cohort_count)           &
                                  ,radscr(ibuff)%light_level_array(1:cohort_count)         &
                                  ,radscr(ibuff)%light_beam_level_array(1:cohort_count)    &
                                  ,radscr(ibuff)%light_diff_level_array(1:cohort_count))
               !---------------------------------------------------------------------------!

            case (1)
               !---------------------------------------------------------------------------!
               !      Multiple-scatter model.                                              !
               !---------------------------------------------------------------------------!

               call sw_multiple_scatter(albedo_ground_par,albedo_ground_nir,cosaoi         &
                                  ,cohort_count                                            &
                                  ,radscr(ibuff)%pft_array(1:cohort_count)                 &
                                  ,radscr(ibuff)%LAI_array(1:cohort_count)                 &
                                  ,radscr(ibuff)%WAI_array(1:cohort_count)                 &
                                  ,radscr(ibuff)%CA_array(1:cohort_count)                  &
                                  ,radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count) &
                                  ,radscr(ibuff)%par_v_beam_array(1:cohort_count)          &
                                  ,radscr(ibuff)%par_v_diffuse_array(1:cohort_count)       &
                                  ,radscr(ibuff)%rshort_v_beam_array(1:cohort_count)       &
                                  ,radscr(ibuff)%rshort_v_diffuse_array(1:cohort_count)    &
                                  ,downward_par_below_beam                                 &
                                  ,downward_par_below_diffuse                              &
                                  ,upward_par_above_diffuse                                &
                                  ,downward_nir_below_beam                                 &
                                  ,downward_nir_below_diffuse                              &
                                  ,upward_nir_above_diffuse                                &
                                  ,radscr(ibuff)%par_level_beam(1:cohort_count)            &
                                  ,radscr(ibuff)%par_level_diffd(1:cohort_count)           &
                                  ,radscr(ibuff)%par_level_diffu(1:cohort_count)           &
                                  ,radscr(ibuff)%light_level_array(1:cohort_count)         &
                                  ,radscr(ibuff)%light_beam_level_array(1:cohort_count)    &
                                  ,radscr(ibuff)%light_diff_level_array(1:cohort_count))

               !---------------------------------------------------------------------------!
            case (2)
               !---------------------------------------------------------------------------!
               call sw_two_stream(albedo_ground_par,albedo_ground_nir,cosaoi               &
                                 ,cohort_count                                             &
                                 ,radscr(ibuff)%pft_array(1:cohort_count)                  &
                                 ,radscr(ibuff)%LAI_array(1:cohort_count)                  &
                                 ,radscr(ibuff)%WAI_array(1:cohort_count)                  &
                                 ,radscr(ibuff)%CA_array(1:cohort_count)                   &
                                 ,radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count)  &
                                 ,radscr(ibuff)%par_v_beam_array(1:cohort_count)           &
                                 ,radscr(ibuff)%par_v_diffuse_array(1:cohort_count)        &
                                 ,radscr(ibuff)%rshort_v_beam_array(1:cohort_count)        &
                                 ,radscr(ibuff)%rshort_v_diffuse_array(1:cohort_count)     &
                                 ,downward_par_below_beam                                  &
                                 ,downward_par_below_diffuse                               &
                                 ,upward_par_above_diffuse                                 &
                                 ,downward_nir_below_beam                                  &
                                 ,downward_nir_below_diffuse                               &
                                 ,upward_nir_above_diffuse                                 &
                                 ,radscr(ibuff)%par_level_beam(1:cohort_count)             &
                                 ,radscr(ibuff)%par_level_diffd(1:cohort_count)            &
                                 ,radscr(ibuff)%par_level_diffu(1:cohort_count)            &
                                 ,radscr(ibuff)%light_level_array(1:cohort_count)          &
                                 ,radscr(ibuff)%light_beam_level_array(1:cohort_count)     &
                                 ,radscr(ibuff)%light_diff_level_array(1:cohort_count))
               ! THE LAST 5 ARRAYS SHOULD BE 1:cohort_count ?
               !---------------------------------------------------------------------------!
            end select
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !    Since there is no horizontal competition, assuming that the maximum       !
            ! possible PAR is just the PAR from the tallest resolvable cohort is good      !
            ! enough.                                                                      !
            !------------------------------------------------------------------------------!
            if (tuco_leaf /= 0) then
               ipft      = radscr(ibuff)%pft_array(tuco_leaf)
               wleaf_vis = sngloff( clumping_factor(ipft)                                  &
                                  * (1.d0 - leaf_scatter_vis(ipft))                        &
                                  * radscr(ibuff)%LAI_array(tuco_leaf)                     &
                                  / ( clumping_factor(ipft)                                &
                                    * (1.d0 - leaf_scatter_vis(ipft))                      &
                                    * radscr(ibuff)%LAI_array(tuco_leaf)                   &
                                    + (1.d0 - wood_scatter_vis(ipft))                      &
                                    * radscr(ibuff)%WAI_array(tuco_leaf) )                 &
                                  , tiny_offset)
               csite%par_l_beam_max(ipa)    = radscr(ibuff)%par_v_beam_array(tuco_leaf)    &
                                            * wleaf_vis   &
                                            / radscr(ibuff)%LAI_array(tuco_leaf)
               csite%par_l_diffuse_max(ipa) = radscr(ibuff)%par_v_diffuse_array(tuco_leaf) &
                                            * wleaf_vis   &
                                            / radscr(ibuff)%LAI_array(tuco_leaf)
            end if
            !------------------------------------------------------------------------------!



            !----- Above-canopy upwelling radiation. --------------------------------------!
            upward_rshort_above_diffuse = upward_par_above_diffuse                         &
                                        + upward_nir_above_diffuse
            !------------------------------------------------------------------------------!



            !----- Below-canopy downwelling radiation. ------------------------------------!
            downward_rshort_below_beam    = downward_par_below_beam                        &
                                          + downward_nir_below_beam
            downward_rshort_below_diffuse = downward_par_below_diffuse                     &
                                          + downward_nir_below_diffuse
            !------------------------------------------------------------------------------!


            !----- Soil+sfcwater+veg albedo (PAR, NIR, and Total Shortwave). --------------!
            csite%albedo_par (ipa) = upward_par_above_diffuse                              &
                                   / sngloff( par_beam_norm + par_diff_norm, tiny_offset )
            csite%albedo_nir (ipa) = upward_nir_above_diffuse                              &
                                   / sngloff( nir_beam_norm + nir_diff_norm, tiny_offset )
            csite%albedo     (ipa) = upward_rshort_above_diffuse
            !------------------------------------------------------------------------------!
         else


            !----- The code expects values for these, even when it is not daytime. --------!
            downward_par_below_beam         = par_beam_norm
            downward_par_below_diffuse      = par_diff_norm
            downward_nir_below_beam         = nir_beam_norm
            downward_nir_below_diffuse      = nir_diff_norm
            downward_rshort_below_beam      = par_beam_norm + nir_beam_norm
            downward_rshort_below_diffuse   = par_diff_norm + nir_diff_norm
            upward_par_above_diffuse        = albedo_ground_par * par_diff_norm
            upward_nir_above_diffuse        = albedo_ground_nir * nir_diff_norm
            upward_rshort_above_diffuse     = upward_par_above_diffuse                     &
                                            + upward_nir_above_diffuse
            csite%albedo_par          (ipa) = upward_par_above_diffuse                     &
                                            / ( par_beam_norm + par_diff_norm )
            csite%albedo_nir          (ipa) = upward_nir_above_diffuse                     &
                                            / ( nir_beam_norm + nir_diff_norm )
            csite%albedo              (ipa) = upward_rshort_above_diffuse
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Absorption rates of PAR, rshort, and rlong of the vegetation.  Here we      !
         ! check again whether we are solving big leaf or size- and age-structure.         !
         !---------------------------------------------------------------------------------!
         select case (ibigleaf)
         case (0)
            !------------------------------------------------------------------------------!
            !    Size- and age-structure.  We copy the results back to each cohort that is !
            ! resolvable.                                                                  !
            !------------------------------------------------------------------------------!
            il = 0
            do ico = cpatch%ncohorts,1,-1
               if (cpatch%leaf_resolvable(ico) .or. cpatch%wood_resolvable(ico)) then
                  il   = il + 1
                  ipft = radscr(ibuff)%pft_array(il)

                  !------------------------------------------------------------------------!
                  !      Find the weight for leaves and branchwood.  This is a weighted    !
                  ! average between the area and absorptance.  We must treat the visible   !
                  ! and near infrared separately.                                          !
                  !------------------------------------------------------------------------!
                  wleaf_vis = sngloff ( ( clumping_factor(ipft)                            &
                                        * (1.d0 - leaf_scatter_vis(ipft))                  &
                                        * radscr(ibuff)%LAI_array(il) )                    &
                                      / ( clumping_factor(ipft)                            &
                                        * (1.d0 - leaf_scatter_vis(ipft))                  &
                                        * radscr(ibuff)%LAI_array(il)                      &
                                        + (1.d0 - wood_scatter_vis(ipft))                  &
                                        * radscr(ibuff)%WAI_array(il) ), tiny_offset )
                  wleaf_nir = sngloff ( ( clumping_factor(ipft)                            &
                                        * (1.d0 - leaf_scatter_nir(ipft))                  &
                                        * radscr(ibuff)%LAI_array(il) )                    &
                                      / ( clumping_factor(ipft)                            &
                                        * (1.d0 - leaf_scatter_nir(ipft))                  &
                                        * radscr(ibuff)%LAI_array(il)                      &
                                        + (1.d0 - wood_scatter_nir(ipft))                  &
                                        * radscr(ibuff)%WAI_array(il) ), tiny_offset  )
                  wleaf_tir = sngloff( ( leaf_emiss_tir(ipft)                              &
                                       * radscr(ibuff)%LAI_array(il) )                     &
                                     / ( leaf_emiss_tir(ipft)                              &
                                       * radscr(ibuff)%LAI_array(il)                       &
                                       + wood_emiss_tir(ipft)                              &
                                       * radscr(ibuff)%WAI_array(il) )                     &
                                     , tiny_offset )
                  wwood_vis = 1. - wleaf_vis
                  wwood_nir = 1. - wleaf_nir
                  wwood_tir = 1. - wleaf_tir
                  !------------------------------------------------------------------------!

                  !----- Find the NIR absorption, so we average things properly. ----------!
                  nir_v_beam    = radscr(ibuff)%rshort_v_beam_array    (il) -              &
                                  radscr(ibuff)%par_v_beam_array    (il)
                  nir_v_diffuse = radscr(ibuff)%rshort_v_diffuse_array (il) -              &
                                  radscr(ibuff)%par_v_diffuse_array (il)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Light profile: we do not split into leaves and branches.            !
                  !------------------------------------------------------------------------!
                  cpatch%rad_profile    (:,ico) = radscr(ibuff)%radprof_array(:,il)
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !    Split the layer radiation between leaf and branchwood.              !
                  !------------------------------------------------------------------------!
                  !------ Visible (PAR), only leaves need this (photsynthesis model). -----!
                  cpatch%par_l_beam       (ico) = radscr(ibuff)%par_v_beam_array(il)       &
                                                * wleaf_vis
                  cpatch%par_l_diffuse    (ico) = radscr(ibuff)%par_v_diffuse_array(il)    &
                                                * wleaf_vis
                  !------ Total short wave radiation (PAR+NIR). ---------------------------!
                  cpatch%rshort_l_beam    (ico) = radscr(ibuff)%par_v_beam_array    (il)   &
                                                * wleaf_vis + nir_v_beam * wleaf_nir
                  cpatch%rshort_l_diffuse (ico) = radscr(ibuff)%par_v_diffuse_array (il)   &
                                                * wleaf_vis + nir_v_diffuse * wleaf_nir
                  cpatch%rshort_w_beam    (ico) = radscr(ibuff)%par_v_beam_array    (il)   &
                                                * wwood_vis + nir_v_beam * wwood_nir
                  cpatch%rshort_w_diffuse (ico) = radscr(ibuff)%par_v_diffuse_array (il)   &
                                                * wwood_vis + nir_v_diffuse * wwood_nir
                  !----- Thermal infra-red (long wave). -----------------------------------!
                  cpatch%rlong_l          (ico) = radscr(ibuff)%lw_v_array(il) * wleaf_tir
                  cpatch%rlong_w          (ico) = radscr(ibuff)%lw_v_array(il) * wwood_tir
                  !------------------------------------------------------------------------!

                  !----- Save the light levels. -------------------------------------------!
                  cpatch%light_level(ico)      =                                           &
                             sngloff(radscr(ibuff)%light_level_array(il)     ,tiny_offset)
                  cpatch%light_level_beam(ico) =                                           &
                             sngloff(radscr(ibuff)%light_beam_level_array(il),tiny_offset)
                  cpatch%light_level_diff(ico) =                                           &
                             sngloff(radscr(ibuff)%light_diff_level_array(il),tiny_offset)
                  cpatch%par_level_beam  (ico) =                                           &
                             sngloff(radscr(ibuff)%par_level_beam(il)        ,tiny_offset)
                  cpatch%par_level_diffu (ico) =                                           &
                             sngloff(radscr(ibuff)%par_level_diffu(il)       ,tiny_offset)
                  cpatch%par_level_diffd (ico) =                                           &
                             sngloff(radscr(ibuff)%par_level_diffd(il)       ,tiny_offset)

                  !------------------------------------------------------------------------!
               end if
               !---------------------------------------------------------------------------!
            end do
            !------------------------------------------------------------------------------!


         case (1)
            !------------------------------------------------------------------------------!
            !     Big leaf solver.  Radiation flux is an extensive variable since its      !
            ! units are J/m2/s.  Therefore we just need to add the fluxes back.  For the   !
            ! light levels, we cheat and assign the value in the middle.                   !
            !------------------------------------------------------------------------------!
            if (cpatch%leaf_resolvable(1) .or. cpatch%wood_resolvable(1)) then



               !---------------------------------------------------------------------------!
               !    Light profile: we do not split into leaves and branches.  Since this   !
               ! is a big leaf, we only save the values beneath the stack.                 !
               !---------------------------------------------------------------------------!
               cpatch%rad_profile       (:,1) = radscr(ibuff)%radprof_array(:,1)
               !---------------------------------------------------------------------------!


               do il=1,cohort_count
                  ipft = radscr(ibuff)%pft_array(il)

                  !------------------------------------------------------------------------!
                  !      Find the weight for leaves and branchwood.  This is a weighted    !
                  ! average between the area and absorptance.  We must treat the visible   !
                  ! and near infrared separately.                                          !
                  !------------------------------------------------------------------------!
                  wleaf_vis = sngloff ( ( clumping_factor(ipft)                            &
                                        * (1.d0 - leaf_scatter_vis(ipft))                  &
                                        * radscr(ibuff)%LAI_array(il) )                    &
                                      / ( clumping_factor(ipft)                            &
                                        * (1.d0 - leaf_scatter_vis(ipft))                  &
                                        * radscr(ibuff)%LAI_array(il)                      &
                                        + (1.d0 - wood_scatter_vis(ipft))                  &
                                        * radscr(ibuff)%WAI_array(il) ), tiny_offset )
                  wleaf_nir = sngloff ( ( clumping_factor(ipft)                            &
                                        * (1.d0 - leaf_scatter_nir(ipft))                  &
                                        * radscr(ibuff)%LAI_array(il) )                    &
                                      / ( clumping_factor(ipft)                            &
                                        * (1.d0 - leaf_scatter_nir(ipft))                  &
                                        * radscr(ibuff)%LAI_array(il)                      &
                                        + (1.d0 - wood_scatter_nir(ipft))                  &
                                        * radscr(ibuff)%WAI_array(il) ), tiny_offset  )
                  wleaf_tir = sngloff( ( leaf_emiss_tir(ipft)                              &
                                       * radscr(ibuff)%LAI_array(il) )                     &
                                     / ( leaf_emiss_tir(ipft)                              &
                                       * radscr(ibuff)%LAI_array(il)                       &
                                       + wood_emiss_tir(ipft)                              &
                                       * radscr(ibuff)%WAI_array(il) )                     &
                                     , tiny_offset    )
                  wwood_vis = 1. - wleaf_vis
                  wwood_nir = 1. - wleaf_nir
                  wwood_tir = 1. - wleaf_tir
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !     Find the near infrared absorption, so we average things properly.  !
                  !------------------------------------------------------------------------!
                  nir_v_beam    = radscr(ibuff)%rshort_v_beam_array   (il)                 &
                                - radscr(ibuff)%par_v_beam_array      (il)
                  nir_v_diffuse = radscr(ibuff)%rshort_v_diffuse_array(il)                 &
                                - radscr(ibuff)%par_v_diffuse_array   (il)
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !    Split the layer radiation between leaf and branchwood.              !
                  !------------------------------------------------------------------------!
                  !------ Visible (PAR), only leaves need this (photsynthesis model). -----!
                  cpatch%par_l_beam       (1) = cpatch%par_l_beam(1)                       &
                                              + radscr(ibuff)%par_v_beam_array(il)         &
                                              * wleaf_vis
                  cpatch%par_l_diffuse    (1) = cpatch%par_l_diffuse(1)                    &
                                              + radscr(ibuff)%par_v_diffuse_array(il)      &
                                              * wleaf_vis
                  !------ Total short wave radiation (PAR+NIR). ---------------------------!
                  cpatch%rshort_l_beam    (1) = cpatch%rshort_l_beam(1)                    &
                                              + radscr(ibuff)%par_v_beam_array(il)         &
                                              * wleaf_vis + nir_v_beam* wleaf_nir
                  cpatch%rshort_l_diffuse (1) = cpatch%rshort_l_diffuse (1)                &
                                              + radscr(ibuff)%par_v_diffuse_array(il)      &
                                              * wleaf_vis + nir_v_diffuse * wleaf_nir
                  cpatch%rshort_w_beam    (1) = cpatch%rshort_w_beam    (1)                &
                                              + radscr(ibuff)%par_v_beam_array(il)         &
                                              * wwood_vis + nir_v_beam * wwood_nir
                  cpatch%rshort_w_diffuse (1) = cpatch%rshort_w_diffuse (1)                &
                                              + radscr(ibuff)%par_v_diffuse_array(il)      &
                                              * wwood_vis + nir_v_diffuse * wwood_nir
                  !----- Thermal infra-red (long wave). -----------------------------------!
                  cpatch%rlong_l          (1) = cpatch%rlong_l (1)                         &
                                              + radscr(ibuff)%lw_v_array    (il)           &
                                              * wleaf_tir
                  cpatch%rlong_w          (1) = cpatch%rlong_w (1)                         &
                                              + radscr(ibuff)%lw_v_array    (il)           &
                                              * wwood_tir
                  !------------------------------------------------------------------------!
               end do
               !---------------------------------------------------------------------------!



               !----- Save the light levels as the median level. --------------------------!
               il = ceiling(real(cohort_count)/2.0)
               cpatch%light_level      (1) =                                               &
                             sngloff(radscr(ibuff)%light_level_array     (il),tiny_offset )
               cpatch%light_level_beam (1) =                                               &
                             sngloff(radscr(ibuff)%light_beam_level_array(il),tiny_offset )
               cpatch%light_level_diff (1) =                                               &
                             sngloff(radscr(ibuff)%light_diff_level_array(il),tiny_offset )
               cpatch%par_level_beam   (1) =                                               &
                             sngloff(radscr(ibuff)%par_level_beam(il)        ,tiny_offset )
               cpatch%par_level_diffu  (1) =                                               &
                             sngloff(radscr(ibuff)%par_level_diffu(il)       ,tiny_offset )
               cpatch%par_level_diffd  (1) =                                               &
                             sngloff(radscr(ibuff)%par_level_diffd(il)       ,tiny_offset )
               !---------------------------------------------------------------------------!
            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!

      else

         !----- This is the case where there is no vegetation. ----------------------------!
         downward_par_below_beam         = par_beam_norm
         downward_par_below_diffuse      = par_diff_norm
         downward_nir_below_beam         = nir_beam_norm
         downward_nir_below_diffuse      = nir_diff_norm
         downward_rshort_below_beam      = par_beam_norm + nir_beam_norm
         downward_rshort_below_diffuse   = par_diff_norm + nir_diff_norm

         upward_par_above_diffuse        = albedo_ground_par * par_diff_norm
         upward_nir_above_diffuse        = albedo_ground_nir * nir_diff_norm
         upward_rshort_above_diffuse     = upward_par_above_diffuse                        &
                                         + upward_nir_above_diffuse

         surface_netabs_longwave         = emissivity * (rlong - stefan * T_surface**4 )

         csite%albedo_par          (ipa) = upward_par_above_diffuse                        &
                                         / ( par_beam_norm + par_diff_norm )
         csite%albedo_nir          (ipa) = upward_nir_above_diffuse                        &
                                         / ( nir_beam_norm + nir_diff_norm )
         csite%albedo              (ipa) = upward_rshort_above_diffuse
         csite%rlongup             (ipa) = rlong - surface_netabs_longwave
         csite%rlong_albedo        (ipa) = csite%rlongup(ipa) / rlong
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!


      !----- Absorption rate of short wave by the soil. -----------------------------------!
      csite%rshort_g_beam   (ipa) = downward_par_below_beam    * abs_ground_par            &
                                  + downward_nir_below_beam    * abs_ground_nir
      csite%rshort_g_diffuse(ipa) = downward_par_below_diffuse * abs_ground_par            &
                                  + downward_nir_below_diffuse * abs_ground_nir
      csite%par_g_beam      (ipa) = downward_par_below_beam    * abs_ground_par
      csite%par_g_diffuse   (ipa) = downward_par_below_diffuse * abs_ground_par
      csite%parup           (ipa) = upward_par_above_diffuse
      csite%nirup           (ipa) = upward_nir_above_diffuse

      !----- Incident radiation at the ground. --------------------------------------------!
      csite%par_b_beam      (ipa) = downward_par_below_beam
      csite%nir_b_beam      (ipa) = downward_nir_below_beam
      csite%par_b_diffuse   (ipa) = downward_par_below_diffuse
      csite%nir_b_diffuse   (ipa) = downward_nir_below_diffuse

      !----- Absorption rate of short wave by the surface water. --------------------------!
      do k=1,csite%nlev_sfcwater(ipa)
         csite%rshort_s_beam   (k,ipa) = downward_par_below_beam    * abs_sfcw_par(k)      &
                                       + downward_nir_below_beam    * abs_sfcw_nir(k)
         csite%rshort_s_diffuse(k,ipa) = downward_par_below_diffuse * abs_sfcw_par(k)      &
                                       + downward_nir_below_beam    * abs_sfcw_nir(k)
         csite%par_s_beam      (k,ipa) = downward_par_below_beam    * abs_sfcw_par(k)
         csite%par_s_diffuse   (k,ipa) = downward_par_below_diffuse * abs_sfcw_par(k)
      end do

      !----- Long wave absorption rate at the surface. ------------------------------------!
      if (csite%nlev_sfcwater(ipa) == 0) then
         csite%rlong_s (ipa) = 0.0
         csite%rlong_g (ipa) = surface_netabs_longwave
      else
         csite%rlong_s (ipa) = surface_netabs_longwave
         csite%rlong_g (ipa) = 0.0
      end if
      !------------------------------------------------------------------------------------!


      return
   end subroutine sfcrad_ed
   !=======================================================================================!
   !=======================================================================================!







   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   subroutine scale_ed_radiation(rshort,rlong,nighttime,csite,ipa)

      use ed_state_vars        , only : sitetype             & ! intent(in)
                                      , patchtype            ! ! intent(in)
      use ed_misc_coms         , only : writing_long         & ! intent(in)
                                      , radfrq               & ! intent(in)
                                      , radfrq_o_frqsum      ! ! intent(in)
      !$ use omp_lib
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)  , target     :: csite
      real            , intent(in) :: rshort
      real            , intent(in) :: rlong
      logical         , intent(in) :: nighttime
      integer         , intent(in) :: ipa
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype) , pointer    :: cpatch
      integer                      :: ico
      integer                      :: k
      !----- This should be false unless you really want to turn off radiation. -----------!
      logical         , parameter  :: skip_rad = .false.
      !----- External functions. ----------------------------------------------------------!
      real            , external   :: sngloff
      !------------------------------------------------------------------------------------!

      !----- Link to current patch. -------------------------------------------------------!
      cpatch => csite%patch(ipa)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     This block skips radiation.  Obviously this only makes sense for very          !
      ! theoretical tests, because plants kind of like light...                            !
      !------------------------------------------------------------------------------------!
      if (skip_rad) then
         skip_cohortloop: do ico = 1, cpatch%ncohorts
            if (cpatch%leaf_resolvable(ico) .or. cpatch%wood_resolvable(ico)) then
               cpatch%par_l_beam         (ico) = 0.0
               cpatch%par_l_diffuse      (ico) = 0.0
               cpatch%par_l              (ico) = 0.0
               cpatch%rshort_l_beam      (ico) = 0.0
               cpatch%rshort_l_diffuse   (ico) = 0.0
               cpatch%rshort_l           (ico) = 0.0
               cpatch%rlong_l            (ico) = 0.0
               cpatch%rshort_w_beam      (ico) = 0.0
               cpatch%rshort_w_diffuse   (ico) = 0.0
               cpatch%rshort_w           (ico) = 0.0
               cpatch%rlong_w            (ico) = 0.0
               cpatch%light_level        (ico) = 0.0
               cpatch%light_level_diff   (ico) = 0.0
               cpatch%light_level_beam   (ico) = 0.0
               cpatch%rad_profile      (:,ico) = 0.0
            end if
         end do skip_cohortloop

         csite%rshort_g_beam   (ipa) = 0.
         csite%rshort_g_diffuse(ipa) = 0.
         csite%rshort_g        (ipa) = 0.
         csite%par_g_beam      (ipa) = 0.
         csite%par_g_diffuse   (ipa) = 0.
         csite%par_g           (ipa) = 0.
         csite%par_b_beam      (ipa) = 0.
         csite%par_b_diffuse   (ipa) = 0.
         csite%par_b           (ipa) = 0.
         csite%nir_b_beam      (ipa) = 0.
         csite%nir_b_diffuse   (ipa) = 0.
         csite%nir_b           (ipa) = 0.
         csite%parup           (ipa) = 0.
         csite%nirup           (ipa) = 0.
         csite%rshortup        (ipa) = 0.
         csite%rnet            (ipa) = 0.
         !----- Absorption rate of short wave by the surface water. -----------------------!
         do k=1,csite%nlev_sfcwater(ipa)
            csite%rshort_s_beam   (k,ipa) = 0.
            csite%rshort_s_diffuse(k,ipa) = 0.
            csite%rshort_s        (k,ipa) = 0.
            csite%par_s_beam      (k,ipa) = 0.
            csite%par_s_diffuse   (k,ipa) = 0.
            csite%par_s           (k,ipa) = 0.
         end do
         !---------------------------------------------------------------------------------!

         csite%rlong_s(ipa)       = 0.
         csite%rlong_g(ipa)       = 0.
         return
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     For normal runs, we add the scale to radiation, and also integrate averages.   !
      !------------------------------------------------------------------------------------!
      !----- Cohort-level variables. ------------------------------------------------------!
      cohortloop: do ico = 1,cpatch%ncohorts

         if (cpatch%leaf_resolvable(ico) .or. cpatch%wood_resolvable(ico)) then

            cpatch%par_l_beam(ico)       = cpatch%par_l_beam(ico)    * rshort
            cpatch%par_l_diffuse(ico)    = cpatch%par_l_diffuse(ico) * rshort
            cpatch%par_l(ico)            = cpatch%par_l_beam(ico)                          &
                                         + cpatch%par_l_diffuse(ico)

            cpatch%rshort_l_beam(ico)    = cpatch%rshort_l_beam(ico)    * rshort
            cpatch%rshort_l_diffuse(ico) = cpatch%rshort_l_diffuse(ico) * rshort
            cpatch%rshort_l(ico)         = cpatch%rshort_l_beam(ico)                       &
                                         + cpatch%rshort_l_diffuse(ico)

            cpatch%rshort_w_beam(ico)    = cpatch%rshort_w_beam(ico)    * rshort
            cpatch%rshort_w_diffuse(ico) = cpatch%rshort_w_diffuse(ico) * rshort
            cpatch%rshort_w(ico)         = cpatch%rshort_w_beam(ico)                       &
                                         + cpatch%rshort_w_diffuse(ico)


            !----- Only short wave radiation must be scaled... ----------------------------!
            cpatch%rad_profile(1:8,ico)  = cpatch%rad_profile(1:8,ico) * rshort
            !------------------------------------------------------------------------------!
         end if
         !---------------------------------------------------------------------------------!
      end do cohortloop
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      !     Patch-level variables.                                                         !
      !------------------------------------------------------------------------------------!
      csite%par_l_beam_max    (ipa)      = csite%par_l_beam_max    (ipa) * rshort
      csite%par_l_diffuse_max (ipa)      = csite%par_l_diffuse_max (ipa) * rshort
      csite%par_l_max         (ipa)      = csite%par_l_beam_max    (ipa)                   &
                                         + csite%par_l_diffuse_max (ipa)

      csite%rshort_g_beam     (ipa)      = csite%rshort_g_beam     (ipa) * rshort
      csite%rshort_g_diffuse  (ipa)      = csite%rshort_g_diffuse  (ipa) * rshort
      csite%rshort_g          (ipa)      = csite%rshort_g_beam     (ipa)                   &
                                         + csite%rshort_g_diffuse  (ipa)

      csite%par_g_beam        (ipa)      = csite%par_g_beam        (ipa) * rshort
      csite%par_g_diffuse     (ipa)      = csite%par_g_diffuse     (ipa) * rshort
      csite%par_g             (ipa)      = csite%par_g_beam        (ipa)                   &
                                         + csite%par_g_diffuse     (ipa)

      csite%par_b_beam        (ipa)      = csite%par_b_beam        (ipa) * rshort
      csite%par_b_diffuse     (ipa)      = csite%par_b_diffuse     (ipa) * rshort
      csite%par_b             (ipa)      = csite%par_b_beam        (ipa)                   &
                                         + csite%par_b_diffuse     (ipa)
      csite%nir_b_beam        (ipa)      = csite%nir_b_beam        (ipa) * rshort
      csite%nir_b_diffuse     (ipa)      = csite%nir_b_diffuse     (ipa) * rshort
      csite%nir_b             (ipa)      = csite%nir_b_beam        (ipa)                   &
                                         + csite%nir_b_diffuse     (ipa)
      csite%parup             (ipa)      = csite%parup             (ipa) * rshort
      csite%nirup             (ipa)      = csite%nirup             (ipa) * rshort
      csite%rshortup          (ipa)      = csite%parup             (ipa)                   &
                                         + csite%nirup             (ipa)
      csite%rnet              (ipa)      = rshort + rlong                                  &
                                         - csite%rshortup          (ipa)                   &
                                         - csite%rlongup           (ipa)
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !     Patch-level variables, but the absorption rate by each temporary pounding/snow !
      ! layer.                                                                             !
      !------------------------------------------------------------------------------------!
      do k=1,csite%nlev_sfcwater(ipa)
         csite%rshort_s_beam   (k,ipa) = csite%rshort_s_beam   (k,ipa) * rshort
         csite%rshort_s_diffuse(k,ipa) = csite%rshort_s_diffuse(k,ipa) * rshort
         csite%rshort_s        (k,ipa) = csite%rshort_s_beam   (k,ipa)                     &
                                       + csite%rshort_s_diffuse(k,ipa)
         csite%par_s_beam      (k,ipa) = csite%par_s_beam      (k,ipa) * rshort
         csite%par_s_diffuse   (k,ipa) = csite%par_s_diffuse   (k,ipa) * rshort
         csite%par_s           (k,ipa) = csite%par_s_beam      (k,ipa)                     &
                                       + csite%par_s_diffuse   (k,ipa)
      end do
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Integrate the mean radiation fluxes.  The average fluxes are normalised in    !
      ! average_utils.f90, at sub-routine normalize_averaged_vars.                         !
      !------------------------------------------------------------------------------------!
      !----- Cohort-level variables. ------------------------------------------------------!
      mean_cohortloop: do ico=1,cpatch%ncohorts
         cpatch%fmean_par_l             (ico) = cpatch%fmean_par_l             (ico)       &
                                              + cpatch%par_l                   (ico)       &
                                              * radfrq_o_frqsum
         cpatch%fmean_par_l_beam        (ico) = cpatch%fmean_par_l_beam        (ico)       &
                                              + cpatch%par_l_beam              (ico)       &
                                              * radfrq_o_frqsum
         cpatch%fmean_par_l_diff        (ico) = cpatch%fmean_par_l_diff        (ico)       &
                                              + cpatch%par_l_diffuse           (ico)       &
                                              * radfrq_o_frqsum
         cpatch%fmean_rshort_l          (ico) = cpatch%fmean_rshort_l          (ico)       &
                                              + cpatch%rshort_l                (ico)       &
                                              * radfrq_o_frqsum
         cpatch%fmean_rlong_l           (ico) = cpatch%fmean_rlong_l           (ico)       &
                                              + cpatch%rlong_l                 (ico)       &
                                              * radfrq_o_frqsum
         cpatch%fmean_rshort_w          (ico) = cpatch%fmean_rshort_w          (ico)       &
                                              + cpatch%rshort_w                (ico)       &
                                              * radfrq_o_frqsum
         cpatch%fmean_rlong_w           (ico) = cpatch%fmean_rlong_w           (ico)       &
                                              + cpatch%rlong_w                 (ico)       &
                                              * radfrq_o_frqsum

         cpatch%fmean_par_level_beam    (ico) = cpatch%fmean_par_level_beam    (ico)       &
                                              + cpatch%par_level_beam(ico)                 &
                                              * radfrq_o_frqsum

         cpatch%fmean_par_level_diffd   (ico) = cpatch%fmean_par_level_diffd   (ico)       &
                                              + cpatch%par_level_diffd(ico)                &
                                              * radfrq_o_frqsum

         cpatch%fmean_par_level_diffu   (ico) = cpatch%fmean_par_level_diffu   (ico)       &
                                              + cpatch%par_level_diffu(ico)                &
                                              * radfrq_o_frqsum

         cpatch%fmean_light_level       (ico) = cpatch%fmean_light_level       (ico)       &
                                              + cpatch%light_level             (ico)       &
                                              * radfrq_o_frqsum
         cpatch%fmean_light_level_beam  (ico) = cpatch%fmean_light_level_beam  (ico)       &
                                              + cpatch%light_level_beam        (ico)       &
                                              * radfrq_o_frqsum
         cpatch%fmean_light_level_diff  (ico) = cpatch%fmean_light_level_diff  (ico)       &
                                              + cpatch%light_level_diff        (ico)       &
                                              * radfrq_o_frqsum
         cpatch%fmean_rad_profile     (:,ico) = cpatch%fmean_rad_profile     (:,ico)       &
                                              + cpatch%rad_profile           (:,ico)       &
                                              * radfrq_o_frqsum
         !----- Light level is integrated only when there is some radiation. --------------!
         if (.not. nighttime .and. writing_long) then
            cpatch%dmean_light_level     (ico) = cpatch%dmean_light_level     (ico)        &
                                               + cpatch%light_level           (ico)        &
                                               * radfrq
            cpatch%dmean_light_level_beam(ico) = cpatch%dmean_light_level_beam(ico)        &
                                               + cpatch%light_level_beam      (ico)        &
                                               * radfrq
            cpatch%dmean_light_level_diff(ico) = cpatch%dmean_light_level_diff(ico)        &
                                               + cpatch%light_level_diff      (ico)        &
                                               * radfrq
         end if
         !---------------------------------------------------------------------------------!
      end do mean_cohortloop
      !----- Patch-level variables. -------------------------------------------------------!
      csite%fmean_rshort_gnd      (ipa) = csite%fmean_rshort_gnd   (ipa)                   &
                                        + csite%rshort_g           (ipa)                   &
                                        * radfrq_o_frqsum
      csite%fmean_par_gnd         (ipa) = csite%fmean_par_gnd      (ipa)                   &
                                        + csite%par_g              (ipa)                   &
                                        * radfrq_o_frqsum
      csite%fmean_rlong_gnd       (ipa) = csite%fmean_rlong_gnd    (ipa)                   &
                                        + csite%rlong_g            (ipa)                   &
                                        * radfrq_o_frqsum
      csite%fmean_parup           (ipa) = csite%fmean_parup        (ipa)                   &
                                        + csite%parup              (ipa)                   &
                                        * radfrq_o_frqsum
      csite%fmean_nirup           (ipa) = csite%fmean_nirup        (ipa)                   &
                                        + csite%nirup              (ipa)                   &
                                        * radfrq_o_frqsum
      csite%fmean_rshortup        (ipa) = csite%fmean_rshortup     (ipa)                   &
                                        + csite%rshortup           (ipa)                   &
                                        * radfrq_o_frqsum
      csite%fmean_rlongup         (ipa) = csite%fmean_rlongup      (ipa)                   &
                                        + csite%rlongup            (ipa)                   &
                                        * radfrq_o_frqsum
      csite%fmean_rnet            (ipa) = csite%fmean_rnet         (ipa)                   &
                                        + csite%rnet               (ipa)                   &
                                        * radfrq_o_frqsum
      csite%fmean_albedo          (ipa) = csite%fmean_albedo       (ipa)                   &
                                        + csite%albedo             (ipa)                   &
                                        * radfrq_o_frqsum
      csite%fmean_albedo_par      (ipa) = csite%fmean_albedo_par   (ipa)                   &
                                        + csite%albedo_par         (ipa)                   &
                                        * radfrq_o_frqsum
      csite%fmean_albedo_nir      (ipa) = csite%fmean_albedo_nir   (ipa)                   &
                                        + csite%albedo_nir         (ipa)                   &
                                        * radfrq_o_frqsum
      csite%fmean_rlong_albedo    (ipa) = csite%fmean_rlong_albedo (ipa)                   &
                                        + csite%rlong_albedo       (ipa)                   &
                                        * radfrq_o_frqsum
      !----- Daily mean of albedo is integrated only when there is some radiation. --------!
      if (.not. nighttime .and. writing_long) then
         csite%dmean_albedo        (ipa)  = csite%dmean_albedo        (ipa)                &
                                          + csite%albedo              (ipa)                &
                                          * radfrq
         csite%dmean_albedo_par    (ipa)  = csite%dmean_albedo_par    (ipa)                &
                                          + csite%albedo_par          (ipa)                &
                                          * radfrq
         csite%dmean_albedo_nir    (ipa)  = csite%dmean_albedo_nir    (ipa)                &
                                          + csite%albedo_nir          (ipa)                &
                                          * radfrq
      end if
      !------------------------------------------------------------------------------------!


      return
   end subroutine scale_ed_radiation
   !=======================================================================================!
   !=======================================================================================!
end module radiate_driver
