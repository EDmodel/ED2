module radiate_driver_module
   contains
!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the two-stream radiation scheme.  This is called every  !
! step, but not every sub-step.                                                            !
!------------------------------------------------------------------------------------------!
subroutine radiate_driver(cgrid)
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
   implicit none
   !----- Argument. -----------------------------------------------------------------------!
   type(edtype)     , target   :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer  :: cpoly
   type(sitetype)   , pointer  :: csite
   type(patchtype)  , pointer  :: cpatch
   real                        :: rshort_tot
   integer                     :: maxcohort
   integer                     :: ipy
   integer                     :: isi
   integer                     :: ipa
   integer                     :: tuco
   logical                     :: daytime
   logical                     :: twilight
   real                        :: hrangl
   real                        :: sloperad
   real                        :: aspectrad
   real                        :: sum_norm
   !----- External functions. -------------------------------------------------------------!
   real             , external :: ed_zen
   !---------------------------------------------------------------------------------------!


   !----- Check whether it is time to update radiative fluxes and heating rates -----------!
   if (mod(current_time%time + .001,radfrq) < dtlsm) then

      !----- Loop over polygons and sites. ------------------------------------------------!

      polyloop: do ipy = 1,cgrid%npolygons

         !----- Find the solar zenith angle [cosz] -----------------------------------------!
         cgrid%cosz(ipy) = ed_zen(cgrid%lon(ipy),cgrid%lat(ipy),current_time)
         !---------------------------------------------------------------------------------!

         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites

            csite => cpoly%site(isi)

            !------------------------------------------------------------------------------!
            !     Update angle of incidence.                                               !
            !------------------------------------------------------------------------------!
            hrangl    = 15. * pio180                                                       &
                      * (mod(current_time%time + cgrid%lon(ipy) / 15. + 24., 24.) - 12.)
            sloperad  = cpoly%slope(isi)  * pio180
            aspectrad = cpoly%aspect(isi) * pio180
            call angle_of_incid(cpoly%cosaoi(isi),cgrid%cosz(ipy),hrangl                   &
                               ,sloperad,aspectrad)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !    Find the two logicals, that will tell which part of the day we are at     !
            ! least.                                                                       !
            !------------------------------------------------------------------------------!
            daytime  = cpoly%cosaoi(isi) > cosz_min .and.                                  &
                       cpoly%met(isi)%rshort > rshort_twilight_min
!            twilight = cpoly%met(isi)%rshort > rshort_twilight_min
            twilight = rshort_twilight_min > 0.
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !      Update the daylight length and nighttime flag.                          !
            !------------------------------------------------------------------------------!
            cpoly%nighttime(isi)              = .not. twilight
            if (twilight) cpoly%daylight(isi) = cpoly%daylight(isi) + radfrq
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            !      In case the angle of incidence is too high (i.e., its cosine is too     !
            ! close to zero), eliminate all direct radiation.   This is different from the !
            ! cosine of the zenith angle because mountains can hide the sun even when it   !
            ! is still above the horizon.                                                  !
            !------------------------------------------------------------------------------!
            if (daytime) then
               rshort_tot    = cpoly%met(isi)%rshort
               par_beam_norm = max( 1.d-5 , dble(cpoly%met(isi)%par_beam   )               &
                                          / dble(cpoly%met(isi)%rshort     ) )
               par_diff_norm = max( 1.d-5 , dble(cpoly%met(isi)%par_diffuse)               &
                                          / dble(cpoly%met(isi)%rshort     ) )
               nir_beam_norm = max( 1.d-5 , dble(cpoly%met(isi)%nir_beam   )               &
                                          / dble(cpoly%met(isi)%rshort     ) )
               nir_diff_norm = max( 1.d-5 , dble(cpoly%met(isi)%nir_diffuse)               &
                                          / dble(cpoly%met(isi)%rshort     ) )
               sum_norm      = par_beam_norm + par_diff_norm                               &
                             + nir_beam_norm + nir_diff_norm
            elseif (twilight) then
               rshort_tot    = cpoly%met(isi)%rshort_diffuse
               par_beam_norm = 1.d-5
               par_diff_norm = 1.d-5
               nir_beam_norm = max(1.d-5, dble(cpoly%met(isi)%nir_beam   )                 &
                                        / dble(cpoly%met(isi)%rshort     ) )
               nir_diff_norm = max(1.d-5, dble(cpoly%met(isi)%nir_diffuse)                 &
                                        / dble(cpoly%met(isi)%rshort     ) )
               sum_norm      = par_beam_norm + par_diff_norm                               &
                             + nir_beam_norm + nir_diff_norm
            else 
               !---------------------------------------------------------------------------!
               !     Night-time, nothing will happen, fill split equally to the 4          !
               ! components.                                                               !
               !---------------------------------------------------------------------------!
               rshort_tot    = 0.0
               par_beam_norm = 2.5d-1
               par_diff_norm = 2.5d-1
               nir_beam_norm = 2.5d-1
               nir_diff_norm = 2.5d-1
               sum_norm      = 1.d0
            end if
            !------------------------------------------------------------------------------!
            !     Because we must tweak the radiation so none of the terms are zero, we    !
            ! must correct the normalised radiation variables so they add up to one.       !
            !------------------------------------------------------------------------------!
            par_beam_norm = par_beam_norm / sum_norm
            par_diff_norm = par_diff_norm / sum_norm
            nir_beam_norm = nir_beam_norm / sum_norm
            nir_diff_norm = nir_diff_norm / sum_norm
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !    Loop over subgrid-scale patches.  These routines can be done as arrays.   !
            !------------------------------------------------------------------------------!
            maxcohort = 1
            do ipa = 1,csite%npatches
               cpatch=>csite%patch(ipa)
               if ( cpatch%ncohorts>maxcohort ) maxcohort = cpatch%ncohorts
            end do
            !------------------------------------------------------------------------------!

            !----- Get unnormalized radiative transfer information. -----------------------!
            call sfcrad_ed(cgrid%cosz(ipy),cpoly%cosaoi(isi),csite,nzg,nzs                 &
                          ,cpoly%ntext_soil(:,isi),cpoly%ncol_soil(isi),maxcohort,tuco     &
                          ,rshort_tot,cpoly%met(isi)%rshort_diffuse,cpoly%met(isi)%rlong   &
                          ,daytime,twilight)
            !------------------------------------------------------------------------------!

            !----- Normalize the absorbed radiations. -------------------------------------!
            call scale_ed_radiation(tuco,rshort_tot,cpoly%met(isi)%rshort_diffuse          &
                                   ,cpoly%met(isi)%rlong,cpoly%nighttime(isi),csite)
            !------------------------------------------------------------------------------!

         end do siteloop
      end do polyloop

      !----- Update the average radiation for phenology. ----------------------------------!
      call update_rad_avg(cgrid)
      !------------------------------------------------------------------------------------!

   end if

   return
end subroutine radiate_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will drive the distribution of radiation among crowns, snow layers,  !
! and soil.                                                                                !
!------------------------------------------------------------------------------------------!
subroutine sfcrad_ed(cosz,cosaoi,csite,mzg,mzs,ntext_soil,ncol_soil,maxcohort,tuco         &
                    ,rshort_tot,rshort_diffuse,rlong,daytime,twilight)

   use ed_state_vars        , only : sitetype             & ! structure
                                   , patchtype            ! ! structure
   use canopy_layer_coms    , only : crown_mod            & ! intent(in)
                                   , tai_lyr_max          ! ! intent(in)
   use canopy_radiation_coms, only : icanrad              & ! intent(in)
                                   , cosz_min             & ! intent(in)
                                   , clumping_factor      & ! intent(in)
                                   , par_beam_norm        & ! intent(in)
                                   , par_diff_norm        & ! intent(in)
                                   , nir_beam_norm        & ! intent(in)
                                   , nir_diff_norm        & ! intent(in)
                                   , leaf_trans_vis       & ! intent(in)
                                   , leaf_reflect_vis     & ! intent(in)
                                   , leaf_scatter_vis     & ! intent(out)
                                   , leaf_backscatter_vis & ! intent(out)
                                   , wood_trans_vis       & ! intent(in)
                                   , wood_reflect_vis     & ! intent(in)
                                   , wood_scatter_vis     & ! intent(out)
                                   , wood_backscatter_vis & ! intent(out)
                                   , leaf_trans_nir       & ! intent(in)
                                   , leaf_reflect_nir     & ! intent(in)
                                   , leaf_scatter_nir     & ! intent(out)
                                   , leaf_backscatter_nir & ! intent(out)
                                   , wood_trans_nir       & ! intent(in)
                                   , wood_reflect_nir     & ! intent(in)
                                   , wood_scatter_nir     & ! intent(out)
                                   , wood_backscatter_nir & ! intent(out)
                                   , leaf_emiss_tir       & ! intent(in)
                                   , wood_emiss_tir       & ! intent(in)
                                   , snow_albedo_vis      & ! intent(in)
                                   , snow_albedo_nir      & ! intent(in)
                                   , snow_emiss_tir       & ! intent(in)
                                   , orient_factor        & ! intent(out)
                                   , phi1                 & ! intent(out)
                                   , phi2                 & ! intent(out)
                                   , mu_bar               & ! intent(out)
                                   , radscr
   use soil_coms            , only : soil                 & ! intent(in)
                                   , soilcol              ! ! intent(in)
   use pft_coms             , only : SLA                  ! ! intent(in)
   use consts_coms          , only : stefan               & ! intent(in)
                                   , lnexp_max            ! ! intent(in)
   use ed_max_dims          , only : n_pft                & ! intent(in)
                                   , n_radprof            ! ! intent(in)
   use allometry            , only : h2crownbh            & ! intent(in)
                                   , size2bl              & ! function
                                   , area_indices         ! ! function 
   use ed_misc_coms         , only : ibigleaf             & ! intent(in)
                                   , radfrq               & ! intent(in)
                                   , current_time         ! ! intent(in)
   !$ use omp_lib

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)                  , target      :: csite
   integer                         , intent(in)  :: mzg
   integer                         , intent(in)  :: mzs
   integer         , dimension(mzg), intent(in)  :: ntext_soil
   integer                         , intent(in)  :: ncol_soil
   real                            , intent(in)  :: rshort_tot
   real                            , intent(in)  :: rshort_diffuse
   real                            , intent(in)  :: rlong
   real                            , intent(in)  :: cosaoi
   real                            , intent(in)  :: cosz
   integer                         , intent(in)  :: maxcohort
   logical                         , intent(in)  :: daytime
   logical                         , intent(in)  :: twilight
   integer                         , intent(out) :: tuco
   !----- Local variables. ----------------------------------------------------------------!
   type(patchtype) , pointer                     :: cpatch
   integer                                       :: il
   integer                                       :: ipa
   integer                                       :: ico
   integer                                       :: ipft
   integer                                       :: cohort_count
   integer                                       :: max_cohort_count
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
   real                                          :: ground_par_check
   real                                          :: ground_nir_check
   integer                                       :: ibuff

   integer                                       :: wlr, lvis, lnir
   integer                                       :: lpft, ii
   integer ,allocatable ,dimension(:)            :: npft
   real    ,allocatable ,dimension(:,:)          :: ltv,lrv
   real    ,allocatable, dimension(:,:)          :: ltn,lrn
   real    ,allocatable, dimension(:)            :: alp, aln
   real    ,allocatable ,dimension(:)            :: wrv, srv
   real    ,allocatable ,dimension(:)            :: wrn, srn
   !----- External function. --------------------------------------------------------------!
   real            , external                    :: sngloff
   !----- Local constants. ----------------------------------------------------------------!
   real(kind=8)    , parameter                   :: tiny_offset = 1.d-20
   !---------------------------------------------------------------------------------------!

   open(20,file='lengths.dat')
   read(20,*) lvis, lnir
   read(20,*) lpft
   allocate(npft(lpft))
   
   do ii = 1, lpft
      read(20,*) npft(ii)
   end do
   close(20)

   allocate(ltv(lvis,lpft))
   allocate(lrv(lvis,lpft))
   allocate(alp(lvis))
   allocate(wrv(lvis))
   allocate(srv(lvis))

   allocate(ltn(lnir,lpft))
   allocate(lrn(lnir,lpft))
   allocate(aln(lnir))
   allocate(wrn(lnir))
   allocate(srn(lnir))

   ltv = 0.
   ltn = 0.
   lrv = 0.
   lrn = 0.
   wrv = 0.
   wrn = 0.
   srv = 0.
   srn = 0.

   alp = 0.
   aln = 0.

   open(21,file='trans_par.dat')
   do ii = 1, lpft
      read(21,*), ltv(:,ii)
   end do
   close(21)

   open(22,file='reflect_par.dat')
   do ii = 1, lpft
      read(22,*), lrv(:,ii)
   end do
   close(22)

   open(23,file='trans_nir.dat')
   do ii = 1, lpft
      read(23,*), ltn(:,ii)
   end do
   close(23)

   open(24,file='reflect_nir.dat')
   do ii = 1, lpft
      read(24,*), lrn(:,ii)
   end do
   close(24)

   open(25,file='wood_reflect_par.dat')
   read(25,*), wrv
   close(25)

!   open(26,file='wood_reflect_nir.dat')
!   read(26,*), wrn
!   close(26)

   open(27,file='soil_reflect_par.dat')
   read(27,*), srv
   close(27)

!   open(28,file='soil_reflect_nir.dat')
!   read(28,*), srn
!   close(28)


   !---------------------------------------------------------------------------------------!
   !     Light extinction coefficients.   These are found following CLM technical manual,  !
   ! and the values fall back to ED-2.0 defaults when orient_factor is zero.               !
   !---------------------------------------------------------------------------------------!
   phi1 = 5.d-1 - orient_factor * ( 6.33d-1 + 3.3d-1 * orient_factor )
   phi2 = 8.77d-1 * (1.d0 - 2.d0 * phi1)

   !---------------------------------------------------------------------------------------!
   !     Find the average inverse diffuse optical depth per unit leaf and stem area.       !
   ! We follow CLM technical manual, equation 3.4 only when the orientation factor is      !
   ! non-zero.   Otherwise, we make it 1.d0, which is the limit of that equation when      !
   ! phi2 approaches zero.                                                                 !
   !---------------------------------------------------------------------------------------!
   do ipft = 1, n_pft
      if (orient_factor(ipft) == 0.d0) then
         mu_bar(ipft) = 1.d0
      else
         mu_bar(ipft) = ( 1.d0                                                             &
                        - phi1(ipft) * log(1.d0 + phi2(ipft) / phi1(ipft)) / phi2(ipft) )  &
                        / phi2(ipft)
      end if
   end do
   !---------------------------------------------------------------------------------------!


   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(         &
   !$OMP ibuff,cpatch,cohort_count,tuco,tuco_leaf,    &
   !$OMP ico,bl_lai_each,bl_wai_each,nsoil,colour,    &
   !$OMP albedo_soil_par,albedo_soil_nir,             &
   !$OMP albedo_damp_par,albedo_damp_nir,fcpct,       &
   !$OMP rad_sfcw_par,rad_sfcw_nir,albedo_ground_par, &
   !$OMP albedo_ground_nir,abs_sfcw_par,abs_sfcw_nir, &
   !$OMP emissivity,T_surface,ksn,albedo_sfcw_par,    &
   !$OMP albedo_sfcw_nir,k,sfcw_odepth,fractrans_par, &
   !$OMP fractrans_nir,downward_lw_below,             &
   !$OMP upward_lw_below,upward_lw_above,             &
   !$OMP surface_netabs_longwave,                     &
   !$OMP downward_par_below_beam,                     &
   !$OMP downward_par_below_diffuse,                  &
   !$OMP upward_par_above_diffuse,                    &
   !$OMP downward_nir_below_beam,                     &
   !$OMP downward_nir_below_diffuse,                  &
   !$OMP upward_nir_above_diffuse,                    &
   !$OMP ipft,wleaf_vis,wleaf_nir,wleaf_tir,          &
   !$OMP wwood_vis,wwood_nir,wwood_tir,               &
   !$OMP upward_rshort_above_diffuse,                 &
   !$OMP downward_rshort_below_beam,                  &
   !$OMP downward_rshort_below_diffuse,               &
   !$OMP il,nir_v_beam,nir_v_diffuse,                 &
   !$OMP abs_ground_par,abs_ground_nir )                

   
   !----- Loop over the patches -----------------------------------------------------------!
   do ipa = 1,csite%npatches
      cpatch => csite%patch(ipa)
      
      ibuff = 1
      !$ ibuff = OMP_get_thread_num()+1

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
      tuco         = 0
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
               !----- This will eventually have the index of the tallest used cohort. -----!
               tuco = ico
               !---------------------------------------------------------------------------!

               cohort_count                          = cohort_count + 1
               radscr(ibuff)%pft_array(cohort_count) = cpatch%pft(ico)

               cpatch%bleaf(ico) = size2bl(cpatch%dbh(ico),cpatch%hite(ico), cpatch%pft(ico))

               cpatch%sla(ico) = SLA(cpatch%pft(ico))

               call area_indices(cpatch, ico)

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
               tuco = 1 !---- Dummy variable. ---------------------------------------------!

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
         !      impact on the general circulation model simulated surface climate.  J.     !
         !      Geophys. Res.-Atmosph., 107(D14), 4221, 10.1029/2001JD000809.              !
         !                                                                                 !
         !  Oleson, K.W., et al., 2010: Technical description of version 4.0 of the        !
         !      Community Land Model (CLM). NCAR Technical Note NCAR/TN-478+STR.           !
         !                                                                                 !
         !---------------------------------------------------------------------------------!
         albedo_sfcw_par = snow_albedo_vis + csite%sfcwater_fracliq(ksn,ipa)               &
                                           * ( albedo_damp_par - snow_albedo_vis )
         albedo_sfcw_nir = snow_albedo_nir + csite%sfcwater_fracliq(ksn,ipa)               &
                                           * ( albedo_damp_nir - snow_albedo_nir )
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
            call old_lw_two_stream(emissivity,T_surface,rlong,cohort_count,                &
                                   radscr(ibuff)%pft_array(1:cohort_count),                &
                                   radscr(ibuff)%LAI_array(1:cohort_count),                &
                                   radscr(ibuff)%WAI_array(1:cohort_count),                & 
                                   radscr(ibuff)%CA_array(1:cohort_count),                 &
                                   radscr(ibuff)%leaf_temp_array(1:cohort_count),          &
                                   radscr(ibuff)%wood_temp_array(1:cohort_count),          &
                                   radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count),&
                                   radscr(ibuff)%lw_v_array(1:cohort_count),               &
                                   downward_lw_below,                                      &
                                   upward_lw_below,                                        &
                                   upward_lw_above)
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
            call lw_multiple_scatter(emissivity,T_surface,rlong,cohort_count,              &
                                   radscr(ibuff)%pft_array(1:cohort_count),                &
                                   radscr(ibuff)%LAI_array(1:cohort_count),                &
                                   radscr(ibuff)%WAI_array(1:cohort_count),                & 
                                   radscr(ibuff)%CA_array(1:cohort_count),                 &
                                   radscr(ibuff)%leaf_temp_array(1:cohort_count),          &
                                   radscr(ibuff)%wood_temp_array(1:cohort_count),          &
                                   radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count),&
                                   radscr(ibuff)%lw_v_array(1:cohort_count),               &
                                   downward_lw_below,                                      &
                                   upward_lw_below,                                        &
                                   upward_lw_above)
            !------------------------------------------------------------------------------!





         case (2) 
            !------------------------------------------------------------------------------!
            !    Updated two-stream model.                                                 !
            !------------------------------------------------------------------------------!
            call lw_two_stream(emissivity,T_surface,rlong,cohort_count,                    &
                                   radscr(ibuff)%pft_array(1:cohort_count),                &
                                   radscr(ibuff)%LAI_array(1:cohort_count),                &
                                   radscr(ibuff)%WAI_array(1:cohort_count),                & 
                                   radscr(ibuff)%CA_array(1:cohort_count),                 &
                                   radscr(ibuff)%leaf_temp_array(1:cohort_count),          &
                                   radscr(ibuff)%wood_temp_array(1:cohort_count),          &
                                   radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count),&
                                   radscr(ibuff)%lw_v_array(1:cohort_count),               &
                                   downward_lw_below,upward_lw_below,upward_lw_above)
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

            upward_rshort_above_diffuse = 0.
            downward_rshort_below_beam = 0.
            downward_rshort_below_diffuse = 0.

            do wlr = 1, size(lrv(:,1))

               do ii = 1, lpft
                  leaf_reflect_vis(npft(ii)) = lrv(wlr,ii)
                  leaf_trans_vis(npft(ii)) = ltv(wlr,ii)
               end do

               wood_reflect_vis = wrv(wlr)
               wood_trans_vis = 0.               

               albedo_ground_par = srv(wlr)


               !---------------------------------------------------------------------------------------!
               !     Scattering coefficients.  Contrary to ED-2.1, these values are based on the       !
               ! description by by Sellers (1985) and the CLM technical manual, which includes the     !
               ! leaf orientation factor in the backscattering.  This DOES NOT reduce to ED-2.1 case   !
               ! when the leaf orientation is random.                                                  !
               !---------------------------------------------------------------------------------------!
               
               !---------------------------------------------------------------------------------------!
               !     Forward scattering.                                                               !
               !---------------------------------------------------------------------------------------!
               !----- Visible (PAR). ------------------------------------------------------------------!
               leaf_scatter_vis = leaf_reflect_vis + leaf_trans_vis
               wood_scatter_vis = wood_reflect_vis + wood_trans_vis
               !----- Near infrared (NIR). ------------------------------------------------------------!
               leaf_scatter_nir = leaf_reflect_nir + leaf_trans_nir
               wood_scatter_nir = wood_reflect_nir + wood_trans_nir
               !---------------------------------------------------------------------------------------!
               
               !---------------------------------------------------------------------------------------!
               !      Back-scattering coefficients following CLM.                                      !
               !---------------------------------------------------------------------------------------!
               !----- Visible (PAR). ------------------------------------------------------------------!
               leaf_backscatter_vis = ( leaf_scatter_vis                                               &
                    + 2.5d-1 * ( leaf_reflect_vis - leaf_trans_vis   )               &
                    * ( 1.d0 + orient_factor) ** 2 )                                 &
                    / ( 2.d0 * leaf_scatter_vis )
               wood_backscatter_vis = ( wood_scatter_vis                                               &
                    + 2.5d-1                                                         &
                    * ( wood_reflect_vis - wood_trans_vis   )                        &
                    * ( 1.d0 + orient_factor) ** 2 )                                 &
                    / ( 2.d0 * wood_scatter_vis )
               !----- Near infrared (NIR). ------------------------------------------------------------!
               leaf_backscatter_nir = ( leaf_scatter_nir                                               &
                    + 2.5d-1                                                         &
                    * ( leaf_reflect_nir - leaf_trans_nir   )                        &
                    * ( 1.d0 + orient_factor) ** 2 )                                 &
                    / ( 2.d0 * leaf_scatter_nir )
               wood_backscatter_nir = ( wood_scatter_nir                                               &
                    + 2.5d-1                                                         &
                    * ( wood_reflect_nir - wood_trans_nir   )                        &
                    * ( 1.d0 + orient_factor) ** 2 )                                 &
                    / ( 2.d0 * wood_scatter_nir )
               
               
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
                       ,cohort_count                                         &
                       ,radscr(ibuff)%pft_array(1:cohort_count)              &
                       ,radscr(ibuff)%LAI_array(1:cohort_count)                 &
                       ,radscr(ibuff)%WAI_array(1:cohort_count)                 &
                       ,radscr(ibuff)%CA_array(1:cohort_count)                  &
                       ,radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count) &
                       ,radscr(ibuff)%par_v_beam_array(1:cohort_count)          &
                       ,radscr(ibuff)%par_v_diffuse_array(1:cohort_count)       &
                       ,radscr(ibuff)%rshort_v_beam_array(1:cohort_count)       &
                       ,radscr(ibuff)%rshort_v_diffuse_array(1:cohort_count)    &
                       ,downward_par_below_beam                              &
                       ,downward_par_below_diffuse                           &
                       ,upward_par_above_diffuse                             &
                       ,downward_nir_below_beam                              &
                       ,downward_nir_below_diffuse                           &
                       ,upward_nir_above_diffuse                             &
                       ,radscr(ibuff)%par_level_beam(1:cohort_count)         &
                       ,radscr(ibuff)%par_level_diffd(1:cohort_count)        &
                       ,radscr(ibuff)%par_level_diffu(1:cohort_count)        &
                       ,radscr(ibuff)%light_level_array(1:cohort_count)      &
                       ,radscr(ibuff)%light_beam_level_array(1:cohort_count) &
                       ,radscr(ibuff)%light_diff_level_array(1:cohort_count))
                  
                  !---------------------------------------------------------------------------!
                  
               case (1)
                  !---------------------------------------------------------------------------!
                  !      Multiple-scatter model.                                              !
                  !---------------------------------------------------------------------------!
                  
                  call sw_multiple_scatter(albedo_ground_par,albedo_ground_nir,cosaoi         &
                       ,cohort_count                                         &
                       ,radscr(ibuff)%pft_array(1:cohort_count)              &
                       ,radscr(ibuff)%LAI_array(1:cohort_count)                 &
                       ,radscr(ibuff)%WAI_array(1:cohort_count)                 &
                       ,radscr(ibuff)%CA_array(1:cohort_count)                  &
                       ,radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count) &
                       ,radscr(ibuff)%par_v_beam_array(1:cohort_count)          &
                       ,radscr(ibuff)%par_v_diffuse_array(1:cohort_count)       &
                       ,radscr(ibuff)%rshort_v_beam_array(1:cohort_count)       &
                       ,radscr(ibuff)%rshort_v_diffuse_array(1:cohort_count)    &
                       ,downward_par_below_beam                              &
                       ,downward_par_below_diffuse                           &
                       ,upward_par_above_diffuse                             &
                       ,downward_nir_below_beam                              &
                       ,downward_nir_below_diffuse                           &
                       ,upward_nir_above_diffuse                             &
                       ,radscr(ibuff)%par_level_beam(1:cohort_count)         &
                       ,radscr(ibuff)%par_level_diffd(1:cohort_count)        &
                       ,radscr(ibuff)%par_level_diffu(1:cohort_count)        &
                       ,radscr(ibuff)%light_level_array(1:cohort_count)      &
                       ,radscr(ibuff)%light_beam_level_array(1:cohort_count) &
                       ,radscr(ibuff)%light_diff_level_array(1:cohort_count))

                  !---------------------------------------------------------------------------!
               case (2)
                  !---------------------------------------------------------------------------!
                  !    Updated two-stream model.                                              !
                  !---------------------------------------------------------------------------!

                  call sw_two_stream(albedo_ground_par,albedo_ground_nir,cosaoi                &
                       ,cohort_count                                         &
                       ,radscr(ibuff)%pft_array(1:cohort_count)              &
                       ,radscr(ibuff)%LAI_array(1:cohort_count)                 &
                       ,radscr(ibuff)%WAI_array(1:cohort_count)                 &
                       ,radscr(ibuff)%CA_array(1:cohort_count)                  &
                       ,radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count) &
                       ,radscr(ibuff)%par_v_beam_array(1:cohort_count)          &
                       ,radscr(ibuff)%par_v_diffuse_array(1:cohort_count)       &
                       ,radscr(ibuff)%rshort_v_beam_array(1:cohort_count)       &
                       ,radscr(ibuff)%rshort_v_diffuse_array(1:cohort_count)    &
                       ,downward_par_below_beam                              &
                       ,downward_par_below_diffuse                           &
                       ,upward_par_above_diffuse                             &
                       ,downward_nir_below_beam                              &
                       ,downward_nir_below_diffuse                           &
                       ,upward_nir_above_diffuse                             &
                       ,radscr(ibuff)%par_level_beam(1:cohort_count)         &
                       ,radscr(ibuff)%par_level_diffd(1:cohort_count)        &
                       ,radscr(ibuff)%par_level_diffu(1:cohort_count)        &
                       ,radscr(ibuff)%light_level_array(1:cohort_count)      &
                       ,radscr(ibuff)%light_beam_level_array(1:cohort_count) &
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
               ! (par_beam_norm,par_diff_norm,nir_beam_norm and nir_diff_norm are site level  !
               csite%albedo_par (ipa) = upward_par_above_diffuse                              &
                    / sngloff( par_beam_norm + par_diff_norm, tiny_offset )
               csite%albedo_nir (ipa) = upward_nir_above_diffuse                              &
                    / sngloff( nir_beam_norm + nir_diff_norm, tiny_offset )
               csite%albedo     (ipa) = upward_rshort_above_diffuse
               !------------------------------------------------------------------------------!
               alp(wlr) = csite%albedo_par(ipa)
            end do

            do wlr = 1, size(lrn(:,1))

               do ii = 1, lpft
                  leaf_reflect_nir(npft(ii)) = lrn(wlr,ii)
                  leaf_trans_nir(npft(ii)) = ltn(wlr,ii)
               end do

               !---------------------------------------------------------------------------------------!
               !     Scattering coefficients.  Contrary to ED-2.1, these values are based on the       !
               ! description by by Sellers (1985) and the CLM technical manual, which includes the     !
               ! leaf orientation factor in the backscattering.  This DOES NOT reduce to ED-2.1 case   !
               ! when the leaf orientation is random.                                                  !
               !---------------------------------------------------------------------------------------!
               
               !---------------------------------------------------------------------------------------!
               !     Forward scattering.                                                               !
               !---------------------------------------------------------------------------------------!
               !----- Visible (PAR). ------------------------------------------------------------------!
               leaf_scatter_vis = leaf_reflect_vis + leaf_trans_vis
               wood_scatter_vis = wood_reflect_vis + wood_trans_vis
               !----- Near infrared (NIR). ------------------------------------------------------------!
               leaf_scatter_nir = leaf_reflect_nir + leaf_trans_nir
               wood_scatter_nir = wood_reflect_nir + wood_trans_nir
               !---------------------------------------------------------------------------------------!
               
               !---------------------------------------------------------------------------------------!
               !      Back-scattering coefficients following CLM.                                      !
               !---------------------------------------------------------------------------------------!
               !----- Visible (PAR). ------------------------------------------------------------------!
               leaf_backscatter_vis = ( leaf_scatter_vis                                               &
                    + 2.5d-1 * ( leaf_reflect_vis - leaf_trans_vis   )               &
                    * ( 1.d0 + orient_factor) ** 2 )                                 &
                    / ( 2.d0 * leaf_scatter_vis )
               wood_backscatter_vis = ( wood_scatter_vis                                               &
                    + 2.5d-1                                                         &
                    * ( wood_reflect_vis - wood_trans_vis   )                        &
                    * ( 1.d0 + orient_factor) ** 2 )                                 &
                    / ( 2.d0 * wood_scatter_vis )
               !----- Near infrared (NIR). ------------------------------------------------------------!
               leaf_backscatter_nir = ( leaf_scatter_nir                                               &
                    + 2.5d-1                                                         &
                    * ( leaf_reflect_nir - leaf_trans_nir   )                        &
                    * ( 1.d0 + orient_factor) ** 2 )                                 &
                    / ( 2.d0 * leaf_scatter_nir )
               wood_backscatter_nir = ( wood_scatter_nir                                               &
                    + 2.5d-1                                                         &
                    * ( wood_reflect_nir - wood_trans_nir   )                        &
                    * ( 1.d0 + orient_factor) ** 2 )                                 &
                    / ( 2.d0 * wood_scatter_nir )
               
               
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
                       ,cohort_count                                         &
                       ,radscr(ibuff)%pft_array(1:cohort_count)              &
                       ,radscr(ibuff)%LAI_array(1:cohort_count)                 &
                       ,radscr(ibuff)%WAI_array(1:cohort_count)                 &
                       ,radscr(ibuff)%CA_array(1:cohort_count)                  &
                       ,radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count) &
                       ,radscr(ibuff)%par_v_beam_array(1:cohort_count)          &
                       ,radscr(ibuff)%par_v_diffuse_array(1:cohort_count)       &
                       ,radscr(ibuff)%rshort_v_beam_array(1:cohort_count)       &
                       ,radscr(ibuff)%rshort_v_diffuse_array(1:cohort_count)    &
                       ,downward_par_below_beam                              &
                       ,downward_par_below_diffuse                           &
                       ,upward_par_above_diffuse                             &
                       ,downward_nir_below_beam                              &
                       ,downward_nir_below_diffuse                           &
                       ,upward_nir_above_diffuse                             &
                       ,radscr(ibuff)%par_level_beam(1:cohort_count)         &
                       ,radscr(ibuff)%par_level_diffd(1:cohort_count)        &
                       ,radscr(ibuff)%par_level_diffu(1:cohort_count)        &
                       ,radscr(ibuff)%light_level_array(1:cohort_count)      &
                       ,radscr(ibuff)%light_beam_level_array(1:cohort_count) &
                       ,radscr(ibuff)%light_diff_level_array(1:cohort_count))
                  
                  !---------------------------------------------------------------------------!
                  
               case (1)
                  !---------------------------------------------------------------------------!
                  !      Multiple-scatter model.                                              !
                  !---------------------------------------------------------------------------!
                  
                  call sw_multiple_scatter(albedo_ground_par,albedo_ground_nir,cosaoi         &
                       ,cohort_count                                         &
                       ,radscr(ibuff)%pft_array(1:cohort_count)              &
                       ,radscr(ibuff)%LAI_array(1:cohort_count)                 &
                       ,radscr(ibuff)%WAI_array(1:cohort_count)                 &
                       ,radscr(ibuff)%CA_array(1:cohort_count)                  &
                       ,radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count) &
                       ,radscr(ibuff)%par_v_beam_array(1:cohort_count)          &
                       ,radscr(ibuff)%par_v_diffuse_array(1:cohort_count)       &
                       ,radscr(ibuff)%rshort_v_beam_array(1:cohort_count)       &
                       ,radscr(ibuff)%rshort_v_diffuse_array(1:cohort_count)    &
                       ,downward_par_below_beam                              &
                       ,downward_par_below_diffuse                           &
                       ,upward_par_above_diffuse                             &
                       ,downward_nir_below_beam                              &
                       ,downward_nir_below_diffuse                           &
                       ,upward_nir_above_diffuse                             &
                       ,radscr(ibuff)%par_level_beam(1:cohort_count)         &
                       ,radscr(ibuff)%par_level_diffd(1:cohort_count)        &
                       ,radscr(ibuff)%par_level_diffu(1:cohort_count)        &
                       ,radscr(ibuff)%light_level_array(1:cohort_count)      &
                       ,radscr(ibuff)%light_beam_level_array(1:cohort_count) &
                       ,radscr(ibuff)%light_diff_level_array(1:cohort_count))

                  !---------------------------------------------------------------------------!
               case (2)
                  !---------------------------------------------------------------------------!
                  !    Updated two-stream model.                                              !
                  !---------------------------------------------------------------------------!

                  call sw_two_stream(albedo_ground_par,albedo_ground_nir,cosaoi                &
                       ,cohort_count                                         &
                       ,radscr(ibuff)%pft_array(1:cohort_count)              &
                       ,radscr(ibuff)%LAI_array(1:cohort_count)                 &
                       ,radscr(ibuff)%WAI_array(1:cohort_count)                 &
                       ,radscr(ibuff)%CA_array(1:cohort_count)                  &
                       ,radscr(ibuff)%radprof_array(1:n_radprof,1:cohort_count) &
                       ,radscr(ibuff)%par_v_beam_array(1:cohort_count)          &
                       ,radscr(ibuff)%par_v_diffuse_array(1:cohort_count)       &
                       ,radscr(ibuff)%rshort_v_beam_array(1:cohort_count)       &
                       ,radscr(ibuff)%rshort_v_diffuse_array(1:cohort_count)    &
                       ,downward_par_below_beam                              &
                       ,downward_par_below_diffuse                           &
                       ,upward_par_above_diffuse                             &
                       ,downward_nir_below_beam                              &
                       ,downward_nir_below_diffuse                           &
                       ,upward_nir_above_diffuse                             &
                       ,radscr(ibuff)%par_level_beam(1:cohort_count)         &
                       ,radscr(ibuff)%par_level_diffd(1:cohort_count)        &
                       ,radscr(ibuff)%par_level_diffu(1:cohort_count)        &
                       ,radscr(ibuff)%light_level_array(1:cohort_count)      &
                       ,radscr(ibuff)%light_beam_level_array(1:cohort_count) &
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
               ! (par_beam_norm,par_diff_norm,nir_beam_norm and nir_diff_norm are site level  !
               csite%albedo_par (ipa) = upward_par_above_diffuse                              &
                    / sngloff( par_beam_norm + par_diff_norm, tiny_offset )
               csite%albedo_nir (ipa) = upward_nir_above_diffuse                              &
                    / sngloff( nir_beam_norm + nir_diff_norm, tiny_offset )
               csite%albedo     (ipa) = upward_rshort_above_diffuse
               !------------------------------------------------------------------------------!

               aln(wlr) = csite%albedo_nir(ipa)
            end do
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

         open(31,file='albedo_par.dat')
         write(31,*) alp
         close(31)

         open(32,file='albedo_nir.dat')
         write(32,*) aln
         close(32)

         open(33,file='LAI.dat')
         write(33,*) cpatch%lai
         close(33)
         
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
                  wleaf_tir = sngloff( ( leaf_emiss_tir(ipft) * radscr(ibuff)%LAI_array(il) ) &
                                     / ( leaf_emiss_tir(ipft) * radscr(ibuff)%LAI_array(il)   &
                                       + wood_emiss_tir(ipft) * radscr(ibuff)%WAI_array(il) ) &
                                     , tiny_offset )
                  wwood_vis = 1. - wleaf_vis
                  wwood_nir = 1. - wleaf_nir
                  wwood_tir = 1. - wleaf_tir
                  !------------------------------------------------------------------------!

                  !----- Find the near infrared absorption, so we average things properly. !
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
                  cpatch%par_l_beam       (ico) = radscr(ibuff)%par_v_beam_array(il) *     &
                                                  wleaf_vis
                  cpatch%par_l_diffuse    (ico) = radscr(ibuff)%par_v_diffuse_array(il) *  &
                                                  wleaf_vis
                  !------ Total short wave radiation (PAR+NIR). ---------------------------!
                  cpatch%rshort_l_beam    (ico) = radscr(ibuff)%par_v_beam_array    (il) * &
                                                  wleaf_vis + nir_v_beam * wleaf_nir
                  cpatch%rshort_l_diffuse (ico) = radscr(ibuff)%par_v_diffuse_array (il) * &
                                                  wleaf_vis + nir_v_diffuse * wleaf_nir
                  cpatch%rshort_w_beam    (ico) = radscr(ibuff)%par_v_beam_array    (il) * &
                                                  wwood_vis + nir_v_beam * wwood_nir
                  cpatch%rshort_w_diffuse (ico) = radscr(ibuff)%par_v_diffuse_array (il) * &
                                                  wwood_vis + nir_v_diffuse * wwood_nir
                  !----- Thermal infra-red (long wave). -----------------------------------!
                  cpatch%rlong_l          (ico) = radscr(ibuff)%lw_v_array(il) * wleaf_tir
                  cpatch%rlong_w          (ico) = radscr(ibuff)%lw_v_array(il) * wwood_tir
                  !------------------------------------------------------------------------!

                  !----- Save the light levels. -------------------------------------------!
                  cpatch%light_level(ico)       = sngloff(radscr(ibuff)%light_level_array(il)            &
                                                         ,tiny_offset )
                  cpatch%light_level_beam(ico)  = sngloff(radscr(ibuff)%light_beam_level_array(il)       &
                                                         ,tiny_offset )
                  cpatch%light_level_diff(ico)  = sngloff(radscr(ibuff)%light_diff_level_array(il)       &
                                                         ,tiny_offset )

                  cpatch%par_level_beam(ico)  = sngloff(radscr(ibuff)%par_level_beam(il),tiny_offset)
                  cpatch%par_level_diffu(ico) = sngloff(radscr(ibuff)%par_level_diffu(il),tiny_offset)
                  cpatch%par_level_diffd(ico) = sngloff(radscr(ibuff)%par_level_diffd(il),tiny_offset)

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
                  ! average between the area and absorptance.  We must treat the           !
                  ! visible and near infrared separately.                                  !
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
                  wleaf_tir = sngloff( ( leaf_emiss_tir(ipft) * radscr(ibuff)%LAI_array(il) ) &
                                     / ( leaf_emiss_tir(ipft) * radscr(ibuff)%LAI_array(il)   &
                                       + wood_emiss_tir(ipft) * radscr(ibuff)%WAI_array(il) ) &
                                     , tiny_offset    )
                  wwood_vis = 1. - wleaf_vis
                  wwood_nir = 1. - wleaf_nir
                  wwood_tir = 1. - wleaf_tir
                  !------------------------------------------------------------------------!
            


                  !------------------------------------------------------------------------!
                  !     Find the near infrared absorption, so we average things            !
                  ! properly.                                                              !
                  !------------------------------------------------------------------------!
                  nir_v_beam    = radscr(ibuff)%rshort_v_beam_array(il) -                  &
                                  radscr(ibuff)%par_v_beam_array(il)
                  nir_v_diffuse = radscr(ibuff)%rshort_v_diffuse_array(il) -               &
                                  radscr(ibuff)%par_v_diffuse_array(il)
                  !------------------------------------------------------------------------!

                  !------------------------------------------------------------------------!
                  !    Split the layer radiation between leaf and branchwood.              !
                  !------------------------------------------------------------------------!
                  !------ Visible (PAR), only leaves need this (photsynthesis model). -----!
                  cpatch%par_l_beam       (1) = cpatch%par_l_beam(1)                       &
                                              + radscr(ibuff)%par_v_beam_array(il)*wleaf_vis
                  cpatch%par_l_diffuse    (1) = cpatch%par_l_diffuse(1)                    &
                                              + radscr(ibuff)%par_v_diffuse_array(il)*wleaf_vis
                  !------ Total short wave radiation (PAR+NIR). ---------------------------!
                  cpatch%rshort_l_beam    (1) = cpatch%rshort_l_beam(1)                &
                                              + radscr(ibuff)%par_v_beam_array(il)*wleaf_vis &
                                              + nir_v_beam                  * wleaf_nir
                  cpatch%rshort_l_diffuse (1) = cpatch%rshort_l_diffuse (1)                &
                                              + radscr(ibuff)%par_v_diffuse_array    (il) * wleaf_vis    &
                                              + nir_v_diffuse               * wleaf_nir
                  cpatch%rshort_w_beam    (1) = cpatch%rshort_w_beam    (1)                &
                                              + radscr(ibuff)%par_v_beam_array       (il) * wwood_vis    &
                                              + nir_v_beam                  * wwood_nir
                  cpatch%rshort_w_diffuse (1) = cpatch%rshort_w_diffuse (1)                &
                                              + radscr(ibuff)%par_v_diffuse_array    (il) * wwood_vis    &
                                              + nir_v_diffuse               * wwood_nir
                  !----- Thermal infra-red (long wave). -----------------------------------!
                  cpatch%rlong_l          (1) = cpatch%rlong_l (1)                         &
                                              + radscr(ibuff)%lw_v_array    (il) * wleaf_tir
                  cpatch%rlong_w          (1) = cpatch%rlong_w (1)                         &
                                              + radscr(ibuff)%lw_v_array    (il) * wwood_tir
                  !------------------------------------------------------------------------!
               end do
               !---------------------------------------------------------------------------!



               !----- Save the light levels as the median level. --------------------------!
               il = ceiling(real(cohort_count)/2.0)
               cpatch%light_level      (1) = sngloff(radscr(ibuff)%light_level_array     (il)            &
                                                      ,tiny_offset )
               cpatch%light_level_beam (1) = sngloff(radscr(ibuff)%light_beam_level_array(il)            &
                                                      ,tiny_offset )
               cpatch%light_level_diff (1) = sngloff(radscr(ibuff)%light_diff_level_array(il)            &
                                                      ,tiny_offset )

               cpatch%par_level_beam(1)  = sngloff(radscr(ibuff)%par_level_beam(il),tiny_offset)
               cpatch%par_level_diffu(1) = sngloff(radscr(ibuff)%par_level_diffu(il),tiny_offset)
               cpatch%par_level_diffd(1) = sngloff(radscr(ibuff)%par_level_diffd(il),tiny_offset)



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
      end if


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

   end do

   !---------------------------------------------------------------------------------------!
   !$OMP END PARALLEL DO

   return
end subroutine sfcrad_ed
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine scale_ed_radiation(tuco,rshort,rshort_diffuse,rlong,nighttime,csite)

   use ed_state_vars        , only : sitetype             & ! intent(in)
                                   , patchtype            ! ! intent(in)
   use ed_misc_coms         , only : writing_long         & ! intent(in)
                                   , radfrq               & ! intent(in)
                                   , frqsum               ! ! intent(in)
   use canopy_radiation_coms, only : cosz_min             ! ! intent(in)
   !$ use omp_lib
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)  , target     :: csite
   integer         , intent(in) :: tuco
   real            , intent(in) :: rshort
   real            , intent(in) :: rshort_diffuse
   real            , intent(in) :: rlong
   logical         , intent(in) :: nighttime
   !----- Local variables. ----------------------------------------------------------------!
   type(patchtype) , pointer    :: cpatch
   integer                      :: ipa,ico, k
   !----- This should be false unless you really want to turn off radiation. --------------!
   logical         , parameter  :: skip_rad = .false.
   !----- External functions. -------------------------------------------------------------!
   real            , external   :: sngloff
   !----- Locally saved variables. --------------------------------------------------------!
   real              , save    :: radfrq_o_frqsum
   logical           , save    :: first_time = .true.
   !---------------------------------------------------------------------------------------!


   !----- Assign the constant scaling factor. ---------------------------------------------!
   if (first_time) then
      first_time      = .false.
      radfrq_o_frqsum = radfrq / frqsum
   end if
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     This block skips radiation.  Obviously this only makes sense for very theoretical !
   ! tests, because plants kind of like light...                                           !
   !---------------------------------------------------------------------------------------!
   if (skip_rad) then
      skip_patchloop: do ipa = 1, csite%npatches
         cpatch => csite%patch(ipa)
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

         csite%rlong_s(ipa)       = 0.
         csite%rlong_g(ipa)       = 0.
      end do skip_patchloop
      return
   end if
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     For normal runs, we add the scale to radiation, and also integrate averages.      !
   !---------------------------------------------------------------------------------------!
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(cpatch,ico,k)

   patchloop: do ipa = 1,csite%npatches
      cpatch => csite%patch(ipa)

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
                                              + cpatch%par_level_beam(ico) * radfrq_o_frqsum

         cpatch%fmean_par_level_diffd   (ico) = cpatch%fmean_par_level_diffd   (ico)       &
                                              + cpatch%par_level_diffd(ico) * radfrq_o_frqsum

         cpatch%fmean_par_level_diffu   (ico) = cpatch%fmean_par_level_diffu   (ico)       &
                                              + cpatch%par_level_diffu(ico) * radfrq_o_frqsum

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
   end do patchloop
   !---------------------------------------------------------------------------------------!
   !$OMP END PARALLEL DO


   return
end subroutine scale_ed_radiation
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine calculates angle of incidence based on local slope and aspect.       !
!------------------------------------------------------------------------------------------!
subroutine angle_of_incid(aoi,cosz,solar_hour_aspect,slope,terrain_aspect)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real, intent(in)  :: cosz              ! Cosine of zenithal angle
   real, intent(in)  :: slope             ! Terrain slope
   real, intent(in)  :: solar_hour_aspect ! horizontal location of the sun defined with the
                                          !    same reference as terrain aspect.
   real, intent(in)  :: terrain_aspect    ! Terrain aspect
   real, intent(out) :: aoi               ! Angle of incidence
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)      :: cosz8             ! Double prec. counterpart of cosz
   real(kind=8)      :: sinz8             ! Sine of zenithal angle
   real(kind=8)      :: slope8            ! Double prec. counterpart of slope
   real(kind=8)      :: sh_asp8           ! Double prec. counterpart of solar_hour_aspect
   real(kind=8)      :: terr_asp8         ! Double prec. counterpart of terrain_aspect
   real(kind=8)      :: aoi8              ! Double prec. counterpart of aoi
   !----- Local parameters. ---------------------------------------------------------------!
   real(kind=8), parameter :: tiny_offset=1.d-20
   !----- External functions. -------------------------------------------------------------!
   real        , external  :: sngloff           
   !---------------------------------------------------------------------------------------!

   cosz8     = dble(cosz)
   sinz8     = sqrt(1.d0-cosz8*cosz8)
   slope8    = dble(slope)
   sh_asp8   = dble(solar_hour_aspect)
   terr_asp8 = dble(terrain_aspect)
   if (cosz8 < 0.d0) then
      aoi8 = 0.d0
   else
      aoi8 = max(0.d0, cosz8*dcos(slope8) + sinz8*dsin(slope8)*dcos(sh_asp8-terr_asp8))
   end if

   aoi = sngloff(aoi8,tiny_offset) 

   return
end subroutine angle_of_incid
!==========================================================================================!
!==========================================================================================!

end module radiate_driver_module

!==========================================================================================!
!==========================================================================================!
!     This function computes the cosine of the zenith angle.                               !
!------------------------------------------------------------------------------------------!
real function ed_zen(plon,plat,when)
   use ed_misc_coms , only : simtime      ! ! structure
   use consts_coms  , only : pio1808      & ! intent(in)
                           , twopi8       & ! intent(in)
                           , hr_sec8      & ! intent(in)
                           , tiny_num8    ! ! intent(in)
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   real(kind=4) , intent(in) :: plon
   real(kind=4) , intent(in) :: plat
   type(simtime), intent(in) :: when
   !------ Local variables. ---------------------------------------------------------------!
   integer                   :: doy     ! Day of year ("Julian" day)
   real(kind=8)              :: declin  ! Declination
   real(kind=8)              :: sdec    ! Sine of declination
   real(kind=8)              :: cdec    ! Cosine of declination
   real(kind=8)              :: dayhr   ! Hour of day 
   real(kind=8)              :: radlat  ! Latitude in radians
   real(kind=8)              :: clat    ! Cosine of latitude 
   real(kind=8)              :: slat    ! Sine of latitude
   real(kind=8)              :: dayhrr  ! Hour of day in radians
   real(kind=8)              :: hrangl  ! Hour angle
   !----- Local constants. ----------------------------------------------------------------!
   real(kind=8), parameter   :: capri    = -2.344d1 ! Tropic of Capricornium latitude
   real(kind=8), parameter   :: ndaysnl  =  3.65d2  ! Number of days of year (no leap years)
   real(kind=8), parameter   :: ndayslp  =  3.66d2  ! Number of days of year (leap years)
   integer     , parameter   :: shsummer = -10      ! DoY of Southern Hemisphere summer
   !----- External functions. -------------------------------------------------------------!
   integer     , external    :: julday  ! Function to find day of year ("Julian" day)
   logical     , external    :: isleap  ! Function to determine whether the year is leap
   real        , external    :: sngloff ! Function to safely convert double to single prec.
   !---------------------------------------------------------------------------------------!



   !----- Find the day of the year. -------------------------------------------------------!
   doy    = julday(when%month, when%date, when%year)
   !---------------------------------------------------------------------------------------!



   !----- Find the hour angle, then get cosine of zenith angle. ---------------------------!
   dayhr = dble(when%time) / hr_sec8
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !   declin is the solar latitude in degrees (also known as declination).                !
   !   sdec - sine of declination                                                          !
   !   cdec - cosine of declination.                                                       !
   !---------------------------------------------------------------------------------------!
   if (isleap(when%year)) then
      declin = capri * cos(twopi8 * dble(doy - shsummer) / ndaysnl) * pio1808
   else
      declin = capri * cos(twopi8 * dble(doy - shsummer) / ndayslp) * pio1808
   end if
   sdec   = dsin(declin)
   cdec   = dcos(declin)
   !---------------------------------------------------------------------------------------!

   !----- Find the latitude in radians. ---------------------------------------------------!
   radlat = dble(plat) * pio1808
   clat   = dcos(radlat)
   slat   = dsin(radlat)
   !---------------------------------------------------------------------------------------!

   !------ Find the hour angle. -----------------------------------------------------------!
   dayhrr = dmod(dayhr+dble(plon)/1.5d1+2.4d1,2.4d1)
   hrangl = 1.5d1 * (dayhrr - 1.2d1) * pio1808

   ed_zen = sngloff(slat * sdec + clat * cdec * dcos(hrangl),tiny_num8)

   return
end function ed_zen
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     This function computes the average secant of the daytime zenith angle.  In case the  !
! period of integration is that accounts for the zenith angle is 0. or less than one time  !
! step, then the average is the actual value.  Night-time periods are ignored and if there !
! is no daytime value, then we set it to 0.                                                !
!------------------------------------------------------------------------------------------!
real function mean_daysecz(plon,plat,whena,dt,tmax)
   use update_derived_props_module
   use ed_misc_coms         , only : simtime     ! ! structure
   use canopy_radiation_coms, only : cosz_min    ! ! intent(in)
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   real(kind=4) , intent(in) :: plon
   real(kind=4) , intent(in) :: plat
   type(simtime), intent(in) :: whena
   real(kind=4) , intent(in) :: dt
   real(kind=4) , intent(in) :: tmax
   !------ Local variables. ---------------------------------------------------------------!
   type(simtime)             :: now          ! Current time
   integer                   :: is           ! Step counter
   integer                   :: nsteps       ! Number of steps to perform the average
   real                      :: dtfit        ! Delta-t that nicely fits within tmax
   real                      :: dtnow        ! Delta-t for this time
   real                      :: cosz         ! Declination
   real                      :: daytot       ! Total time that was daytime
   real                      :: mean_daycosz ! Average cosine of zenith angle
   !----- External functions. -------------------------------------------------------------!
   real(kind=4), external    :: ed_zen   ! Function to find day of year ("Julian" day)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Check whether tmax is less than the time step.  In case it is, we only have one   !
   ! time, so we don't need to do the average.                                             !
   !---------------------------------------------------------------------------------------!
   if (dt >= tmax) then
      !----- Less than one time, no average necessary. ------------------------------------!
      cosz = ed_zen(plon,plat,whena)
      if (cosz > cosz_min) then
         mean_daysecz = 1.0 / cosz
      else
         !----- Night-time, set the mean to zero. -----------------------------------------!
         mean_daysecz = 0.0
      end if

   else
      !------------------------------------------------------------------------------------!
      !     Several times, first find the number of steps, then the delta-t that fits      !
      ! nicely within the time span.                                                       !
      !------------------------------------------------------------------------------------!
      nsteps = ceiling(tmax / dt)
      dtfit  = tmax / real(nsteps)
      !------------------------------------------------------------------------------------!

      mean_daycosz = 0.0
      daytot       = 0.0
      do is=1,nsteps
         !----- Get the current time. -----------------------------------------------------!
         now   = whena
         dtnow = dtfit * (real(is) - 0.5)
         call update_model_time_dm(now,dtnow)

         !----- Get the cosine of the zenith angle. ---------------------------------------!
         cosz = ed_zen(plon,plat,now)

         !----- Add to the integral only if it this value is valid. -----------------------!
         if (cosz > cosz_min) then
            mean_daycosz = mean_daycosz + dtfit * cosz 
            daytot       = daytot       + dtfit
         end if
         !---------------------------------------------------------------------------------!
      end do
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Find the normalisation factor.                                                 !
      !------------------------------------------------------------------------------------!
      if (daytot > 0.0 .and. mean_daycosz > 0.0) then
         mean_daycosz = mean_daycosz / daytot
         mean_daysecz = 1.0 / mean_daycosz
      else
         mean_daysecz = 0.0
      end if
      !------------------------------------------------------------------------------------!
   end if

   return
end function mean_daysecz
!==========================================================================================!
!==========================================================================================!

