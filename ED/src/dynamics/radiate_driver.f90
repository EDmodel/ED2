!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the two-stream radiation scheme.  This is called every  !
! step, but not every sub-step.                                                            !
!------------------------------------------------------------------------------------------!
subroutine radiate_driver(cgrid)
   use ed_misc_coms             , only : current_time          & ! intent(in)
                                    , radfrq                & ! intent(in)
                                    , dtlsm                 ! ! intent(in)
   use ed_state_vars         , only : edtype                & ! structure
                                    , polygontype           & ! structure
                                    , sitetype              & ! structure
                                    , patchtype             ! ! structure
   use canopy_radiation_coms , only : visible_fraction_dir  & ! intent(in)
                                    , visible_fraction_dif  & ! intent(in)
                                    , rlong_min             ! ! intent(in)
   use consts_coms           , only : pio180                ! ! intent(in)
   implicit none
   !----- Argument. -----------------------------------------------------------------------!
   type(edtype)     , target  :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer :: cpoly
   type(sitetype)   , pointer :: csite
   type(patchtype)  , pointer :: cpatch
   real                       :: total_beam_radiation
   integer                    :: maxcohort
   integer                    :: ipy,isi,ipa
   real                       :: hrangl
   !---------------------------------------------------------------------------------------!


   !----- Check whether it is time to update radiative fluxes and heating rates -----------!
   if (mod(current_time%time + .001,radfrq) < dtlsm) then

      !----- Compute solar zenith angle [cosz] --------------------------------------------!
      call solar_zenith(cgrid)

      !----- Loop over polygons and sites. ------------------------------------------------!

      polyloop: do ipy = 1,cgrid%npolygons

         cpoly => cgrid%polygon(ipy)

         siteloop: do isi = 1,cpoly%nsites

            csite => cpoly%site(isi)

            !------------------------------------------------------------------------------!
            !     Compute the visible fraction of diffuse and beam radiation needed by the !
            ! radiative transfer routine.                                                  !
            !------------------------------------------------------------------------------!
            total_beam_radiation = cpoly%met(isi)%nir_beam + cpoly%met(isi)%par_beam

            if (total_beam_radiation > 0.0) then
               visible_fraction_dir = cpoly%met(isi)%par_beam / total_beam_radiation
            else
               visible_fraction_dir = 0.5
            end if

            if(cpoly%met(isi)%rshort_diffuse > 0.0)then
               visible_fraction_dif = cpoly%met(isi)%par_diffuse                           &
                                    / cpoly%met(isi)%rshort_diffuse
            else
               visible_fraction_dif = 0.5
            end if

            !----- Update angle of incidence ----------------------------------------------!
            hrangl = 15.0 * pio180                                                         &
                   * (mod(current_time%time + cgrid%lon(ipy) / 15.0 + 24.0, 24.0) - 12.0)
            
            call angle_of_incid(cpoly%cosaoi(isi), cgrid%cosz(ipy), hrangl                 &
                               ,cpoly%slope(isi) * pio180, cpoly%aspect(isi) * pio180)


            !------------------------------------------------------------------------------!
            !    Loop over subgrid-scale patches.  These routines can be done as arrays.   !
            !------------------------------------------------------------------------------!
            maxcohort = 1
            do ipa = 1,csite%npatches
               cpatch=>csite%patch(ipa)
               if ( cpatch%ncohorts>maxcohort ) maxcohort = cpatch%ncohorts
            end do


            !----- Get unnormalized radiative transfer information. -----------------------!
            call sfcrad_ed(cgrid%cosz(ipy),cpoly%cosaoi(isi),csite,maxcohort            &
                             ,cpoly%met(isi)%rshort)

            !----- Normalize the absorbed radiations. -------------------------------------!
            call scale_ed_radiation(cpoly%met(isi)%rshort,cpoly%met(isi)%rshort_diffuse &
                                      ,cpoly%met(isi)%rlong,csite)

            !----- Update all albedos and rlongup. ----------------------------------------!
            call ed2land_radiation(cpoly,isi)
         end do siteloop
      end do polyloop

      !KIM---- update the average radiation for phenology----------------------------------!
      call update_rad_avg(cgrid)

   end if

   !---------------------------------------------------------------------------------------!
   !     At this point, all meteorologic driver data for the land surface model has been   !
   ! updated for the current timestep.  Perform the time space average for the output      !
   ! diagnostic.                                                                           !
   !---------------------------------------------------------------------------------------!
   call int_met_avg(cgrid)

   return
end subroutine radiate_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine will drive the distribution of radiation among crowns, snow layers,  !
! and soil.                                                                                !
!------------------------------------------------------------------------------------------!
subroutine sfcrad_ed(cosz, cosaoi, csite, maxcohort, rshort)

   use ed_state_vars        , only : sitetype        & ! structure  
                                   , patchtype       ! ! structure  
   use canopy_radiation_coms, only : Watts2Ein       & ! intent(in) 
                                   , crown_mod       ! ! intent(in)
   use grid_coms            , only : nzg             & ! intent(in)
                                   , nzs             ! ! intent(in)
   use soil_coms            , only : soil            & ! intent(in)
                                   , emisg           ! ! intent(in)
   use consts_coms          , only : stefan          ! ! intent(in)
   use ed_max_dims             , only : n_pft           ! ! intent(in)
   use allometry            , only : dbh2ca          ! ! function

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)  , target     :: csite
   real            , intent(in) :: rshort
   real            , intent(in) :: cosaoi
   real            , intent(in) :: cosz
   integer         , intent(in) :: maxcohort
   !----- Local variables. ----------------------------------------------------------------!
   type(patchtype) , pointer              :: cpatch
   integer         , dimension(maxcohort) :: pft_array
   integer                                :: il
   integer                                :: ipa,ico
   integer                                :: cohort_count
   integer                                :: k
   real                                   :: fcpct
   real                                   :: alg
   real                                   :: als
   real                                   :: rad
   real                                   :: fractrans
   real            , dimension(nzs)       :: fracabs
   real                                   :: absg
   real                                   :: algs
   real                                   :: downward_par_below_beam
   real                                   :: upward_par_above_beam
   real                                   :: downward_nir_below_beam
   real                                   :: upward_nir_above_beam
   real(kind=8)    , dimension(maxcohort) :: veg_temp_array
   real(kind=8)    , dimension(maxcohort) :: tai_array
   real(kind=8)    , dimension(maxcohort) :: CA_array
   real            , dimension(maxcohort) :: par_v_beam_array
   real            , dimension(maxcohort) :: rshort_v_beam_array
   real            , dimension(maxcohort) :: par_v_diffuse_array
   real            , dimension(maxcohort) :: rshort_v_diffuse_array
   real            , dimension(maxcohort) :: lw_v_surf_array
   real            , dimension(maxcohort) :: lw_v_incid_array
   real                                   :: downward_par_below_diffuse
   real                                   :: upward_par_above_diffuse
   real                                   :: downward_nir_below_diffuse
   real                                   :: upward_nir_above_diffuse 
   real                                   :: T_surface
   real                                   :: emissivity
   real                                   :: downward_lw_below_surf
   real                                   :: downward_lw_below_incid
   real                                   :: upward_lw_below_surf
   real                                   :: upward_lw_below_incid
   real                                   :: upward_lw_above_surf
   real                                   :: upward_lw_above_incid
   real                                   :: downward_rshort_below_beam
   real                                   :: downward_rshort_below_diffuse
   real                                   :: surface_absorbed_longwave_surf
   real                                   :: surface_absorbed_longwave_incid
   real                                   :: crown_area
   !----- Loop over the patches -----------------------------------------------------------!

   do ipa = 1,csite%npatches
      cpatch => csite%patch(ipa)
      

      !------------------------------------------------------------------------------------!
      !     Cohort_count is the number of cohorts that affect the radiation balance (i.e.  !
      ! those which are flagged as solvable.                                               !
      !------------------------------------------------------------------------------------!
      cohort_count = 0

      !----- Recalc the maximum photosynthetic rates next time around. --------------------!
      csite%old_stoma_data_max(1:n_pft,ipa)%recalc = 1

      !------------------------------------------------------------------------------------!
      !     Loop over cohorts.  Unusually, we here start at the shortest. Required by      !
      ! radiation schemes.                                                                 !
      !------------------------------------------------------------------------------------!
      do ico = cpatch%ncohorts,1,-1
         
         !----- Initialize values. --------------------------------------------------------!
         cpatch%par_v(ico)                 = 0.0
         cpatch%par_v_beam(ico)            = 0.0
         cpatch%par_v_diffuse(ico)         = 0.0
         
         cpatch%rshort_v(ico)              = 0.0
         cpatch%rshort_v_beam(ico)         = 0.0
         cpatch%rshort_v_diffuse(ico)      = 0.0
         
         cpatch%rlong_v(ico)               = 0.0
         cpatch%rlong_v_incid(ico)         = 0.0
         cpatch%rlong_v_surf(ico)          = 0.0
         cpatch%old_stoma_data(ico)%recalc = 1

         !------ Transfer information from linked lists to arrays. ------------------------!
         if (cpatch%solvable(ico)) then
            cohort_count                         = cohort_count + 1
            pft_array(cohort_count)              = cpatch%pft(ico)
            tai_array(cohort_count)              = dble(cpatch%lai(ico))                   &
                                                 + dble(cpatch%wai(ico))
            veg_temp_array(cohort_count)         = dble(cpatch%veg_temp(ico))
            rshort_v_beam_array(cohort_count)    = 0.0
            par_v_beam_array(cohort_count)       = 0.0
            rshort_v_diffuse_array(cohort_count) = 0.0
            par_v_diffuse_array(cohort_count)    = 0.0
            select case (crown_mod)
            case (0)
               CA_array(cohort_count) = 1.
            case (1)
               !----- Crown area allom from Dietze and Clark (2008). ----------------------!
               crown_area              = cpatch%nplant(ico)                                &
                                       * dbh2ca(cpatch%dbh(ico),cpatch%pft(ico))
               CA_array(cohort_count)  = min(1.d0,dble(crown_area))
            end select
         end if

      end do

      csite%rshort_s_diffuse(:,ipa) = 0.0
      csite%rshort_s_beam(:,ipa)    = 0.0

      !----- Soil water fraction. ---------------------------------------------------------!
      fcpct = csite%soil_water(nzg,ipa) / soil(csite%ntext_soil(nzg,ipa))%slmsts 

      !----- Finding the ground albedo as a function of soil water relative moisture. -----!
      alg = max(.14,.31-.34*fcpct)
      
      !------------------------------------------------------------------------------------!
      !     Deciding what is our surface temperature.  When the soil is exposed, then that !
      ! is the surface temperature.  Otherwise, we pick the temporary surface water or     !
      ! snow layer.                                                                        !
      !------------------------------------------------------------------------------------!
      rad  = 1.0
      algs = 1.0
      if(csite%nlev_sfcwater(ipa) == 0)then
         emissivity = emisg(csite%ntext_soil(nzg,ipa))
         T_surface  = csite%soil_tempk(nzg,ipa)
      else
         !---------------------------------------------------------------------------------!
         !      Sfcwater albedo ALS ranges from wet-soil value .14 for all-liquid to .5    !
         ! for all-ice.                                                                    !
         !---------------------------------------------------------------------------------!
         als = 0.5 - 0.36 * csite%sfcwater_fracliq(csite%nlev_sfcwater(ipa),ipa)

         !----- Fraction shortwave absorbed into sfcwater + soil. -------------------------!
         rad = 1.0 - als
         
         do k = csite%nlev_sfcwater(ipa),1,-1
            
            !------------------------------------------------------------------------------!
            !      Fractrans is fraction of shortwave entering each sfcwater layer that    !
            ! gets transmitted through that layer.                                         !
            !------------------------------------------------------------------------------!
            fractrans = exp(-20.0 * csite%sfcwater_depth(k,ipa))
            
            !------------------------------------------------------------------------------!
            !      Fracabs(k) is fraction of total incident shortwave (at top of top       !
            ! sfcwater layer) that is absorbed in each sfcwater layer.                     !
            !------------------------------------------------------------------------------!
            fracabs(k) = rad * (1.0 - fractrans)
            
            !------------------------------------------------------------------------------!
            !      Rad is fraction of total incident shortwave (at top of top sfcwater     !
            ! layer) that remains at bottom of current sfcwater layer.                     !
            !------------------------------------------------------------------------------!
            rad = rad * fractrans
            
            !------------------------------------------------------------------------------!
            !      Algs will ultimately be the albedo of the soil+sfcwater.  So subtract   !
            ! out whatever is getting absorbed by sfcwater.                                !
            !------------------------------------------------------------------------------!
            algs = algs - fracabs(k)
         end do

         !----- Long wave parameter if sfcwater exists. -----------------------------------!
         emissivity = 1.0
         T_surface  = csite%sfcwater_tempk(csite%nlev_sfcwater(ipa),ipa)
      end if
      
      csite%snowfac(ipa) = min(.99                                                         &
                              ,csite%total_snow_depth(ipa)/max(.001,csite%veg_height(ipa)))
      
      !------------------------------------------------------------------------------------!
      !     This is the fraction of below-canopy radiation that is absorbed by the ground. !
      !------------------------------------------------------------------------------------!
      absg = (1.0 - alg) * rad

      !----- Subtract off ground absorption to obtain the soil+sfcwater albedo. -----------!
      algs = algs - absg

      !----- Call the radiation parameterizations if there is vegetation. -----------------!
      if(cohort_count > 0)then

         !----- Long wave first. ----------------------------------------------------------!
         call lw_twostream(cohort_count,emissivity,T_surface,pft_array(1:cohort_count)     &
                          ,TAI_array(1:cohort_count),CA_array(1:cohort_count)              &
                          ,veg_temp_array(1:cohort_count),lw_v_surf_array(1:cohort_count)  &
                          ,lw_v_incid_array(1:cohort_count), downward_lw_below_surf        &
                          ,downward_lw_below_incid, upward_lw_below_surf                   &
                          ,upward_lw_below_incid, upward_lw_above_surf                     &
                          ,upward_lw_above_incid)

         !----- Upwelling long wave radiation at the top of the canopy. -------------------!
         csite%rlongup(ipa)      = upward_lw_above_surf
         csite%rlong_albedo(ipa) = upward_lw_above_incid
         
         !----- Long wave absorbed by either soil or sfcwater. ----------------------------!
         surface_absorbed_longwave_surf  = downward_lw_below_surf  - upward_lw_below_surf
         surface_absorbed_longwave_incid = downward_lw_below_incid - upward_lw_below_incid
         
         !----- Compute short wave if it is daytime. --------------------------------------!
         if (rshort > 0.5) then
            call sw_twostream_clump(algs,cosz,cosaoi,cohort_count                          &
                                   ,pft_array(1:cohort_count),TAI_array(1:cohort_count)    &
                                   ,CA_array(1:cohort_count)                               &
                                   ,par_v_beam_array(1:cohort_count)                       &
                                   ,par_v_diffuse_array(1:cohort_count)                    &
                                   ,rshort_v_beam_array(1:cohort_count)                    &
                                   ,rshort_v_diffuse_array(1:cohort_count)                 &
                                   ,downward_par_below_beam,downward_par_below_diffuse     &
                                   ,upward_par_above_beam,upward_par_above_diffuse         &
                                   ,downward_nir_below_beam,downward_nir_below_diffuse     &
                                   ,upward_nir_above_beam,upward_nir_above_diffuse)        

            !----- Below-canopy downwelling radiation. ------------------------------------!
            downward_rshort_below_beam    = downward_par_below_beam                        &
                                          + downward_nir_below_beam
            downward_rshort_below_diffuse = downward_par_below_diffuse                     &
                                          + downward_nir_below_diffuse

            !----- Soil+sfcwater+veg albedo (different for diffuse and beam radiation). ---!
            csite%albedo_beam(ipa)    = upward_par_above_beam    + upward_nir_above_beam
            csite%albedo_diffuse(ipa) = upward_par_above_diffuse + upward_nir_above_diffuse
         else

            !----- The code expects values for these, even when it is not daytime. --------!
            downward_rshort_below_beam    = 1.0
            downward_rshort_below_diffuse = 1.0
            csite%albedo_beam(ipa)        = algs
            csite%albedo_diffuse(ipa)     = algs

         end if

         !----- Absorption rates of PAR, rshort, and rlong of the vegetation. -------------!
         il = 0

         do ico = cpatch%ncohorts,1,-1
            if (cpatch%solvable(ico)) then
               il = il + 1
               cpatch%rshort_v_beam(ico)    = rshort_v_beam_array(il)
               cpatch%rshort_v_diffuse(ico) = rshort_v_diffuse_array(il)
               cpatch%par_v_beam(ico)       = par_v_beam_array(il) * Watts2Ein
               cpatch%par_v_diffuse(ico)    = par_v_diffuse_array(il) * Watts2Ein
               cpatch%rlong_v_surf(ico)     = lw_v_surf_array(il)
               cpatch%rlong_v_incid(ico)    = lw_v_incid_array(il)
            end if
         end do

      else
         
         !----- This is the case where there is no vegetation. ----------------------------!
         downward_rshort_below_beam      = 1.0
         downward_rshort_below_diffuse   = 1.0
         surface_absorbed_longwave_surf  = - emissivity * stefan * T_surface**4
         surface_absorbed_longwave_incid = emissivity
         csite%albedo_beam(ipa)          = algs
         csite%albedo_diffuse(ipa)       = algs
         csite%rlongup(ipa)              = - surface_absorbed_longwave_surf
         csite%rlong_albedo(ipa)         = 1.0 - surface_absorbed_longwave_incid
      
      endif
      
      !----- Absorption rate of short wave by the soil. -----------------------------------!
      csite%rshort_g_beam(ipa)    = downward_rshort_below_beam    * absg
      csite%rshort_g_diffuse(ipa) = downward_rshort_below_diffuse * absg

      !----- Absorption rate of short wave by the surface water. --------------------------!
      do k=1,csite%nlev_sfcwater(ipa)
         csite%rshort_s_beam(k,ipa)    = downward_rshort_below_beam    * fracabs(k)
         csite%rshort_s_diffuse(k,ipa) = downward_rshort_below_diffuse * fracabs(k)
      end do

      !----- Long wave absorption rate at the surface. ------------------------------------!
      if (csite%nlev_sfcwater(ipa) == 0) then
         csite%rlong_s_surf(ipa)  = 0.0
         csite%rlong_s_incid(ipa) = 0.0
         csite%rlong_g_surf(ipa)  = surface_absorbed_longwave_surf
         csite%rlong_g_incid(ipa) = surface_absorbed_longwave_incid
      else
         csite%rlong_s_surf(ipa)  = surface_absorbed_longwave_surf
         csite%rlong_s_incid(ipa) = surface_absorbed_longwave_incid
         csite%rlong_g_surf(ipa)  = 0.0
         csite%rlong_g_incid(ipa) = 0.0
      end if
      
   end do
   return
end subroutine sfcrad_ed
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine solar_zenith(cgrid)
   use ed_misc_coms , only : current_time ! ! intent(in)
   use consts_coms  , only : pio1808      & ! intent(in)
                           , twopi8       & ! intent(in)
                           , hr_sec8      ! ! intent(in)
   use ed_state_vars, only : edtype       ! ! intent(in)

   implicit none
   !----- Argument. -----------------------------------------------------------------------!
   type(edtype), target   :: cgrid  ! Current ED grid
   !----- Local variables. ----------------------------------------------------------------!
   integer                :: ipy    ! Polygon counter
   integer                :: jday   ! Day of year ("Julian" day)
   real(kind=8)           :: declin ! Declination
   real(kind=8)           :: sdec   ! Sine of declination
   real(kind=8)           :: cdec   ! Cosine of declination
   real(kind=8)           :: dayhr  ! Hour of day 
   real(kind=8)           :: radlat ! Latitude in radians
   real(kind=8)           :: cslcsd ! Cosine of latitude times cosine of declination
   real(kind=8)           :: snlsnd ! Sine of latitude times sine of declination
   real(kind=8)           :: dayhrr ! Hour of day in radians
   real(kind=8)           :: hrangl ! Hour angle
   !----- External functions. -------------------------------------------------------------!
   integer     , external :: julday  ! Function to find day of year ("Julian" day)
   real        , external :: sngloff ! Function to safely convert double to single prec.
   !---------------------------------------------------------------------------------------!

   jday  = julday(current_time%month, current_time%date, current_time%year)

   !----- sdec - sine of declination, cdec - cosine of declination. -----------------------!
   declin = -2.35d1 * cos(twopi8 / 3.65d2 * dble(jday + 9)) * pio1808
   sdec   = dsin(declin)
   cdec   = dcos(declin)

   !----- Find the hour angle, then get cosine of zenith angle. ---------------------------!
   dayhr = dble(current_time%time) / hr_sec8

   do ipy = 1, cgrid%npolygons

      radlat = dble(cgrid%lat(ipy)) * pio1808
      if (radlat == declin) radlat = radlat + 1.d-5
      cslcsd = dcos(radlat) * cdec
      snlsnd = dsin(radlat) * sdec
      dayhrr = dmod(dayhr+dble(cgrid%lon(ipy))/1.5d1+2.4d1,2.4d1)
      hrangl = 1.5d1 * (dayhrr - 1.2d1) * pio1808

      cgrid%cosz(ipy) = sngloff(snlsnd + cslcsd * dcos(hrangl),1.d-20)

   end do

   return
end subroutine solar_zenith
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed2land_radiation(cpoly,isi)
   use ed_state_vars  , only : polygontype & ! intent(in)
                             , sitetype    & ! intent(in)
                             , patchtype   ! ! intent(in)
   implicit none
   !----- Argument. -----------------------------------------------------------------------!
   type(polygontype) , target     :: cpoly
   integer           , intent(in) :: isi
   !----- Local variables. ----------------------------------------------------------------!
   type(sitetype)    ,  pointer   :: csite
   integer                        :: ipa
   !---------------------------------------------------------------------------------------!


   !----- Initialize arrays to zero. ------------------------------------------------------!
   cpoly%albedo_beam(isi)    = 0.0
   cpoly%albedo_diffuse(isi) = 0.0
   cpoly%rlongup(isi)        = 0.0
   cpoly%rlong_albedo(isi)   = 0.0


   csite => cpoly%site(isi)
   !----- Loop over patches. --------------------------------------------------------------!
   do ipa = 1,csite%npatches
      !------------------------------------------------------------------------------------!
      !     Compute cell-level albedo and upward longwave, weighting patches by fractional !
      ! area.                                                                              !
      !------------------------------------------------------------------------------------!
      cpoly%albedo_beam(isi)    = cpoly%albedo_beam(isi)                                   &
                                + csite%albedo_beam(ipa) * csite%area(ipa)
      cpoly%albedo_diffuse(isi) = cpoly%albedo_diffuse(isi)                                &
                                + csite%albedo_diffuse(ipa) * csite%area(ipa)
      cpoly%rlongup(isi)        = cpoly%rlongup(isi) + csite%rlongup(ipa) * csite%area(ipa)
      cpoly%rlong_albedo(isi)   = cpoly%rlong_albedo(isi)                                  &
                                + csite%rlong_albedo(ipa) * csite%area(ipa)
   end do

   return
end subroutine ed2land_radiation
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine scale_ed_radiation(rshort, rshort_diffuse, rlong, csite)

   use ed_state_vars , only : sitetype  & ! intent(in)
                            , patchtype ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(sitetype)  , target     :: csite
   real            , intent(in) :: rshort
   real            , intent(in) :: rshort_diffuse
   real            , intent(in) :: rlong
   !----- Local variables. ----------------------------------------------------------------!
   type(patchtype) , pointer    :: cpatch
   integer                      :: ipa,ico, k
   real                         :: beam_radiation
   !----- This should be false unless you really want to turn off radiation. --------------!
   logical         , parameter  :: skip_rad = .false.
   !----- External functions. -------------------------------------------------------------!
   real            , external   :: sngloff
   !---------------------------------------------------------------------------------------!

   if (skip_rad) then
      do ipa = 1, csite%npatches
         cpatch => csite%patch(ipa)
         do ico = 1, cpatch%ncohorts
            if (cpatch%solvable(ico)) then
               cpatch%rshort_v_beam(ico)    = 0.
               cpatch%rshort_v_diffuse(ico) = 0.
               cpatch%rshort_v(ico)         = 0.
               cpatch%par_v_beam(ico)       = 0.
               cpatch%par_v_diffuse(ico)    = 0.
               cpatch%par_v(ico)            = 0.
               cpatch%rlong_v_incid(ico)    = 0.
               cpatch%rlong_v_surf(ico)     = 0.
               cpatch%rlong_v(ico)          = 0.
            end if
         end do
         
         csite%rshort_g_beam(ipa)    = 0.
         csite%rshort_g_diffuse(ipa) = 0.
         csite%rshort_g(ipa)         = 0.
         !----- Absorption rate of short wave by the surface water. -----------------------!
         do k=1,csite%nlev_sfcwater(ipa)
            csite%rshort_s_beam(k,ipa)    = 0.
            csite%rshort_s_diffuse(k,ipa) = 0.
            csite%rshort_s(k,ipa)         = 0.
         end do
         csite%rlong_s_incid(ipa) = 0.
         csite%rlong_g_incid(ipa) = 0.
         csite%rlong_s_surf(ipa)  = 0.
         csite%rlong_g_surf(ipa)  = 0.
      
         csite%rlong_s(ipa)       = 0.
         csite%rlong_g(ipa)       = 0.
      end do
      return
   end if

   beam_radiation = rshort - rshort_diffuse

   do ipa = 1,csite%npatches

      cpatch => csite%patch(ipa)
      do ico = 1,cpatch%ncohorts
         
         if (cpatch%solvable(ico)) then
            cpatch%rshort_v_beam(ico)    = cpatch%rshort_v_beam(ico)    * beam_radiation
            cpatch%rshort_v_diffuse(ico) = cpatch%rshort_v_diffuse(ico) * rshort_diffuse
            cpatch%rshort_v(ico)         = cpatch%rshort_v_beam(ico)                       &
                                         + cpatch%rshort_v_diffuse(ico)
            
            cpatch%par_v_beam(ico)       = cpatch%par_v_beam(ico) * beam_radiation
            cpatch%par_v_diffuse(ico)    = cpatch%par_v_diffuse(ico) * rshort_diffuse
            cpatch%par_v(ico)            = cpatch%par_v_beam(ico)                          &
                                         + cpatch%par_v_diffuse(ico)
            
            cpatch%rlong_v_incid(ico)    = cpatch%rlong_v_incid(ico) * rlong
            
            cpatch%rlong_v(ico)          = cpatch%rlong_v_incid(ico)                       &
                                         + cpatch%rlong_v_surf(ico)
         end if
      end do

      csite%rshort_g_beam(ipa)    = csite%rshort_g_beam(ipa)    * beam_radiation
      csite%rshort_g_diffuse(ipa) = csite%rshort_g_diffuse(ipa) * rshort_diffuse
      csite%rshort_g(ipa)         = csite%rshort_g_beam(ipa) + csite%rshort_g_diffuse(ipa)
      
      !----- Absorption rate of short wave by the surface water. --------------------------!
      do k=1,csite%nlev_sfcwater(ipa)
         csite%rshort_s_beam(k,ipa)    = csite%rshort_s_beam(k,ipa) * beam_radiation
         csite%rshort_s_diffuse(k,ipa) = csite%rshort_s_diffuse(k,ipa) * rshort_diffuse
         csite%rshort_s(k,ipa)         = csite%rshort_s_beam(k,ipa)                        &
                                       + csite%rshort_s_diffuse(k,ipa)
      end do

      csite%rlong_s_incid(ipa) = csite%rlong_s_incid(ipa) * rlong
      csite%rlong_g_incid(ipa) = csite%rlong_g_incid(ipa) * rlong
      
      csite%rlong_s(ipa)       = csite%rlong_s_surf(ipa) + csite%rlong_s_incid(ipa)
      csite%rlong_g(ipa)       = csite%rlong_g_surf(ipa) + csite%rlong_g_incid(ipa)
   end do

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






!==========================================================================================!
!==========================================================================================!
subroutine short2diff(swdown,sunang,radvdc)

   implicit none
   !------------------------------------------------------------------
   !---
   !               radiation radive code to use the downward sw at
   ! bottom
   !               and the formulation to estimate radvbc,radvdc,
   ! radndc,
   !               radndc
   !  This subroutine was adapted from the sib2 model (sellars et al.)
   !  
   !------------------------------------------------------------------
   !---

   real swdown
   real sunang, stemp
   real radvdc
   real c1,c2,c3,c4,c5,cloud,difrat
   real vnrat

   ! Arguments:
   ! nsib:
   ! swdown: surface incident shortwave radiation
   ! sunang: cosine of solar zenith angle


   c1 = 580.
   c2 = 464.
   c3 = 499.
   c4 = 963.
   c5 = 1160.


   sunang = max( 0.001 , sunang )
   stemp = swdown
   stemp = max(stemp,0.01 )
   cloud = (c5 * sunang - stemp) / (c4 * sunang)
   cloud = max(cloud,0.)
   cloud = min(cloud,1.)
   !         cloud = max(0.58,cloud)
   
   !z  use the real clouds here!
   !         cloud = cldtot(i)
   !         CLOUD = AMAX1(CLOUD,0.)
   !         CLOUD = AMIN1(CLOUD,1.)
   
   difrat = 0.0604 / ( sunang -0.0223 ) + 0.0683
   if ( difrat .lt. 0. ) difrat = 0.
   if ( difrat .gt. 1. ) difrat = 1.
   
   difrat = difrat + ( 1. - difrat ) * cloud
   vnrat = ( c1 - cloud*c2 ) / ( ( c1 - cloud*c3 ) + ( c1 - cloud*c2 ))
   
   radvdc = difrat*vnrat*stemp
   
   return
end subroutine short2diff



