subroutine radiate_driver_ar(cgrid)
  
  use misc_coms, only: current_time, radfrq, dtlsm
  use ed_state_vars, only: edtype,polygontype,sitetype,patchtype
  use canopy_radiation_coms, only: visible_fraction_dir, visible_fraction_dif, rlong_min
  use consts_coms, only: pio180

  implicit none

  type(edtype),target :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer :: csite
  type(patchtype),pointer :: cpatch
  real :: total_beam_radiation
  integer :: maxcohort
  integer :: ipy,isi,ipa
  real :: hrangl

  ! Check whether it is time to update radiative fluxes and heating rates
  
  if (mod(current_time%time + .001,radfrq) < dtlsm) then
     
     ! Compute solar zenith angle [cosz]

     call solar_zenith_ar(cgrid)

     ! Loop over polygons

     do ipy = 1,cgrid%npolygons

        cpoly => cgrid%polygon(ipy)

        do isi = 1,cpoly%nsites

           csite => cpoly%site(isi)

           
           if(cpoly%met(isi)%rlong < rlong_min ) then
              print*,"STRANGE DATA",cpoly%met(isi)%rlong,ipy,isi,int(real(isi)/4.0),rlong_min
              print*,cpoly%met(isi)
              call fatal_error('Rlong is too low!','radiate_driver_ar'&
                              &,'radiate_driver.f90')
           endif



           ! Compute the visible fraction of diffuse and beam radiation
           ! needed by the radiative transfer routine.

           total_beam_radiation = cpoly%met(isi)%nir_beam +   &
                cpoly%met(isi)%par_beam
           
           if(total_beam_radiation > 0.0)then
              visible_fraction_dir = cpoly%met(isi)%par_beam /   &
                   total_beam_radiation
           else
              visible_fraction_dir = 0.5
           endif
           
           if(cpoly%met(isi)%rshort_diffuse > 0.0)then
              visible_fraction_dif = cpoly%met(isi)%par_diffuse /   &
                   cpoly%met(isi)%rshort_diffuse
           else
              visible_fraction_dif = 0.5
           endif
           
           ! Update angle of incidence
           hrangl = 15.0 * (mod(current_time%time + cgrid%lon(ipy) / 15.0 +   &
                24.0, 24.0) - 12.0) * pio180
           
           call angle_of_incid(cpoly%cosaoi(isi), cgrid%cosz(ipy), hrangl,   &
                cpoly%slope(isi) * pio180, cpoly%aspect(isi) * pio180)
           
           
!           print*,"AR ",ipy,cpoly%cosaoi(isi), cgrid%cosz(ipy), hrangl,   &
!                cpoly%slope(isi),cpoly%aspect(isi),cpoly%met(isi)%rlong, &
!                cpoly%met(isi)%par_beam,cpoly%met(isi)%nir_beam,cpoly%met(isi)%par_diffuse, &
!                cpoly%met(isi)%rshort_diffuse


           ! Loop over subgrid-scale patches.
           ! THESE ROUTINES CAN BE DONE AS ARRAYS

           maxcohort = 1
           do ipa = 1,csite%npatches
              cpatch=>csite%patch(ipa)
              if ( cpatch%ncohorts>maxcohort) then
                 maxcohort = cpatch%ncohorts
              endif
           enddo


           ! Get unnormalized radiative transfer information.
           call sfcrad_ed_ar(cgrid%cosz(ipy), cpoly%cosaoi(isi), csite, &
                maxcohort, cpoly%met(isi)%rshort)

           ! Normalize the absorbed radiations.
           call scale_ed_radiation_ar(cpoly%met(isi)%rshort,   &
                cpoly%met(isi)%rshort_diffuse, cpoly%met(isi)%rlong, csite)
                      
           ! Update all albedos and rlongup.
           call ed2land_radiation_ar(cpoly,isi)
     
  
        enddo
        
     enddo

  endif
       
  ! At this point, all meteorologic driver data for the land
  ! surface model has been updated for the current timestep
  ! Perform the time space average for the output diagnostic
  call int_met_avg(cgrid)

  return
end subroutine radiate_driver_ar

!==========================================================================

subroutine sfcrad_ed_ar(cosz, cosaoi, csite, maxcohort, rshort)

  use ed_state_vars,only:sitetype,patchtype
  use canopy_radiation_coms, only: lai_min, Watts2Ein
  use grid_coms, only: nzg, nzs
  use soil_coms, only: soil, emisg
  use consts_coms, only: stefan
  use max_dims, only: n_pft
  use pft_coms,only:sla
  use allometry, only : dbh2ca

  implicit none

  real, intent(in) :: rshort
  real, intent(in) :: cosaoi
  integer, intent(in) :: maxcohort
  type(sitetype),target   :: csite
  type(patchtype),pointer :: cpatch
  integer :: ipa,ico
  integer :: cohort_count
  real :: fcpct
  real :: alg
  real :: als
  real :: rad
  integer :: k
  real :: fractrans
  real, dimension(nzs) :: fracabs
  real :: absg
  real :: algs
  real :: cosz
  real, dimension(maxcohort) :: par_v_beam_array
  real, dimension(maxcohort) :: rshort_v_beam_array
  real, dimension(maxcohort) :: par_v_diffuse_array
  real, dimension(maxcohort) :: rshort_v_diffuse_array
  real(kind=8), dimension(maxcohort) :: veg_temp_array
  real(kind=8), dimension(maxcohort) ::  lai_array
  real(kind=8), dimension(maxcohort) ::  CA_array
  integer, dimension(maxcohort) :: pft_array
  real :: downward_par_below_beam
  real :: upward_par_above_beam
  real :: downward_nir_below_beam
  real :: upward_nir_above_beam
  integer :: il
  real :: downward_par_below_diffuse
  real :: upward_par_above_diffuse
  real :: downward_nir_below_diffuse
  real :: upward_nir_above_diffuse 
  real :: T_surface
  real :: emissivity
  real, dimension(maxcohort) :: lw_v_surf_array
  real, dimension(maxcohort) :: lw_v_incid_array
  real :: downward_lw_below_surf
  real :: downward_lw_below_incid
  real :: upward_lw_below_surf
  real :: upward_lw_below_incid
  real :: upward_lw_above_surf
  real :: upward_lw_above_incid
  real :: downward_rshort_below_beam
  real :: downward_rshort_below_diffuse
  real :: surface_absorbed_longwave_surf
  real :: surface_absorbed_longwave_incid

  ! Loop over the patches

  do ipa = 1,csite%npatches

     cpatch => csite%patch(ipa)
     
     ! Calculate the total snow depth
     if(csite%nlev_sfcwater(ipa) == 0)then
        csite%total_snow_depth(ipa) = 0.0
     else
        csite%total_snow_depth(ipa) = sum(csite%sfcwater_depth(1:csite%nlev_sfcwater(ipa),ipa))
     endif

     ! cohort_count is the number of cohorts with leaves that are not 
     ! covered by snow.
     cohort_count = 0
     
     ! recalc the maximum photosynthetic rates next time around.
     csite%old_stoma_data_max(1:n_pft,ipa)%recalc = 1

     ! loop over cohorts
     do ico = cpatch%ncohorts,1,-1

        ! Unusually, we here start at the shortest. Required by radiation schemes.
        
        ! initialize values
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

        ! transfer information from linked lists to arrays
        if(cpatch%lai(ico) > lai_min .and. cpatch%hite(ico) > &
             csite%total_snow_depth(ipa))then
           cohort_count = cohort_count + 1
           pft_array(cohort_count) = cpatch%pft(ico)
           lai_array(cohort_count) = dble(cpatch%lai(ico))
           !! crown area allom from Dietze and Clark 2008
           CA_array(cohort_count)  = dble(min(1.0,&
                cpatch%nplant(ico)*dbh2ca(cpatch%dbh(ico),cpatch%pft(ico))))
           veg_temp_array(cohort_count) = dble(cpatch%veg_temp(ico))
           rshort_v_beam_array(cohort_count) = 0.0
           par_v_beam_array(cohort_count) = 0.0
           rshort_v_diffuse_array(cohort_count) = 0.0
           par_v_diffuse_array(cohort_count) = 0.0
           ! long wave doesn't need to be initialized.
        endif

     enddo
     csite%rshort_s_diffuse(:,ipa) = 0.0
     csite%rshort_s_beam(:,ipa) = 0.0

     fcpct = real(csite%soil_water(nzg,ipa) / dble(soil(csite%ntext_soil(nzg,ipa))%slmsts)) ! soil water fraction

     if (fcpct > .5) then
        alg = .14                ! ground albedo
     else
        alg = .31 - .34 * fcpct  ! ground albedo
     endif
     
     rad = 1.0
     algs = 1.0
     if(csite%nlev_sfcwater(ipa) == 0)then
        emissivity = emisg(csite%ntext_soil(nzg,ipa))
        T_surface = csite%soil_tempk(nzg,ipa)
     else

        ! Sfcwater albedo ALS ranges from wet-soil value .14 for all-liquid
        ! to .5 for all-ice
        als = 0.5 - 0.36 * csite%sfcwater_fracliq(csite%nlev_sfcwater(ipa),ipa)

        rad = 1.0 - als    ! fraction shortwave absorbed into sfcwater + soil
        
        do k = csite%nlev_sfcwater(ipa),1,-1
           
           ! fractrans is fraction of shortwave entering each sfcwater layer that
           ! gets transmitted through that layer
           fractrans = exp(-20.0 * csite%sfcwater_depth(k,ipa))
           
           ! fracabs(k) is fraction of total incident shortwave (at top of 
           ! top sfcwater layer) that is absorbed in each sfcwater layer
           fracabs(k) = rad * (1.0 - fractrans)
           
           ! rad is fraction of total incident shortwave (at top of top sfcwater 
           ! layer) that remains at bottom of current sfcwater layer
           rad = rad * fractrans
           
           ! algs will ultimately be the albedo of the soil+sfcwater.  So 
           ! subtract out whatever is getting absorbed by sfcwater.
           algs = algs - fracabs(k)
        enddo
        ! long wave parameter if sfcwater exists
        emissivity = 1.0
        T_surface = csite%sfcwater_tempk(csite%nlev_sfcwater(ipa),ipa)
     endif
     
     csite%snowfac(ipa) = min(.99, csite%total_snow_depth(ipa) / max(.001,csite%veg_height(ipa)))
     
     ! This is the fraction of below-canopy radiation that is absorbed by 
     ! the ground
     absg = (1.0 - alg) * rad
     
     ! Subtract off ground absorption to obtain the soil+sfcwater albedo.
     algs = algs - absg
     
     ! Call the radiation parameterizations if there is vegetation
     if(cohort_count > 0)then
        
        ! Long wave first.
        call lw_twostream(cohort_count, emissivity, T_surface,  &
             pft_array(1:cohort_count), lai_array(1:cohort_count), & 
             CA_array(1:cohort_count),  &
             veg_temp_array(1:cohort_count), lw_v_surf_array(1:cohort_count),  &
             lw_v_incid_array(1:cohort_count), downward_lw_below_surf,  &
             downward_lw_below_incid, upward_lw_below_surf,   &
             upward_lw_below_incid, upward_lw_above_surf, upward_lw_above_incid)
        
        ! Upwelling long wave radiation at the top of the canopy
        csite%rlongup(ipa) = upward_lw_above_surf
        csite%rlong_albedo(ipa) = upward_lw_above_incid
        
        ! long wave absorbed by either soil or sfcwater
        surface_absorbed_longwave_surf = downward_lw_below_surf -   &
             upward_lw_below_surf
        surface_absorbed_longwave_incid = downward_lw_below_incid -   &
             upward_lw_below_incid
        
        ! Compute short wave if it is daytime.
        if(rshort > 0.5)then
           !     if(cosz > 0.03)then
           ! call the two-stream approximation
           call sw_twostream_clump(algs,  &
                cosz,  &
                cosaoi, &
                csite%lai(ipa),  &
                cohort_count,   &
                pft_array(1:cohort_count),  &
                lai_array(1:cohort_count),   &
                CA_array(1:cohort_count),    &
                par_v_beam_array(1:cohort_count),  &
                par_v_diffuse_array(1:cohort_count),  &
                rshort_v_beam_array(1:cohort_count), &
                rshort_v_diffuse_array(1:cohort_count), &
                downward_par_below_beam,  &
                downward_par_below_diffuse,  &
                upward_par_above_beam,  &
                upward_par_above_diffuse,  &
                downward_nir_below_beam,  &
                downward_nir_below_diffuse,  &
                upward_nir_above_beam,   &
                upward_nir_above_diffuse)        
        
           ! below-canopy downwelling radiation
           downward_rshort_below_beam = downward_par_below_beam +   &
                downward_nir_below_beam
           downward_rshort_below_diffuse = downward_par_below_diffuse +   &
                downward_nir_below_diffuse

           ! soil+sfcwater+veg albedo (different for diffuse and beam radiation)
           csite%albedo_beam(ipa) = upward_par_above_beam + upward_nir_above_beam
           csite%albedo_diffuse(ipa) = upward_par_above_diffuse + upward_nir_above_diffuse
           
        else

           ! code expects values for these, even if it is not day.
           downward_rshort_below_beam = 1.0
           downward_rshort_below_diffuse = 1.0
           csite%albedo_beam(ipa) = algs
           csite%albedo_diffuse(ipa) = algs
           
        endif
     
        ! Absorption rates of PAR, rshort, and rlong of the vegetation
        il = 0
        do ico = cpatch%ncohorts,1,-1

           if(cpatch%lai(ico) > lai_min .and. cpatch%hite(ico) > csite%total_snow_depth(ipa))then
              il = il + 1
              cpatch%rshort_v_beam(ico) = rshort_v_beam_array(il)
              cpatch%rshort_v_diffuse(ico) = rshort_v_diffuse_array(il)
              cpatch%par_v_beam(ico) = par_v_beam_array(il) * Watts2Ein
              cpatch%par_v_diffuse(ico) = par_v_diffuse_array(il) * Watts2Ein
              cpatch%rlong_v_surf(ico) = lw_v_surf_array(il)
              cpatch%rlong_v_incid(ico) = lw_v_incid_array(il)

              if (cpatch%rlong_v_surf(ico).ne.cpatch%rlong_v_surf(ico)) then
                 print*,"LWSURF NAN",ico,il,cpatch%ncohorts,cohort_count
                 print*,""
                 print*,lw_v_surf_array
                 call fatal_error('NaN found in LWSURF' &
                                 &,'sfcrad_ed_ar','radiate_driver.f90')
              endif

              if (cpatch%rlong_v_incid(ico).ne.cpatch%rlong_v_incid(ico)) then
                 print*,"LWINCID NAN",ico,il,cpatch%ncohorts,cohort_count
                 print*,""
                 print*,lw_v_incid_array
                 call fatal_error('NaN found in LW_INCID' &
                                 &,'sfcrad_ed_ar','radiate_driver.f90')
                 
              endif
           endif
        enddo

     else
        
        ! This is the case where there is no vegetation
        downward_rshort_below_beam = 1.0
        downward_rshort_below_diffuse = 1.0
        surface_absorbed_longwave_surf = - emissivity * stefan * T_surface**4
        surface_absorbed_longwave_incid = emissivity
        csite%albedo_beam(ipa) = algs
        csite%albedo_diffuse(ipa) = algs
        csite%rlongup(ipa) = - surface_absorbed_longwave_surf
        csite%rlong_albedo(ipa) = 1.0 - surface_absorbed_longwave_incid
     
     endif
     
     ! Absorption rate of short wave by the soil
     csite%rshort_g_beam(ipa) = downward_rshort_below_beam * absg
     csite%rshort_g_diffuse(ipa) = downward_rshort_below_diffuse * absg

     ! Absorption rate of short wave by the surface water
     do k=1,csite%nlev_sfcwater(ipa)
        csite%rshort_s_beam(k,ipa) = downward_rshort_below_beam * fracabs(k)
        csite%rshort_s_diffuse(k,ipa) = downward_rshort_below_diffuse * fracabs(k)
     enddo

     ! Long wave absorption rate at the surface
     if(csite%nlev_sfcwater(ipa) == 0)then
        csite%rlong_s_surf(ipa) = 0.0
        csite%rlong_s_incid(ipa) = 0.0
        csite%rlong_g_surf(ipa) = surface_absorbed_longwave_surf
        csite%rlong_g_incid(ipa) = surface_absorbed_longwave_incid
     else
        csite%rlong_s_surf(ipa) = surface_absorbed_longwave_surf
        csite%rlong_s_incid(ipa) = surface_absorbed_longwave_incid
        csite%rlong_g_surf(ipa) = 0.0
        csite%rlong_g_incid(ipa) = 0.0
     endif
     
  enddo

  return
end subroutine sfcrad_ed_ar

!******************************************************************************

subroutine solar_zenith_ar(cgrid)

  use misc_coms, only: current_time
  use consts_coms, only: pio180,twopi
  use ed_state_vars, only: edtype

  implicit none
  
  type(edtype),target       :: cgrid
  integer :: ipy
  integer :: jday
  integer, external :: julday
  real :: declin
  real :: sdec
  real :: cdec
  real :: d0
  real :: d02
  real :: solfac
  real :: dayhr
  real :: radlat
  real :: cslcsd
  real :: snlsnd
  real :: dayhrr
  real :: hrangl


  jday  = julday(current_time%month, current_time%date, current_time%year)

  ! sdec - sine of declination, cdec - cosine of declination
  
  declin = -23.5 * cos(twopi / 365. * real(jday + 9)) * pio180
  sdec = sin(declin)
  cdec = cos(declin)

  ! Find the factor, solfac, to multiply the solar constant to correct
  ! for Earth's varying distance to the sun.
  
  d0 = twopi * float(jday-1) / 365.
  d02 = d0 * 2.
  solfac = 1.000110 + 0.034221 * cos (d0) + 0.001280 * sin(d0)  &
       + 0.000719 * cos(d02) + 0.000077 * sin(d02)

  ! Find the hour angle, then get cosine of zenith angle.
  
  dayhr = current_time%time / 3600.

  do ipy = 1,cgrid%npolygons

     radlat = cgrid%lat(ipy) * pio180
     if (radlat == declin) radlat = radlat + 1.e-5
     cslcsd = cos(radlat) * cdec
     snlsnd = sin(radlat) * sdec
     dayhrr = mod(dayhr+cgrid%lon(ipy)/15.+24.,24.)
     hrangl = 15. * (dayhrr - 12.) * pio180
     
     cgrid%cosz(ipy) = snlsnd + cslcsd * cos(hrangl)

  enddo

  return
end subroutine solar_zenith_ar

!=================================================================

subroutine ed2land_radiation_ar(cpoly,isi)
  
  use ed_state_vars, only: polygontype,sitetype,patchtype

  implicit none
  
  type(polygontype),target :: cpoly
  type(sitetype),  pointer :: csite
  integer :: isi,ipa

  ! initialize arrays to zero
  cpoly%albedo_beam(isi) = 0.0
  cpoly%albedo_diffuse(isi) = 0.0
  cpoly%rlongup(isi) = 0.0
  cpoly%rlong_albedo(isi) = 0.0

  csite => cpoly%site(isi)

  ! loop over patches
  do ipa = 1,csite%npatches
     ! compute cell-level albedo and upward longwave, weighting patches 
     ! by fractional area
     cpoly%albedo_beam(isi) = cpoly%albedo_beam(isi) + csite%albedo_beam(ipa) * csite%area(ipa)
     cpoly%albedo_diffuse(isi) = cpoly%albedo_diffuse(isi) + csite%albedo_diffuse(ipa) * csite%area(ipa)
     cpoly%rlongup(isi) = cpoly%rlongup(isi) + csite%rlongup(ipa) * csite%area(ipa)
     cpoly%rlong_albedo(isi) = cpoly%rlong_albedo(isi) + csite%rlong_albedo(ipa) * csite%area(ipa)
  enddo

  ! Update the radiation diagnostic quantities
  !  cs%avg_rlong     = cs%avg_rlong    + cpoly%metinput%rlong          * cs%area
  !  cs%avg_rlongup   = cs%avg_rlongup  + cs%rlongup           * cs%area
  
  !  cs%avg_rshort    = cs%avg_rshort   + cpoly%metinput%rshort         * cs%area
  
  !  cs%avg_rshortup  = cs%avg_rshortup + &
  !       (cpoly%metinput%rshort-cpoly%metinput%rshort_diffuse)*cs%albedo_beam*cs%area &
  !       + cpoly%metinput%rshort_diffuse*cs%albedo_diffuse  * cs%area
  
  !  cs%avg_rshortd   = cs%avg_rshortd  + cpoly%metinput%rshort_diffuse * cs%area
  

  return
end subroutine ed2land_radiation_ar

!==========================================================================

subroutine scale_ed_radiation_ar(rshort, rshort_diffuse, rlong, csite)

  use ed_state_vars,only:sitetype,patchtype
  use canopy_radiation_coms, only: lai_min

  implicit none
  type(sitetype),target :: csite
  type(patchtype), pointer :: cpatch
  integer :: ipa,ico
  real :: beam_radiation
  integer :: k
  real, intent(in) :: rshort
  real, intent(in) :: rshort_diffuse
  real, intent(in) :: rlong

  beam_radiation = rshort - rshort_diffuse

  do ipa = 1,csite%npatches

     cpatch => csite%patch(ipa)

     do ico = 1,cpatch%ncohorts
        
        if(cpatch%lai(ico) > lai_min .and. cpatch%hite(ico) > csite%total_snow_depth(ipa))then
           cpatch%rshort_v_beam(ico) = cpatch%rshort_v_beam(ico) * beam_radiation
           cpatch%rshort_v_diffuse(ico) = cpatch%rshort_v_diffuse(ico) * rshort_diffuse
           cpatch%rshort_v(ico) = cpatch%rshort_v_beam(ico) + cpatch%rshort_v_diffuse(ico)
           
           cpatch%par_v_beam(ico) = cpatch%par_v_beam(ico) * beam_radiation
           cpatch%par_v_diffuse(ico) = cpatch%par_v_diffuse(ico) * rshort_diffuse
           cpatch%par_v(ico) = cpatch%par_v_beam(ico) + cpatch%par_v_diffuse(ico)
           
           cpatch%rlong_v_incid(ico) = cpatch%rlong_v_incid(ico) * rlong
           cpatch%rlong_v(ico) = cpatch%rlong_v_incid(ico) + cpatch%rlong_v_surf(ico)

           if (cpatch%rlong_v(ico) .ne. cpatch%rlong_v(ico)) then
              print*,"cpatch%long_v(ico) is nan"
              print*,cpatch%rlong_v(ico),cpatch%rlong_v_incid(ico),cpatch%rlong_v_surf(ico),rlong
              print*,"ico:",ico
              call fatal_error('Rlong_v is NaN','scale_ed_radiation_ar' &
                              &,'radiate_driver.f90')
           end if

        endif

        if (cpatch%rlong_v(ico) .ne. cpatch%rlong_v(ico)) then
           print*,"cpatch%rlong_v(ico) is nan but not calculated"
           print*,cpatch%rlong_v(ico),cpatch%rlong_v_incid(ico),cpatch%rlong_v_surf(ico)
           call fatal_error('cpatch%rlong_v(ico) is nan but not calculated!'&
                           &,'scale_ed_radiation_ar','radiate_driver.f90')
        end if


     enddo
     
     csite%rshort_g_beam(ipa) = csite%rshort_g_beam(ipa) * beam_radiation
     csite%rshort_g_diffuse(ipa) = csite%rshort_g_diffuse(ipa) * rshort_diffuse
     csite%rshort_g(ipa) = csite%rshort_g_beam(ipa) + csite%rshort_g_diffuse(ipa)
     
     ! Absorption rate of short wave by the surface water
     do k=1,csite%nlev_sfcwater(ipa)
        csite%rshort_s_beam(k,ipa) = csite%rshort_s_beam(k,ipa) * beam_radiation
        csite%rshort_s_diffuse(k,ipa) = csite%rshort_s_diffuse(k,ipa) * rshort_diffuse
        csite%rshort_s(k,ipa) = csite%rshort_s_beam(k,ipa) + csite%rshort_s_diffuse(k,ipa)
     enddo
     
     csite%rlong_s_incid(ipa) = csite%rlong_s_incid(ipa) * rlong
     csite%rlong_g_incid(ipa) = csite%rlong_g_incid(ipa) * rlong
     csite%rlong_s(ipa) = csite%rlong_s_surf(ipa) + csite%rlong_s_incid(ipa)
     csite%rlong_g(ipa) = csite%rlong_g_surf(ipa) + csite%rlong_g_incid(ipa)

  enddo

  return
end subroutine scale_ed_radiation_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine angle_of_incid(aoi,cosz,solar_hour_aspect,slope,terrain_aspect)
  !!calculates angle of incidence based on local slope and aspect
  !!zenith is the angle of the sun with vertical
  !!solar_hour_aspect is the horizontal location of the sun defined with the same reference as terrain aspect
  implicit none
  real :: cosz,zenith,slope,solar_hour_aspect,terrain_aspect,aoi
  zenith = acos(cosz)
  aoi = cos(zenith)*cos(slope)+sin(zenith)*sin(slope)*cos(solar_hour_aspect-terrain_aspect)
  if(cosz .lt. 0.0) aoi = 0.0
  aoi = max(aoi,0.0)
  
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



