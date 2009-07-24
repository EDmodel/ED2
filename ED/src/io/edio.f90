
 
subroutine ed_output(analysis_time,new_day,dail_analy_time,mont_analy_time,annual_time&
                    ,writing_dail,writing_mont,history_time,the_end)
  
  use ed_state_vars,only:filltab_alltypes,edgrid_g,filltables

  use grid_coms, only: ngrids,nzg  ! INTENT(IN)

  use ed_node_coms,only : mynum,nnodetot

  use ed_misc_coms, only: dtlsm, current_time, &
       idoutput,    &
       imoutput,    &
       iyoutput,    &
       isoutput,    &
       ifoutput,    &
       iprintpolys, &
       frqsum

  implicit none
    
  logical, intent(in)  :: the_end,analysis_time,dail_analy_time
  logical, intent(in)  :: writing_dail,writing_mont
  logical, intent(in) :: mont_analy_time,history_time,new_day,annual_time
  integer :: ifm


  ! If there is any IO, then we need to check if the pointer tables
  ! need to be rehashed, they will need to be rehashed if their has been 
  ! a change in the number of cohorts or patches, ie if a monthly event had
  ! just happened.

  if(analysis_time .or. history_time .or. dail_analy_time .or. mont_analy_time .or. annual_time ) then
     if(filltables) then
        
        ! Rehash the tables
        call filltab_alltypes
        ! Reset the rehash flag
        filltables=.false.

     endif
  endif



  if(analysis_time .or. history_time .or. (new_day .and. (writing_dail .or. writing_mont))) then
     do ifm=1,ngrids
        call normalize_averaged_vars(edgrid_g(ifm),frqsum,dtlsm)
     enddo

     !  Perform averaging and data preparation
     call spatial_averages
     
     if (writing_dail .or. writing_mont) then
        do ifm=1,ngrids
           call integrate_ed_daily_output_flux(edgrid_g(ifm))
        end do
     end if
  endif
  
  
  if (analysis_time) then
     
    
     !  Write out analysis fields - mostly polygon averages
     if (ifoutput.eq.3) then
        call h5_output('INST')
     endif
     
     ! If printpolys is on then print this info to
     ! the screen
     
     if (iprintpolys.eq.1) then
        do ifm=1,ngrids
           call print_fields(ifm,edgrid_g(ifm))
        enddo
     endif

  endif

  ! Daily analysis output and monthly integration
  if (new_day .and. (writing_dail .or. writing_mont)) then

     do ifm=1,ngrids
        call normalize_ed_daily_output_vars(edgrid_g(ifm))
        if (writing_mont) call integrate_ed_monthly_output_vars(edgrid_g(ifm))
     end do

     if (dail_analy_time) call h5_output('DAIL')

     do ifm=1,ngrids
       call zero_ed_daily_output_vars(edgrid_g(ifm))
     end do

  end if

  ! Monthly analysis output
  if (mont_analy_time) then

     do ifm=1,ngrids
        call normalize_ed_monthly_output_vars(edgrid_g(ifm))
     end do

     call h5_output('MONT')

     do ifm=1,ngrids
        call zero_ed_monthly_output_vars(edgrid_g(ifm))
     end do
  end if

  if (annual_time) then

     call h5_output('YEAR')

  endif

  ! History files should only be output at a frequency which
  ! divides by frqanl, thus the integrated fast-time variables
  ! are valid, but representative of the last frqanl period, not
  ! the last frqhist period.


  if(history_time) then

     if (isoutput /= 0) then
        call h5_output('HIST')
     end if
  endif



  return
end subroutine ed_output
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     The following subroutine performs several spatial averaging functions and temporal   !
! integrations.  Specifically, it area averages patch level quantities to the site, and    !
! site level quantities to the polygon.                                                    !
!------------------------------------------------------------------------------------------!
subroutine spatial_averages
   use ed_state_vars         , only : edtype            & ! structure
                                    , polygontype       & ! structure
                                    , sitetype          & ! structure
                                    , patchtype         & ! structure
                                    , edgrid_g          ! ! structure
   use grid_coms             , only : ngrids            & ! intent(in)
                                    , nzg               & ! intent(in)
                                    , nzs               ! ! intent(in)
   use consts_coms           , only : alvl              & ! intent(in)
                                    , wdns              ! ! intent(in)
   use ed_misc_coms             , only : frqsum            ! ! intent(in)
   use therm_lib             , only : qwtk              & ! subroutine
                                    , qtk               ! ! subroutine
   use soil_coms             , only : min_sfcwater_mass & ! intent(in)
                                    , soil              ! ! intent(in)
   use c34constants          , only : n_stoma_atts
   use ed_max_dims              , only : n_pft
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
   real                           :: lai_patch
   real                           :: laiarea_site
   real                           :: laiarea_poly
   real                           :: site_area_i
   real                           :: poly_area_i
   real                           :: site_avg_veg_hcap
   real                           :: poly_avg_veg_hcap
   real                           :: frqsumi
   real                           :: snowarea
   !---------------------------------------------------------------------------------------!

   !----- Time scale for output.  We will use the inverse more often. ---------------------!
   frqsumi = 1.0 / frqsum

   gridloop: do igr=1,ngrids
      cgrid => edgrid_g(igr)

      !------------------------------------------------------------------------------------!
      !    WARNING! cgrid variables should never be initialized outside the                !
      !             "do ipy=1,cgrid%npolygons" loop. npolygons is often 0 for some nodes   !
      !             on coupled runs (sudomains entirely over ocean), and initializing here !
      !             will cause either a crash or even worse, a memory leak.                !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy=1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         !----- Initialise some integrated variables --------------------------------------!
         area_sum           = 0.0
         laiarea_poly       = 0.0
         poly_avg_soil_hcap = 0.0
         poly_avg_veg_hcap  = 0.0

         !---------------------------------------------------------------------------------!
         !    Here is the safe place to initialize cgrid-level variables.                  !
         !---------------------------------------------------------------------------------!
         cgrid%avg_lai_ebalvars(:,:,ipy) = 0.0
         cgrid%avg_balive(ipy)      = 0.0
         cgrid%avg_bdead(ipy)       = 0.0
         cgrid%lai(ipy)             = 0.0
         cgrid%avg_gpp(ipy)         = 0.0
         cgrid%avg_leaf_resp(ipy)   = 0.0
         cgrid%avg_root_resp(ipy)   = 0.0
         cgrid%avg_plant_resp(ipy)  = 0.0
         cgrid%avg_htroph_resp(ipy) = 0.0
         cgrid%max_veg_temp(ipy)    = -huge(1.)
         cgrid%min_veg_temp(ipy)    =  huge(1.)
         cgrid%max_soil_temp(ipy)   = -huge(1.)
         cgrid%min_soil_temp(ipy)   =  huge(1.)

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
            cpoly%avg_vapor_vc(isi) = sum(csite%avg_vapor_vc * csite%area ) * site_area_i
            cpoly%avg_dew_cg(isi)   = sum(csite%avg_dew_cg   * csite%area ) * site_area_i
            cpoly%avg_vapor_gc(isi) = sum(csite%avg_vapor_gc * csite%area ) * site_area_i
            cpoly%avg_wshed_vg(isi) = sum(csite%avg_wshed_vg * csite%area ) * site_area_i
            cpoly%avg_vapor_ac(isi) = sum(csite%avg_vapor_ac * csite%area ) * site_area_i
            cpoly%avg_transp(isi)   = sum(csite%avg_transp   * csite%area ) * site_area_i
            cpoly%avg_evap(isi)     = sum(csite%avg_evap     * csite%area ) * site_area_i
            cpoly%avg_drainage(isi) = sum(csite%avg_drainage * csite%area ) * site_area_i
            cpoly%avg_runoff(isi)   = sum(csite%avg_runoff   * csite%area ) * site_area_i
            cpoly%aux(isi)          = sum(csite%aux          * csite%area ) * site_area_i
            cpoly%avg_sensible_vc(isi)   = sum(csite%avg_sensible_vc    * csite%area )     &
                                         * site_area_i
            cpoly%avg_sensible_2cas(isi) = sum(csite%avg_sensible_2cas  * csite%area )     &
                                         * site_area_i
            cpoly%avg_qwshed_vg(isi)     = sum(csite%avg_qwshed_vg      * csite%area )     &
                                         * site_area_i
            cpoly%avg_sensible_gc(isi)   = sum(csite%avg_sensible_gc    * csite%area )     &
                                         * site_area_i
            cpoly%avg_sensible_ac(isi)   = sum(csite%avg_sensible_ac    * csite%area )     &
                                         * site_area_i
            cpoly%avg_sensible_tot(isi)  = sum(csite%avg_sensible_tot   * csite%area )     &
                                         * site_area_i


            !----- Extra variables for NACP intercomparision (MCD) ------------------------!
            cpoly%avg_fsc(isi)    = sum(csite%fast_soil_C       * csite%area ) * site_area_i
            cpoly%avg_ssc(isi)    = sum(csite%slow_soil_C       * csite%area ) * site_area_i
            cpoly%avg_stsc(isi)   = sum(csite%structural_soil_C * csite%area ) * site_area_i


            !------------------------------------------------------------------------------!
            !------------------------------------------------------------------------------!
            !     Finding average soil properties.  The average soil temperature and       !
            ! liquid fraction is not directly computed, since the total mass and therefore !
            ! the total heat capacity (dry + water/ice) is different from each patch.      !
            !------------------------------------------------------------------------------!
            do k=cpoly%lsl(isi),nzg
               cpoly%avg_sensible_gg(k,isi) = sum(csite%avg_sensible_gg(k,:) * csite%area) &
                                            * site_area_i
               cpoly%avg_smoist_gg(k,isi)   = sum(csite%avg_smoist_gg(k,:)   * csite%area) &
                                            * site_area_i
               cpoly%avg_smoist_gc(k,isi)   = sum(csite%avg_smoist_gc(k,:)   * csite%area) &
                                            * site_area_i
               cpoly%aux_s(k,isi)           = sum(csite%aux_s(k,:)           * csite%area) &
                                            * site_area_i

               cpoly%avg_soil_energy(k,isi) = sum(csite%soil_energy(k,:)     * csite%area) &
                                            * site_area_i
               cpoly%avg_soil_water(k,isi)  = sum(csite%soil_water(k,:)      * csite%area) &
                                            * site_area_i

               !---------------------------------------------------------------------------!
               !    Finding the mean heat capacity. This shouldn't matter in the current   !
               ! version as all patches within the same site have the same soil texture.   !
               ! Since heat capacity is given in J/m3/K, it's safe to average it even if   !
               ! soil texture changed.                                                     !
               !---------------------------------------------------------------------------!
               site_avg_soil_hcap(k) = 0.
               do ipa=1,csite%npatches
                  site_avg_soil_hcap(k) = site_avg_soil_hcap(k)                            &
                                        + soil(csite%ntext_soil(k,ipa))%slcpd              &
                                        * csite%area(ipa) * site_area_i
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


            !------------------------------------------------------------------------------!
            !------------------------------------------------------------------------------!
            !    Temporary water/snow layer must be averaged differently, since each patch !
            ! may or may not have the layer, and the number of layers may vary from patch  !
            ! to patch.  Here we average the integrated energy, mass and depth, and        !
            ! compute the average temperature and liquid fraction from the average.        !
            !------------------------------------------------------------------------------!
            cpoly%avg_snowdepth(isi)   = 0.0
            cpoly%avg_snowenergy(isi)  = 0.0
            cpoly%avg_snowmass(isi)    = 0.0

            do ipa=1,csite%npatches
               ksn = csite%nlev_sfcwater(ipa)
               !----- Check whether mass is present.  If so, add this patch. --------------!
               if (ksn > 0) then
                  do k=1,ksn
                     cpoly%avg_snowdepth(isi)   = cpoly%avg_snowdepth(isi)                 &
                                                + csite%sfcwater_depth(k,ipa)              &
                                                * csite%area(ipa) * site_area_i
                     cpoly%avg_snowmass(isi)    = cpoly%avg_snowmass(isi)                  &
                                                + csite%sfcwater_mass(k,ipa)               &
                                                * csite%area(ipa) * site_area_i
                     !---------------------------------------------------------------------!
                     !     Internal energy.  For averaging we use the extensive one (J/m2) !
                     ! we will switch back to J/kg after the mean mass is found.           !
                     !---------------------------------------------------------------------!
                     cpoly%avg_snowenergy(isi)  = cpoly%avg_snowenergy(isi)                &
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
            if (cpoly%avg_snowmass(isi) > min_sfcwater_mass) then
               cpoly%avg_snowenergy(isi) = cpoly%avg_snowenergy(isi)                       &
                                         / cpoly%avg_snowmass(isi)
               call qtk(cpoly%avg_snowenergy(isi),cpoly%avg_snowtempk(isi)                 &
                       ,cpoly%avg_snowfracliq(isi))
            else
               cpoly%avg_snowmass(isi)    = 0.
               cpoly%avg_snowdepth(isi)   = 0.
               cpoly%avg_snowenergy(isi)  = 0.
               cpoly%avg_snowtempk(isi)   = cpoly%avg_soil_temp(nzg,isi)
               cpoly%avg_snowfracliq(isi) = cpoly%avg_soil_fracliq(nzg,isi)
            end if
            !------------------------------------------------------------------------------!


            !----- Average over patches. --------------------------------------------------!
            longpatchloop: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               !---------------------------------------------------------------------------!
               !     Adding cohort "extensive" variables. Those that are not must be       !
               ! scaled by nplant. Just make sure that we have at least one cohort.        !
               !---------------------------------------------------------------------------!
               if (cpatch%ncohorts > 0) then
                  lai_patch = sum(cpatch%lai, cpatch%solvable)
                  csite%avg_veg_energy(ipa) = sum(cpatch%veg_energy)
                  csite%avg_veg_water(ipa)  = sum(cpatch%veg_water)
                  csite%hcapveg(ipa)        = sum(cpatch%hcapveg)
                  call qwtk(csite%avg_veg_energy(ipa),csite%avg_veg_water(ipa)             &
                           ,csite%hcapveg(ipa),csite%avg_veg_temp(ipa)                     &
                           ,csite%avg_veg_fliq(ipa))

                  cgrid%avg_gpp(ipy)        = cgrid%avg_gpp(ipy)                           &
                                            + sum(cpatch%mean_gpp)                         &
                                            * csite%area(ipa)*cpoly%area(isi)              &
                                            * site_area_i * poly_area_i

                  cgrid%avg_leaf_resp(ipy)  = cgrid%avg_leaf_resp(ipy)                     &
                                            + sum(cpatch%mean_leaf_resp)                   &
                                            * csite%area(ipa)*cpoly%area(isi)              &
                                            * site_area_i * poly_area_i

                  cgrid%avg_root_resp(ipy)  = cgrid%avg_root_resp(ipy)                     &
                                            + sum(cpatch%mean_root_resp)                   &
                                            * csite%area(ipa)*cpoly%area(isi)              &
                                            * site_area_i * poly_area_i

                  cgrid%avg_balive(ipy)     = cgrid%avg_balive(ipy)                        &
                                            + sum(cpatch%balive*cpatch%nplant)             &
                                            * csite%area(ipa)*cpoly%area(isi)              &
                                            * site_area_i * poly_area_i

                  cgrid%avg_bdead(ipy)      = cgrid%avg_bdead(ipy)                         &
                                            + sum(cpatch%bdead*cpatch%nplant)              &
                                            * csite%area(ipa)*cpoly%area(isi)              &
                                            * site_area_i * poly_area_i

                  !----- Check the extremes and update if necessary. ----------------------!
                  if (maxval(cpatch%veg_temp) > cgrid%max_veg_temp(ipy)) then
                     cgrid%max_veg_temp(ipy) = maxval(cpatch%veg_temp)
                  end if
                  if (minval(cpatch%veg_temp) < cgrid%min_veg_temp(ipy)) then
                     cgrid%min_veg_temp(ipy) = minval(cpatch%veg_temp)
                  end if
               else
                  lai_patch                 = 0.
                  csite%avg_veg_energy(ipa) = 0.
                  csite%avg_veg_water(ipa)  = 0.
                  csite%hcapveg(ipa)        = 0.
                  csite%avg_veg_temp(ipa)   = csite%can_temp(ipa)
                  csite%avg_veg_fliq(ipa)   = 0.
               end if

               if (lai_patch > 0.) then
                  csite%laiarea(ipa) = csite%area(ipa)
               else
                  csite%laiarea(ipa) = 0.
               end if

               !---------------------------------------------------------------------------!
               !    Updating some other flux variables that need to be scaled by frqsum.   !
               !---------------------------------------------------------------------------!
               cgrid%avg_plant_resp(ipy)  = cgrid%avg_plant_resp(ipy)                      &
                                          + csite%co2budget_plresp(ipa)                    &
                                          * csite%area(ipa)*cpoly%area(isi)                &
                                          * site_area_i * poly_area_i * frqsumi

               cgrid%avg_htroph_resp(ipy) = cgrid%avg_htroph_resp(ipy)                     &
                                          + csite%co2budget_rh(ipa)                        &
                                          * csite%area(ipa)*cpoly%area(isi)                &
                                          * site_area_i * poly_area_i * frqsumi

               !----- Not sure what these variables do. -----------------------------------!
               lai_index = min(3,max(1, floor(csite%lai(ipa)/2.0) + 1)  )
               area_sum(lai_index) = area_sum(lai_index) + csite%area(ipa)

               !----- Net radiation. ------------------------------------------------------!
               cgrid%avg_lai_ebalvars(lai_index,1,ipy) =                                   &
                      cgrid%avg_lai_ebalvars(lai_index,1,ipy)                              &
                    + csite%avg_netrad(ipa)                                                &
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

            !----- Site average of canopy thermodynamic state -----------------------------!
            cpoly%avg_can_temp(isi) = sum(csite%can_temp  * csite%area) * site_area_i
            cpoly%avg_can_shv(isi)  = sum(csite%can_shv   * csite%area) * site_area_i
            cpoly%avg_can_co2(isi)  = sum(csite%can_co2   * csite%area) * site_area_i

            !------------------------------------------------------------------------------!
            !   Site average of leaf properties.  Again, we average "extensive" properties !
            ! and find the average temperature based on the average leaf internal energy.  !
            ! mass, and heat capacity.                                                     !
            !------------------------------------------------------------------------------!
            cpoly%avg_veg_energy(isi) = sum(csite%avg_veg_energy * csite%area)*site_area_i
            cpoly%avg_veg_water(isi)  = sum(csite%avg_veg_water  * csite%area)*site_area_i
            site_avg_veg_hcap         = sum(csite%hcapveg        * csite%area)*site_area_i
            !----- Also integrate the polygon-level mean leaf heat capacity. --------------!
            poly_avg_veg_hcap         = poly_avg_veg_hcap                                  &
                                      + site_avg_veg_hcap * cpoly%area(isi) * poly_area_i
            !------------------------------------------------------------------------------!
            !    Unless this is a bare site or it is absolute leafless, there should be    !
            ! some heat capacity, then compute the average leaf temperature. Otherwise,    !
            ! assign mean canopy temperature.                                              !
            !------------------------------------------------------------------------------!
            if (laiarea_site > 0.) then
               call qwtk(cpoly%avg_veg_energy(isi),cpoly%avg_veg_water(isi)                &
                        ,site_avg_veg_hcap,cpoly%avg_veg_temp(isi),cpoly%avg_veg_fliq(isi))
            else
               cpoly%avg_veg_temp(isi) = cpoly%avg_can_temp(isi)
               cpoly%avg_veg_fliq(isi) = 0.
            end if

         end do siteloop
     

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
         cgrid%avg_vapor_vc(ipy)     = sum(cpoly%avg_vapor_vc    * cpoly%area)*poly_area_i
         cgrid%avg_dew_cg(ipy)       = sum(cpoly%avg_dew_cg      * cpoly%area)*poly_area_i
         cgrid%avg_vapor_gc(ipy)     = sum(cpoly%avg_vapor_gc    * cpoly%area)*poly_area_i
         cgrid%avg_wshed_vg(ipy)     = sum(cpoly%avg_wshed_vg    * cpoly%area)*poly_area_i
         cgrid%avg_fsc(ipy)          = sum(cpoly%avg_fsc         * cpoly%area)*poly_area_i
         cgrid%avg_stsc(ipy)         = sum(cpoly%avg_stsc        * cpoly%area)*poly_area_i
         cgrid%avg_ssc(ipy)          = sum(cpoly%avg_ssc         * cpoly%area)*poly_area_i
         cgrid%avg_runoff_heat(ipy)  = sum(cpoly%avg_runoff_heat * cpoly%area)*poly_area_i
         cgrid%avg_runoff(ipy)       = sum(cpoly%avg_runoff      * cpoly%area)*poly_area_i
         cgrid%avg_drainage(ipy)     = sum(cpoly%avg_drainage    * cpoly%area)*poly_area_i
         cgrid%avg_vapor_ac(ipy)     = sum(cpoly%avg_vapor_ac     *cpoly%area)*poly_area_i
         cgrid%avg_transp(ipy)       = sum(cpoly%avg_transp       *cpoly%area)*poly_area_i
         cgrid%avg_evap(ipy)         = sum(cpoly%avg_evap         *cpoly%area)*poly_area_i
         cgrid%aux(ipy)              = sum(cpoly%aux              *cpoly%area)*poly_area_i
         cgrid%avg_sensible_vc(ipy)  = sum(cpoly%avg_sensible_vc  *cpoly%area)*poly_area_i
         cgrid%avg_sensible_2cas(ipy)= sum(cpoly%avg_sensible_2cas*cpoly%area)*poly_area_i
         cgrid%avg_qwshed_vg(ipy)    = sum(cpoly%avg_qwshed_vg    *cpoly%area)*poly_area_i
         cgrid%avg_sensible_gc(ipy)  = sum(cpoly%avg_sensible_gc  *cpoly%area)*poly_area_i
         cgrid%avg_sensible_ac(ipy)  = sum(cpoly%avg_sensible_ac  *cpoly%area)*poly_area_i
         cgrid%avg_sensible_tot(ipy) = sum(cpoly%avg_sensible_tot *cpoly%area)*poly_area_i
         cgrid%avg_can_temp(ipy)     = sum(cpoly%avg_can_temp     *cpoly%area)*poly_area_i
         cgrid%avg_can_shv(ipy)      = sum(cpoly%avg_can_shv      *cpoly%area)*poly_area_i
         cgrid%avg_can_co2(ipy)      = sum(cpoly%avg_can_co2      *cpoly%area)*poly_area_i

         !---------------------------------------------------------------------------------!
         !    Similar to the site level, average mass, heat capacity and energy then find  !
         ! the average temperature and liquid water fraction.                              !
         !---------------------------------------------------------------------------------!
         do k=cgrid%lsl(ipy),nzg
            cgrid%avg_sensible_gg(k,ipy) = sum(cpoly%avg_sensible_gg(k,:)*cpoly%area)      &
                                         * poly_area_i
            cgrid%avg_smoist_gg(k,ipy)   = sum(cpoly%avg_smoist_gg(k,:)  *cpoly%area)      &
                                         * poly_area_i
            cgrid%avg_smoist_gc(k,ipy)   = sum(cpoly%avg_smoist_gc(k,:)  *cpoly%area)      &
                                         * poly_area_i
            cgrid%aux_s(k,ipy)           = sum(cpoly%aux_s(k,:)          *cpoly%area)      &
                                         * poly_area_i
            cgrid%avg_soil_energy(k,ipy) = sum(cpoly%avg_soil_energy(k,:)*cpoly%area)      &
                                         * poly_area_i
            cgrid%avg_soil_water(k,ipy)  = sum(cpoly%avg_soil_water(k,:) *cpoly%area)      &
                                         * poly_area_i

            !------------------------------------------------------------------------------!
            !     Finding the average temperature and liquid fraction.  The polygon-level  !
            ! mean heat capacity was already found during the site loop.                   !
            !------------------------------------------------------------------------------!
            call qwtk(cgrid%avg_soil_energy(k,ipy),cgrid%avg_soil_water(k,ipy)*wdns        &
                     ,poly_avg_soil_hcap(k),cgrid%avg_soil_temp(k,ipy)                     &
                     ,cgrid%avg_soil_fracliq(k,ipy))
         end do


         !---------------------------------------------------------------------------------!
         !    Also using the same idea as the site-level: average energy, mass, and depth, !
         ! but find temperature and liquid fraction with the averaged values. Again, use   !
         ! the energy in J/m2 instead of J/kg to do the polygon averaging.                 !
         !---------------------------------------------------------------------------------!
         cgrid%avg_snowmass(ipy)     = sum(cpoly%avg_snowmass  * cpoly%area) * poly_area_i
         cgrid%avg_snowdepth(ipy)    = sum(cpoly%avg_snowdepth * cpoly%area) * poly_area_i
         cgrid%avg_snowenergy(ipy)   = sum( cpoly%avg_snowenergy                           &
                                          * cpoly%avg_snowmass * cpoly%area) * poly_area_i

         !----- Scale energy and find temp and fracliq if there is enogh mass -------------!
         if (cgrid%avg_snowmass(ipy) > min_sfcwater_mass) then
            cgrid%avg_snowenergy(ipy) = cgrid%avg_snowenergy(ipy) / cgrid%avg_snowmass(ipy)
            call qtk(cgrid%avg_snowenergy(ipy),cgrid%avg_snowtempk(ipy)                    &
                    ,cgrid%avg_snowfracliq(ipy))
         else
            cgrid%avg_snowmass(ipy)    = 0.
            cgrid%avg_snowdepth(ipy)   = 0.
            cgrid%avg_snowenergy(ipy)  = 0.
            cgrid%avg_snowtempk(ipy)   = cgrid%avg_soil_temp(nzg,ipy)
            cgrid%avg_snowfracliq(ipy) = cgrid%avg_soil_fracliq(nzg,ipy)
         end if

         !---------------------------------------------------------------------------------!
         !    Similar to site level, compute mean leaf internal energy and water mass.     !
         ! Also find the mean heat capacity.  If there is enough LAI, then find the mean   !
         ! temperature, otherwise just assume to be the canopy temperature.                !
         !---------------------------------------------------------------------------------!
         cgrid%avg_veg_energy(ipy) = sum(cpoly%avg_veg_energy * cpoly%area) * poly_area_i
         cgrid%avg_veg_water(ipy)  = sum(cpoly%avg_veg_water  * cpoly%area) * poly_area_i

         if (laiarea_poly > 0.) then
            call qwtk(cgrid%avg_veg_energy(ipy),cgrid%avg_veg_water(ipy),poly_avg_veg_hcap &
                     ,cgrid%avg_veg_temp(ipy),cgrid%avg_veg_fliq(ipy))
         else
            cgrid%avg_veg_temp(ipy) = cgrid%avg_can_temp(ipy)
            cgrid%avg_veg_fliq(ipy) = 0.
         end if

      end do polyloop
   end do gridloop

   return
end subroutine spatial_averages

! ==============================

subroutine get2d(m1,m2,temp2,in_ptr)
  
  implicit none
  integer :: m1,m2,i,j
  real,dimension(m1,m2) :: temp2
  real,dimension(m1,m2) :: in_ptr
  do i = 1,m1
     do j = 1,m2
        temp2(i,j) = in_ptr(i,j)
     enddo
  enddo
  return
end subroutine get2d

! =============================

subroutine get3d(m1,m2,m3,temp3,in_ptr)
  
  implicit none
  integer :: m1,m2,m3,i,j,k
  real,dimension(m1,m2,m3) :: temp3
  real,dimension(m1,m2,m3) :: in_ptr
  do i = 1,m1
     do j = 1,m2
        do k = 1,m3
           temp3(i,j,k) = in_ptr(i,j,k)
        enddo
     enddo
  enddo
  return
end subroutine get3d
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine print_fields(ifm,cgrid)
  
  !------------------------------------------------------
  ! PRINT OUT FIELDS OF INTEREST
  !
  ! This subroutine prints patch, cohort, polygon or site level
  ! data, upscales that data to the site level and stores
  ! the data in its spatial array coordinate.
  ! The data is then printed to the screen, based on a
  ! specified window of data.
  ! Note, this may be printing windows on various nodes,
  ! or this may be called from a master process.  Be
  ! conscious of this; as it will dictate what part of the
  ! domain you are accessing variables from, and whether
  ! or not the variable of interest is stored in memory
  ! at that time.  For instance, many variables are stored
  ! only on the slave nodes, and need not be passed back
  ! to the master.  Likewise, many slave node data will
  ! accumulate after each lsm call, until they are passed back
  ! to the master, where they are normalized. These variables
  ! will be immediately zeroed on the slaves after being
  ! sent to the master.
  ! Dont forget to adjust the number precision on the 
  ! format string at the end.  The X.Xf
  !--------------------------------------------------------
  
  use ed_node_coms,only: mynum,nnodetot,sendnum,recvnum,master_num,machs
  use ed_state_vars,only: edtype,polygontype
  use ed_misc_coms, only: &
            printvars,  &
            ipmax,      &
            ipmin,      &
            pfmtstr,    &
            iprintpolys
  
  use ed_var_tables,only:vt_info,num_var


  implicit none

  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE) :: status
  integer  :: ifm,nv,np,i,ip
  integer  :: ping
  integer  :: npolys
  integer  :: g_idmin,g_idmax,l_idmin,l_idmax,ierr
  integer  :: node_idmin,node_idmax
  integer  :: mast_idmin,mast_idmax
  integer  :: g_id,g_ln,nm
  integer  :: ncols,row,maxrows,col
  integer,parameter :: maxcols = 10

  real,pointer,dimension(:) :: pvar_l
  real,pointer,dimension(:) :: pvar_g
  
  character(len=30)     :: fmtstr
  character(len=32)     :: pvar_name
  
  ! Linked structures
  type(edtype),target     :: cgrid
  
  logical :: pvartrue
  logical :: ptr_recv
  logical :: ptr_send

  ping = 8675309

  
  npolys = ipmax - ipmin + 1
  
  ! Adjust the format string according to the chosen
  ! Variables
  ! ------------------------------------------------

  ! CHeck the window size
  ! ---------------------
  
  if (ipmax.gt.cgrid%npolygons_global) then
     print*,"========================================="
     print*,"You have specified a print index"
     print*,"greater than the total number of polygons"
     print*,"You must reduce this number. Stopping"
     stop
  end if


  ! Allocate the print and scratch vector
  allocate(pvar_l(npolys))

  if (mynum .eq. nnodetot .or. nnodetot .eq. 1) allocate(pvar_g(npolys))

  
  ! Loop through the printvar entries from the namelist

  ip = 0
  
  count_pvars: do

     ip = ip+1
     pvar_name = printvars(ip)
     if (len_trim(pvar_name).eq.32) then
        
        exit count_pvars
     endif

!     call MPI_Barrier(MPI_COMM_WORLD,ierr)

     pvartrue = .false.
     do nv = 1,num_var(ifm)

        if(trim(vt_info(nv,ifm)%name) .eq. trim(pvar_name)) then
           pvartrue = .true.
           
           ! If this is root, then collect the sends, keep reading to find out what it is
           ! receiving
           if (nnodetot.gt.1) then
              
              if (mynum .eq. nnodetot) then
                 
                 pvar_g = -99.9
                 ! Loop through the variable table to match the print variable
                 print*,""
                 print*," ============ ",trim(pvar_name)," =================="
                 print*,""

                 do nm = 1,nnodetot-1

                    call MPI_Recv(ptr_recv,1,MPI_LOGICAL,machs(nm),120,MPI_COMM_WORLD,status,ierr)

                    if (ptr_recv) then
                       call MPI_Recv(mast_idmin,1,MPI_INTEGER,machs(nm),121,MPI_COMM_WORLD,status,ierr)
                       call MPI_Recv(mast_idmax,1,MPI_INTEGER,machs(nm),122,MPI_COMM_WORLD,status,ierr)
                       call MPI_Recv(pvar_g(mast_idmin:mast_idmax),mast_idmax-mast_idmin+1,MPI_REAL,&
                            machs(nm),123,MPI_COMM_WORLD,status,ierr)
                    end if
                 enddo
              endif
           else
              
              ! Loop through the variable table to match the print variable
              print*,""
              print*," ============ ",trim(pvar_name)," =================="
              print*,""
              pvar_g = -99.9
           endif
        
              
           ! The namelist print entry has been matched with the var_table
           ! entry.  Now lets cycle through our machines and determine
           ! if those machines hold data that should be printed.  If so, 
           ! then send that data to the master (machine 0).

           ! This is scratch space that all machines will use

           pvar_l =-99.9
           
           ! Set the blocking recieve to allow ordering, start with machine 1
           if (mynum /= 1) call MPI_Recv(ping,1,MPI_INTEGER,recvnum,93,MPI_COMM_WORLD,status,ierr)
           
           ! Cycle through this node's pointers for the current variable.  If the index
           ! falls within the printable range. Save that pointer to a local array.
           ! Once all the pointers have been cycled, send the local array to the 
           ! master to populate a global array and print
           
           ptr_send = .false.
           node_idmin = -1
           node_idmax = -1

           do np = 1,vt_info(nv,ifm)%nptrs
              g_id = vt_info(nv,ifm)%vt_vector(np)%globid+1
              g_ln = vt_info(nv,ifm)%vt_vector(np)%varlen
              
              ! Determine if any of this segment falls within the 
              ! range desired for output, globid+1 is the global index

              if (g_id .le. ipmax .and. g_id+g_ln-1 .ge. ipmin ) then

                 ! OK, this segment is good, set the send flag
                 ptr_send = .true.

                 ! These are the indices of the data in the current segment to use
                 ! and the indices in the global array they will be sent to
                 if (g_id >= ipmin) then
                    l_idmin = 1
                    g_idmin = g_id - ipmin + 1
                 else
                    l_idmin = ipmin - g_id + 1
                    g_idmin = 1
                 endif

                 if (g_id+g_ln-1 < ipmax) then
                    l_idmax = g_ln
                    g_idmax = g_id + g_ln - ipmin
                 else
                    l_idmax = ipmax - g_id + 1
                    g_idmax = ipmax - ipmin + 1
                 endif
                 
                 ! These should be the same size, if not...
                 if (l_idmax - l_idmin .ne. g_idmax - g_idmin ) then
                    print*,"NOT THE SAME LENGTHS"
                    print*,l_idmax - l_idmin,g_idmax - g_idmin
                    print*,l_idmax,l_idmin,g_idmax,g_idmin
                    stop
                 endif

                 ! Shift the global dataset so that it is applied to the
                 ! first index

                 call fillvar_l(pvar_l,vt_info(nv,ifm)%vt_vector(np)%var_rp,npolys,g_ln,g_idmin,l_idmin,l_idmax)

                 ! Determine the minimum and maximum indices that will be sent
                 if (g_idmin < node_idmin .or. node_idmin.eq.-1) node_idmin = g_idmin
                 if (g_idmax > node_idmax .or. node_idmax.eq.-1) node_idmax = g_idmax

              end if
              
           enddo
           
           if (nnodetot.gt.1) then
              
              ! The local array for this machine has been created. Send it off to the master

              if (mynum /= nnodetot) then

                 call MPI_Send(ptr_send,1,MPI_LOGICAL,machs(nnodetot),120,MPI_COMM_WORLD,ierr)
                 if (ptr_send) then
                    call MPI_Send(node_idmin,1,MPI_INTEGER,machs(nnodetot),121,MPI_COMM_WORLD,ierr)
                    call MPI_Send(node_idmax,1,MPI_INTEGER,machs(nnodetot),122,MPI_COMM_WORLD,ierr)
                    call MPI_Send(pvar_l(node_idmin:node_idmax),node_idmax-node_idmin+1, &
                         MPI_REAL,machs(nnodetot),123,MPI_COMM_WORLD,ierr)
                 end if
                 
              
                 ! When this node is finished, send the blocking MPI_Send to the next machine

                 call MPI_Send(ping,1,MPI_INTEGER,sendnum,93,MPI_COMM_WORLD,ierr)
                 
                 ! If this is root, then just copy the array to the global
              else
                            
                 if (ptr_send) then

                    pvar_g(node_idmin:node_idmax) = pvar_l(node_idmin:node_idmax)
           
                 end if
                 
              endif


              
           else
              
              pvar_g = pvar_l
              
           endif
           
           
           ! The data over the desired range of indices have been collected
           ! if this is the only machine or the root machine, then print it 
           ! to standard output
           
           if (mynum .eq. nnodetot .or. nnodetot .eq. 1) then

              ! Print out a maximum of 10 variables per row...

              maxrows = ceiling(real(npolys)/real(maxcols))

              do row = 1,maxrows
                 
                 ncols = min( maxcols,npolys-((row-1)*maxcols)   )
                 col   = ( (row-1)*maxcols)+1
                 
                 write(fmtstr,'(i3)')ncols
                 fmtstr = '('// trim(fmtstr)  // '(2x,' // trim(pfmtstr(ip)) // '))'
                 
                 print(trim(fmtstr)),(pvar_g(i),i=col,col+ncols-1)
              enddo
              print*,""
              print*,""
              
           endif
           
        endif
        
     enddo

     ! Check to see if we matched the variable
     if(.not.pvartrue) then
        print*,"The diagnostic variable named:",trim(pvar_name)
        print*,"does not match any of the var_table variables"
        print*,"Check you namelist entries, and the variable "
        print*,"registry and/or remove this"
        print*,"diagnostic variable."
        stop
     endif

  enddo count_pvars

  ! Dont proceed until everything is written out
  ! if this is not done, then there will be other writing
  ! contaminating the output, and thats icky

!  call MPI_Barrier(MPI_COMM_WORLD,ierr)

  
  return
end subroutine print_fields

! =======================================================

subroutine fillvar_l(pvar_l,vt_ptr,npts_out,npts_in,out1,in1,in2)
  
  implicit none
  real,dimension(npts_out)  :: pvar_l
  real,dimension(npts_in)   :: vt_ptr
  integer,intent(in)      :: npts_in,npts_out
  integer,intent(in)      :: out1,in1,in2
  integer :: i,j
  
  j = out1
  do i = in1,in2     
     pvar_l(j) = vt_ptr(i)
     j = j + 1
  enddo
  return
end subroutine fillvar_l

