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
   integer                      :: ico
   integer                      :: k
   integer                      :: ksn
   integer                      :: ipft
   !---------------------------------------------------------------------------------------!


   !----- Find the total snow depth. ------------------------------------------------------!
   ksn = csite%nlev_sfcwater(ipa)
   csite%total_sfcw_depth(ipa) = 0.
   do k=1,ksn
      csite%total_sfcw_depth(ipa) = csite%total_sfcw_depth(ipa)                            &
                                  + csite%sfcwater_depth(k,ipa)
   end do
   !---------------------------------------------------------------------------------------!




   !----- Reset properties. ---------------------------------------------------------------!
   csite%veg_height(ipa)       = 0.0
   csite%lai(ipa)              = 0.0
   csite%wpa(ipa)              = 0.0
   csite%wai(ipa)              = 0.0
   weight_sum                  = 0.0
   csite%opencan_frac(ipa)     = 1.0
   csite%plant_ag_biomass(ipa) = 0.0
   !---------------------------------------------------------------------------------------!


   cpatch => csite%patch(ipa)

   !----- Loop over cohorts and integrate the patch-level properties. ---------------------!
   do ico = 1,cpatch%ncohorts

      ipft = cpatch%pft(ico)

      !----- Update the patch-level area indices. -----------------------------------------!
      csite%lai(ipa)  = csite%lai(ipa)  + cpatch%lai(ico)
      csite%wpa(ipa)  = csite%wpa(ipa)  + cpatch%wpa(ico)
      csite%wai(ipa)  = csite%wai(ipa)  + cpatch%wai(ico)
      !------------------------------------------------------------------------------------!



      !----- Compute the patch-level above-ground biomass
      csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)                            &
                                  + ed_biomass(cpatch%bdead(ico),cpatch%balive(ico)        &
                                              ,cpatch%bleaf(ico),cpatch%pft(ico)           &
                                              ,cpatch%hite(ico),cpatch%bstorage(ico)       &
                                              ,cpatch%bsapwood(ico))                       &
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



   !----- Find the fraction of vegetation buried in snow. ---------------------------------!
   csite%snowfac(ipa) = min(0.99, csite%total_sfcw_depth(ipa)/csite%veg_height(ipa))
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
!      This subroutine will take care of some diagnostic thermodynamic properties, namely  !
! the canopy air density and temperature.                                                  !
!------------------------------------------------------------------------------------------!
subroutine update_patch_thermo_props(csite,ipaa,ipaz,mzg,mzs,ntext_soil)
  
   use ed_state_vars, only : sitetype      ! ! structure
   use therm_lib    , only : idealdenssh   & ! function
                           , qwtk          & ! function
                           , qtk           ! ! function
   use consts_coms  , only : p00i          & ! intent(in)
                           , rocp          & ! intent(in)
                           , t00           & ! intent(in)
                           , wdns          ! ! intent(in)
   use soil_coms    , only : soil          ! ! intent(in)
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
   !---------------------------------------------------------------------------------------!


   do ipa=ipaa,ipaz

      csite%can_temp(ipa)     = csite%can_theta(ipa) * (p00i * csite%can_prss(ipa)) ** rocp
      csite%can_rhos(ipa)     = idealdenssh(csite%can_prss(ipa),csite%can_temp(ipa)        &
                                           ,csite%can_shv(ipa))

      !----- Update soil temperature and liquid water fraction. ---------------------------!
      do k = 1, mzg
         nsoil    = ntext_soil(k)
         soilhcap = soil(nsoil)%slcpd
         call qwtk(csite%soil_energy(k,ipa),csite%soil_water(k,ipa)*wdns,soilhcap          &
                  ,csite%soil_tempk(k,ipa),csite%soil_fracliq(k,ipa))
      end do

      !----- Update temporary surface water temperature and liquid water fraction. --------!
      ksn = csite%nlev_sfcwater(ipa)
      csite%total_sfcw_depth(ipa) = 0.
      do k = 1, ksn
         call qtk(csite%sfcwater_energy(k,ipa),csite%sfcwater_tempk(k,ipa)                 &
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
   end do

   return
end subroutine update_patch_thermo_props
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
   use allometry     , only : ed_biomass   ! ! function
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
   real                           :: ba
   !---------------------------------------------------------------------------------------!
   
   !----- Initialise the variables before looping. ----------------------------------------!
   cpoly%basal_area(:,:,isi) = 0.0
   cpoly%agb(:,:,isi)        = 0.0
   
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
!     This subroutine will update the derived properties at the polygon level.             !
!------------------------------------------------------------------------------------------!
subroutine update_polygon_derived_props(cgrid)

   use ed_state_vars , only : edtype      & ! structure
                            , polygontype ! ! structure
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)      , target  :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer :: cpoly
   integer                     :: ipy
   integer                     :: isi
   !---------------------------------------------------------------------------------------!

   do ipy=1,cgrid%npolygons

      cgrid%agb(:,:,ipy)        = 0.0
      cgrid%basal_area(:,:,ipy) = 0.0
      
      cpoly => cgrid%polygon(ipy)
      do isi = 1,cpoly%nsites
         cgrid%agb(:,:,ipy)        = cgrid%agb(:,:,ipy)                                    &
                                   + cpoly%area(isi) * cpoly%agb(:,:,isi)
         cgrid%basal_area(:,:,ipy) = cgrid%basal_area(:,:,ipy)                             &
                                   + cpoly%area(isi) * cpoly%basal_area(:,:,isi)
      end do
   end do

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
   use consts_coms   , only : cliqvlme     & ! intent(in)
                            , cicevlme     & ! intent(in)
                            , t3ple        & ! intent(in)
                            , tsupercool   ! ! intent(in)
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
                           if(csite%soil_tempk(k,ipa) > t3ple)then
                              csite%soil_energy(k,ipa) = soil(ntext)%slcpd                 &
                                                       * csite%soil_tempk(k,ipa)           &
                                                       + csite%soil_water(k,ipa)           &
                                                       * cliqvlme*(csite%soil_tempk(k,ipa) &
                                                                 - tsupercool)
                              csite%soil_fracliq(k,ipa) = 1.0
                           else
                              csite%soil_energy(k,ipa) = soil(ntext)%slcpd                 &
                                                       * csite%soil_tempk(k,ipa)           &
                                                       + csite%soil_water(k,ipa)           &
                                                       * cicevlme*csite%soil_tempk(k,ipa)
                              csite%soil_fracliq(k,ipa) = 0.0
                           end if
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
                                       ,csite%can_prss(ipa)                                &
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
!    This subroutine updates the 10-day running average of radiation, which is used for    !
! phenology.                                                                               !
!------------------------------------------------------------------------------------------!
subroutine update_rad_avg(cgrid)
   use ed_state_vars , only : edtype      & ! structure
                            , polygontype & ! structure
                            , sitetype    ! ! structure
   use ed_misc_coms     , only : radfrq      ! ! intent(in)
   use consts_coms   , only : day_sec     ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   type(edtype)     , target    :: cgrid
   !----- Local variables. ----------------------------------------------------------------!
   type(polygontype), pointer   :: cpoly
   type(sitetype)   , pointer   :: csite
   integer                      :: ipy
   integer                      :: isi
   real                         :: tfact
   !----- Local constants. ----------------------------------------------------------------!
   real             , parameter :: tendays_sec = 10.*day_sec
   !---------------------------------------------------------------------------------------!

   tfact = radfrq/tendays_sec


         
   polyloop: do ipy = 1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      siteloop: do isi = 1,cpoly%nsites
         cpoly%rad_avg(isi) = cpoly%rad_avg(isi) * (1.0 - tfact)                           &
                            + cpoly%met(isi)%rshort * tfact
      end do siteloop
   end do polyloop

   return
end subroutine update_rad_avg
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
