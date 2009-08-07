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
           call update_patch_derived_props(csite,cpoly%lsl(isi),cpoly%met(isi)%rhos,ipa)
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
subroutine update_patch_derived_props(csite,lsl,rhos,ipa)
  
   use ed_state_vars       , only : sitetype                   & ! structure
                                  , patchtype                  ! ! structure
   use allometry           , only : ed_biomass                 ! ! function
   use fusion_fission_coms , only : ff_ndbh                    ! ! intent(in)
   use fuse_fiss_utils     , only : patch_pft_size_profile     ! ! subroutine
   use canopy_air_coms     , only : veg_height_min             & ! intent(in)
                                  , minimum_canopy_depth       ! ! intent(in)
   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)  , target     :: csite
   integer         , intent(in) :: ipa
   integer         , intent(in) :: lsl
   real            , intent(in) :: rhos
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype) , pointer    :: cpatch
   real                         :: norm_fac, ba, old_height
   integer                      :: ico
   !----- External functions. -------------------------------------------------------------!
   real            , external   :: compute_water_storage
   real            , external   :: compute_energy_storage
   real            , external   :: compute_co2_storage
   !---------------------------------------------------------------------------------------!



   !----- Reset properties. ---------------------------------------------------------------!
   csite%veg_height(ipa)       = 0.0
   csite%lai(ipa)              = 0.0
   csite%wpa(ipa)              = 0.0
   csite%wai(ipa)              = 0.0
   norm_fac                    = 0.0
   csite%plant_ag_biomass(ipa) = 0.0

   cpatch => csite%patch(ipa)

   !----- Loop over cohorts. --------------------------------------------------------------!
   do ico = 1,cpatch%ncohorts
      
      !----- Compute contribution to height. ----------------------------------------------!
      ba                    = cpatch%nplant(ico) * cpatch%dbh(ico) * cpatch%dbh(ico)
      norm_fac              = norm_fac + ba
      csite%veg_height(ipa) = csite%veg_height(ipa) + cpatch%hite(ico) * ba
      
      !----- Update the patch-level area indices. -----------------------------------------!
      csite%lai(ipa)  = csite%lai(ipa)  + cpatch%lai(ico)
      csite%wpa(ipa)  = csite%wpa(ipa)  + cpatch%wpa(ico)
      csite%wai(ipa)  = csite%wai(ipa)  + cpatch%wai(ico)

      !----- Compute the patch-level above-ground biomass
      csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)                            &
                                  + ed_biomass(cpatch%bdead(ico),cpatch%balive(ico)        &
                                              ,cpatch%bleaf(ico),cpatch%pft(ico)           &
                                              ,cpatch%hite(ico),cpatch%bstorage(ico))      &
                                  * cpatch%nplant(ico)           
  
   end do
  
   !----- Find the PFT-dependent size distribution of this patch. -------------------------!
   call patch_pft_size_profile(csite,ipa,ff_ndbh)

  
   !----- Update vegetation height of this patch. -----------------------------------------!
   old_height = csite%veg_height(ipa)

   if (norm_fac > tiny(1.0)) then
      csite%veg_height(ipa) = max(veg_height_min,csite%veg_height(ipa) / norm_fac)
   else
      !---- In case we are dealing with an empty patch. -----------------------------------!
      csite%veg_height(ipa) = veg_height_min
   end if
   
   !----- Updating the canopy depth . -----------------------------------------------------!
   csite%can_depth(ipa) = max(csite%veg_height(ipa), minimum_canopy_depth)
   
   !----- Finding the patch roughness due to vegetation. ----------------------------------!
   csite%veg_rough(ipa) = 0.13 * csite%veg_height(ipa)
  
   !----- Computing the storage terms for CO2, energy, and water budgets. -----------------!
   csite%co2budget_initialstorage(ipa) = compute_co2_storage(csite,rhos,ipa)
   csite%wbudget_initialstorage(ipa)   = compute_water_storage(csite,lsl,rhos,ipa)
   csite%ebudget_initialstorage(ipa)   = compute_energy_storage(csite,lsl,rhos,ipa)

   !----- Updating the cohort count (may be redundant as well...) -------------------------!
   csite%cohort_count(ipa) = cpatch%ncohorts

   return
end subroutine update_patch_derived_props
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
   cpoly%agb_lu(:,isi)       = 0.0
   
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
            ba   = cpatch%nplant(ico) * cpatch%dbh(ico) * cpatch%dbh(ico)

            cpoly%basal_area(ipft,bdbh,isi) = cpoly%basal_area(ipft, bdbh,isi)             &
                                            + csite%area(ipa) * ba * pio4
            cpoly%agb(ipft,bdbh,isi) = cpoly%agb(ipft, bdbh,isi)                           &
                                     + ed_biomass(cpatch%bdead(ico),cpatch%balive(ico)     &
                                                 ,cpatch%bleaf(ico),cpatch%pft(ico)        &
                                                 ,cpatch%hite(ico),cpatch%bstorage(ico))   &
                                     * cpatch%nplant(ico) * 10.0 * csite%area(ipa)

            cpoly%agb_lu(ilu,isi) = cpoly%agb_lu(ilu,isi)                                  &
                                  + ed_biomass(cpatch%bdead(ico),cpatch%balive(ico)        &
                                              ,cpatch%bleaf(ico),cpatch%pft(ico)           &
                                              ,cpatch%hite(ico),cpatch%bstorage(ico))      &
                                  * cpatch%nplant(ico) * 10.0 * csite%area(ipa)
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
subroutine read_soil_moist_temp(cgrid)

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
   type(edtype)      , target    :: cgrid          ! Alias for current ED grid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer   :: cpoly          ! Alias for current polygon
   type(sitetype)    , pointer   :: csite          ! Alias for current site
   type(patchtype)   , pointer   :: cpatch         ! Alias for current patch
   integer                       :: ntext          !
   integer                       :: ilat           !
   integer                       :: ilon           !
   integer                       :: ilatf          !
   integer                       :: ilonf          !
   integer                       :: nls            !
   integer                       :: nlsw1          !
   integer                       :: k              !
   integer                       :: ipy            !
   integer                       :: isi            !
   integer                       :: ipa            !
   logical                       :: l1             !
   real                          :: glat           !
   real                          :: surface_temp   !
   real                          :: surface_fliq   !
   real                          :: glon           !
   real                          :: soil_tempaux   !
   real                          :: tmp1           !
   real                          :: tmp2           !
   real                          :: soilw1         !
   real                          :: soilw2         !
   !----- Local constants.  ---------------------------------------------------------------!
   logical           , parameter :: harvard_override = .false.
   integer           , parameter :: nlon = 144
   integer           , parameter :: nlat = 73
   real              , parameter :: dlon = 2.5
   real              , parameter :: dlat = 2.5
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
                           ntext = csite%ntext_soil(k,ipa)

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


                       !----- Initial condition is with no snow/pond. ----------------------!
                       csite%nlev_sfcwater(ipa) = 0
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
                        
                        nls = 1
                        call ed_grndvap(nls,csite%ntext_soil(nzg,ipa)                      &
                                       ,csite%soil_water(nzg,ipa)                          &
                                       ,csite%soil_energy(nzg,ipa)                         &
                                       ,csite%sfcwater_energy(nlsw1,ipa)                   &
                                       ,cpoly%met(isi)%rhos,csite%can_shv(ipa)             &
                                       ,csite%ground_shv(ipa),csite%surface_ssh(ipa)       &
                                       ,surface_temp,surface_fliq)

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
