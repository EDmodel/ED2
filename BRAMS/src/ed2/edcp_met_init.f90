!==========================================================================================!
!==========================================================================================!
!     This subroutine will initialise the upwelling radiation (surface emission) and       !
! albedo.                                                                                  !
!------------------------------------------------------------------------------------------!
subroutine ed_init_radiation()

   use mem_radiate  , only : radiate_g ! ! intent(inout)
   use mem_leaf     , only : leaf_g    ! ! intent(inout)
   use mem_grid     , only : ngrids    ! ! intent(in)
   use rconstants   , only : stefan    ! ! intent(in)
   use ed_state_vars, only : edgrid_g  & ! intent(inout)
                           , edtype    ! ! structure

   implicit none
   !----- Local variables. ----------------------------------------------------------------!
  
   type(edtype)     , pointer :: cgrid
   integer                    :: igr
   integer                    :: ipy
   integer                    :: ix
   integer                    :: iy
   !---------------------------------------------------------------------------------------!

   gridloop: do igr=1,ngrids
      cgrid => edgrid_g(igr)

      !----- First, do a catch all using sst and crude albedo. ----------------------------!
      radiate_g(igr)%rlongup  = 0.98 * stefan * leaf_g(igr)%seatp**4
      radiate_g(igr)%albedt   = 0.3

      !------------------------------------------------------------------------------------!
      !      Then use the air temperature to initialize a brightness temperature.  This is !
      ! not ideal, but the air temperature just initialized the land-surface variables.    !
      ! These radiation parameters are just to prevent the model from crashing.            !
      !------------------------------------------------------------------------------------!
      polyloop: do ipy = 1,cgrid%npolygons
         ix = cgrid%ilon(ipy)
         iy = cgrid%ilat(ipy)

         radiate_g(igr)%rlongup(ix,iy) = 0.97 * stefan * cgrid%met(ipy)%atm_tmp**4
         radiate_g(igr)%albedt (ix,iy) = 0.3
      end do polyloop
   end do gridloop

   return
end subroutine ed_init_radiation
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This sub-routine copies the soil moisture information from LEAF-3 to ED.              !
!------------------------------------------------------------------------------------------!
subroutine leaf2ed_soil_moist_energy(cgrid,ifm)
   use ed_state_vars, only : edtype       & ! structure
                           , polygontype  & ! structure
                           , sitetype     & ! structure
                           , patchtype    ! ! structure
   use grid_coms    , only : nzg          & ! intent(in)
                           , nzs          ! ! intent(in)
   use ed_therm_lib , only : ed_grndvap   ! ! subroutine
   use therm_lib    , only : uextcm2tl    & ! subroutine
                           , cmtl2uext    ! ! function
   use rconstants   , only : wdns         ! ! intent(in)
   use mem_leaf     , only : leaf_g       ! ! structure
   use leaf_coms    , only : slcpd        ! ! intent(in)
   use soil_coms    , only : soil         ! ! intent(in)
   
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)         , target      :: cgrid           ! Alias for current ED grid
   integer              , intent(in)  :: ifm             ! Current grid
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype)    , pointer     :: cpoly           ! Alias for current polygon
   type(sitetype)       , pointer     :: csite           ! Alias for current site
   type(patchtype)      , pointer     :: cpatch          ! Alias for current patch
   integer                            :: ntext           ! Alias for ED-2 soil text. class
   integer, dimension(:), allocatable :: lsoil_text      ! LEAF-3 soil texture class
   integer                            :: ipy             ! Counters for polygons
   integer                            :: isi             ! Counters for sites
   integer                            :: ipa             ! Counters for patches
   integer                            :: ico             ! Counters for cohorts
   integer                            :: k               ! Counter for soil layers
   integer                            :: ix              ! Counter for longitude
   integer                            :: iy              ! Counter for latitude
   integer                            :: ilp             ! Counter for LEAF-3 patch
   integer                            :: ksn             ! Alias for # of pond/snow layers
   integer                            :: ksnw1           ! Alias for # of pond/snow layers
   real   , dimension(:), allocatable :: lsoil_temp      ! LEAF-3 soil temperature
   real   , dimension(:), allocatable :: lsoil_fliq      ! LEAF-3 soil liquid fraction
   real                               :: fice            ! soil ice fraction
   !---------------------------------------------------------------------------------------!


   !----- Allocate the scratch arrays. ----------------------------------------------------!
   allocate (lsoil_text(nzg))
   allocate (lsoil_temp(nzg))
   allocate (lsoil_fliq(nzg))

   !----- Loop over land points -----------------------------------------------------------!
   polyloop: do ipy=1,cgrid%npolygons
      ix = cgrid%ilon(ipy)
      iy = cgrid%ilat(ipy)

      cpoly => cgrid%polygon(ipy)

      !----- Loop over sites --------------------------------------------------------------!
      siteloop: do isi=1,cpoly%nsites

         !---------------------------------------------------------------------------------!
         !     ILP is the LEAF-3 "patch", which is offset by one because the first one is  !
         ! always water.                                                                   !
         !---------------------------------------------------------------------------------!
         ilp = isi + 1
         !---------------------------------------------------------------------------------!


         !------ Alias for current site. --------------------------------------------------!
         csite => cpoly%site(isi)
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !    Determine initial soil temperature and liquid fraction.  The reason we find  !
         ! this for LEAF-3 instead of simply copying the soil energy and water to ED-2 is  !
         ! that depending on the way the user set up his or her RAMSIN, soil types may not !
         ! match and this could put ED-2.1 in an inconsistent initial state.               !
         !---------------------------------------------------------------------------------!
         do k=1,nzg
            lsoil_text(k) =nint(leaf_g(ifm)%soil_text(k,ix,iy,ilp))
            call uextcm2tl(leaf_g(ifm)%soil_energy(k,ix,iy,ilp)                            &
                          ,leaf_g(ifm)%soil_water(k,ix,iy,ilp)*wdns                        &
                          ,slcpd(lsoil_text(k)),lsoil_temp(k),lsoil_fliq(k))
         end do
         !---------------------------------------------------------------------------------!


         !----- Loop over patches ---------------------------------------------------------!
         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)

            do k=1,nzg
            
               ntext = cpoly%ntext_soil(k,isi)
               !---------------------------------------------------------------------------!
               !   Soil water.  Ensuring that the initial condition is within the accept-  !
               ! able range.                                                               !
               !---------------------------------------------------------------------------!
               csite%soil_water(k,ipa) = max(soil(ntext)%soilcp                            &
                                            ,min(soil(ntext)%slmsts                        &
                                                ,leaf_g(ifm)%soil_water(k,ix,iy,ilp) ) )

               !---------------------------------------------------------------------------!
               !   Soil temperature and liquid fraction. Simply use what we found a few    !
               ! lines above.                                                              !
               !---------------------------------------------------------------------------!
               csite%soil_tempk(k,ipa)   = lsoil_temp(k)
               csite%soil_fracliq(k,ipa) = lsoil_fliq(k)
               fice = 1.-lsoil_fliq(k)
               
               
               !---------------------------------------------------------------------------!
               !   Soil energy. Now that temperature, moisture and liquid partition are    !
               ! set, simply use the definition of internal energy to find it.             !
               !---------------------------------------------------------------------------!
               csite%soil_energy(k,ipa) = cmtl2uext( soil(ntext)%slcpd                     &
                                                   , csite%soil_water  (k,ipa) * wdns      &
                                                   , csite%soil_tempk  (k,ipa)             &
                                                   , csite%soil_fracliq(k,ipa)             )
            end do
            
            !----- Initialising surface snow/pond layers with nothing as default. ---------!
            csite%nlev_sfcwater(ipa) = 0
            do k=1,nzs
               csite%sfcwater_energy (k,ipa) = 0.
               csite%sfcwater_depth  (k,ipa) = 0.
               csite%sfcwater_mass   (k,ipa) = 0.
               csite%sfcwater_tempk  (k,ipa) = csite%soil_tempk(nzg,ipa)
               csite%sfcwater_fracliq(k,ipa) = csite%soil_fracliq(nzg,ipa)
            end do




            ntext = cpoly%ntext_soil(nzg,isi)
            ksn   = csite%nlev_sfcwater(ipa)
            ksnw1 = max(ksn,1)
            call ed_grndvap(ksn,ntext,csite%soil_water(nzg,ipa),csite%soil_tempk(nzg,ipa)  &
                           ,csite%soil_fracliq(nzg,ipa),csite%sfcwater_tempk(ksnw1,ipa)    &
                           ,csite%sfcwater_fracliq(ksnw1,ipa),csite%can_prss(ipa)          &
                           ,csite%can_shv(ipa),csite%ground_shv(ipa)                       &
                           ,csite%ground_ssh(ipa),csite%ground_temp(ipa)                   &
                           ,csite%ground_fliq(ipa),csite%ggsoil(ipa))

         end do patchloop
      end do siteloop
   end do polyloop
  
   !----- Allocate the scratch arrays. ----------------------------------------------------!
   deallocate (lsoil_text)
   deallocate (lsoil_temp)
   deallocate (lsoil_fliq)

   return
end subroutine leaf2ed_soil_moist_energy
!==========================================================================================!
!==========================================================================================!
