!==========================================================================================!
!==========================================================================================!
!     This subroutine initialises several ED fields that depend on atmospheric initial     !
! conditions plus some soil parameters.                                                    !
!------------------------------------------------------------------------------------------!
subroutine ed_init_atm()
   use stable_cohorts
   use update_derived_props_module
   use ed_misc_coms          , only : runtype                & ! intent(in)
                                    , ibigleaf               & ! intent(in)
                                    , ied_init_mode          ! ! intent(in)
   use ed_state_vars         , only : edtype                 & ! structure
                                    , polygontype            & ! structure
                                    , sitetype               & ! structure
                                    , patchtype              & ! structure
                                    , edgrid_g               ! ! structure
   use soil_coms             , only : soil_rough             & ! intent(in)
                                    , isoilstateinit         & ! intent(in)
                                    , soil                   & ! intent(in)
                                    , slmstr                 & ! intent(in)
                                    , stgoff                 & ! intent(in)
                                    , ed_soil_idx2water      & ! function
                                    , matric_potential       ! ! function
   use consts_coms           , only : wdns                   & ! intent(in)
                                    , t3ple                  ! ! intent(in)
   use grid_coms             , only : nzs                    & ! intent(in)
                                    , nzg                    & ! intent(in)
                                    , ngrids                 ! ! intent(in)
   use fuse_fiss_utils       , only : fuse_patches           & ! subroutine
                                    , rescale_patches        & ! subroutine
                                    , fuse_cohorts           & ! subroutine
                                    , terminate_cohorts      & ! subroutine
                                    , split_cohorts          ! ! subroutine
   use ed_therm_lib          , only : calc_veg_hcap          & ! subroutine
                                    , ed_grndvap             ! ! subroutine
   use canopy_layer_coms     , only : canstr                 & ! intent(out)
                                    , alloc_canopy_layer_mbs ! ! subroutine
   use therm_lib             , only : thetaeiv               & ! function
                                    , vpdefil                & ! function
                                    , idealdenssh            & ! function
                                    , qslif                  & ! function
                                    , reducedpress           & ! function
                                    , press2exner            & ! function
                                    , extheta2temp           & ! function
                                    , cmtl2uext              ! ! function
   use met_driver_coms       , only : met_driv_state         ! ! structure
   use canopy_struct_dynamics, only : canopy_turbulence      ! ! subroutine
   use budget_utils          , only : update_budget          ! ! subroutine
   !$ use omp_lib
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   type(edtype)        , pointer  :: cgrid
   type(polygontype)   , pointer  :: cpoly
   type(met_driv_state), pointer  :: cmet
   type(sitetype)      , pointer  :: csite
   type(patchtype)     , pointer  :: cpatch
   integer                        :: igr
   integer                        :: ipy
   integer                        :: isi
   integer                        :: ipa
   integer                        :: ico
   integer                        :: k
   integer                        :: nsoil
   integer                        :: nls
   integer                        :: nlsw1
   integer                        :: ncohorts
   integer                        :: npatches
   integer                        :: ping
   real                           :: site_area_i
   real                           :: poly_area_i
   real                           :: poly_lai
   real                           :: poly_nplant
   real                           :: elim_nplant
   real                           :: elim_lai
   real                           :: can_exner
   real                           :: rvaux
   integer                        :: ibuff
   integer                        :: nbuff
   !----- Add the MPI common block. -------------------------------------------------------!
#if defined(RAMS_MPI)
   include 'mpif.h'
#endif
   !---------------------------------------------------------------------------------------!

   !----- This is just any integer used to control the MPI sending/receiving tools. -------!
   ping = 6
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      This is as good a place as any to initialize canopy variables, particularly      !
   ! those that require separate buffers.                                                  !
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Initialize the canopy structure arrays.                                            !
   !---------------------------------------------------------------------------------------!
   nbuff = 1
   !$ nbuff= OMP_get_max_threads()
   allocate(canstr(nbuff))
   do ibuff=1,nbuff
      call alloc_canopy_layer_mbs(canstr(ibuff))
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   gridloop: do igr = 1,ngrids
      
      cgrid => edgrid_g(igr)

      !------------------------------------------------------------------------------------!
      !     If this is a standard ED2 restart, we will read these fields in from a history !
      ! file and therefore not worry about setting them here.                              !
      !------------------------------------------------------------------------------------!
      if(trim(runtype) == 'HISTORY' )return
      !------------------------------------------------------------------------------------!

      !----- Loop over polygons, sites and patches. ---------------------------------------!
      polyloop1: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop1: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)
            cmet  => cpoly%met(isi)

            patchloop1: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)


               !---------------------------------------------------------------------------!
               !      This first call is just to have the vegetation height so we can      !
               ! compute the initial canopy pressure...  It must be called again to have   !
               ! the storage right.                                                        !
               !---------------------------------------------------------------------------!
               call update_patch_derived_props(csite,ipa)
               !---------------------------------------------------------------------------!

               !---------------------------------------------------------------------------!
               !      Initialise some thermodynamic properties of the canopy air space.    !
               !---------------------------------------------------------------------------!
               csite%can_theta(ipa) = cmet%atm_theta
               csite%can_shv  (ipa) = cmet%atm_shv
               csite%can_co2  (ipa) = cmet%atm_co2
               csite%can_prss (ipa) = reducedpress(cmet%prss,cmet%atm_theta,cmet%atm_shv   &
                                                  ,cmet%geoht,csite%can_theta(ipa)         &
                                                  ,csite%can_shv(ipa),csite%can_depth(ipa))
               can_exner            = press2exner(csite%can_prss(ipa))
               csite%can_temp (ipa) = extheta2temp(can_exner,csite%can_theta(ipa))
               csite%can_temp_pv(ipa)=csite%can_temp(ipa)
               rvaux                = csite%can_shv(ipa) / (1. - csite%can_shv(ipa))
               

               csite%can_theiv(ipa) = thetaeiv(csite%can_theta(ipa),csite%can_prss(ipa)    &
                                              ,csite%can_temp(ipa),rvaux,rvaux)
               csite%can_vpdef(ipa) = vpdefil (csite%can_prss(ipa),csite%can_temp(ipa)     &
                                              ,csite%can_shv (ipa),.true.)
               csite%can_rhos (ipa) = idealdenssh(csite%can_prss(ipa)                      &
                                                 ,csite%can_temp(ipa),csite%can_shv(ipa))

               !----- Initialise the ground radiation parameters. -------------------------!
               csite%rshort_g(ipa) = 0.
               csite%par_g   (ipa) = 0.
               csite%rlong_g (ipa) = 0.

               csite%rough(ipa) = soil_rough

               !---------------------------------------------------------------------------!
               !      This value temporarily functions as a flag.  Please, don't change it !
               ! here.  The correct value will be assigned later on this subroutine.       !
               !---------------------------------------------------------------------------!
               csite%soil_tempk(1,ipa) = -100.0

               cohortloop1: do ico = 1,cpatch%ncohorts
                  !------------------------------------------------------------------------!
                  !      Here we initialise both the leaf and wood properties, assuming    !
                  ! thermal equilibrium with the canopy air space and no intercepted       !
                  ! water sitting on top of leaves and branches.                           !
                  !------------------------------------------------------------------------!
                  cpatch%leaf_water   (ico) = 0.0
                  cpatch%wood_water   (ico) = 0.0
                  cpatch%leaf_temp    (ico) = csite%can_temp   (ipa)
                  cpatch%wood_temp    (ico) = csite%can_temp   (ipa)
                  cpatch%leaf_temp_pv (ico) = csite%can_temp_pv(ipa)
                  cpatch%wood_temp_pv (ico) = csite%can_temp_pv(ipa)
                  cpatch%leaf_vpdef   (ico) = csite%can_vpdef  (ipa)
                  if (csite%can_temp(ipa) == t3ple) then
                     cpatch%leaf_fliq   (ico) = 0.5
                     cpatch%wood_fliq   (ico) = 0.5
                  elseif (csite%can_temp(ipa) > t3ple) then
                     cpatch%leaf_fliq   (ico) = 1.0
                     cpatch%wood_fliq   (ico) = 1.0
                  else
                     cpatch%leaf_fliq   (ico) = 0.0
                     cpatch%wood_fliq   (ico) = 0.0
                  end if
                  
                  
                  call calc_veg_hcap( cpatch%bleaf     (ico) , cpatch%bdead    (ico)       &
                                    , cpatch%bsapwooda (ico) , cpatch%nplant   (ico)       &
                                    , cpatch%pft       (ico) , cpatch%leaf_hcap(ico)       &
                                    , cpatch%wood_hcap (ico) )

                  cpatch%leaf_energy (ico) = cmtl2uext( cpatch%leaf_hcap   (ico)           &
                                                      , cpatch%leaf_water  (ico)           &
                                                      , cpatch%leaf_temp   (ico)           &
                                                      , cpatch%leaf_fliq   (ico) )
                  cpatch%wood_energy (ico) = cmtl2uext( cpatch%wood_hcap   (ico)           &
                                                      , cpatch%wood_water  (ico)           &
                                                      , cpatch%wood_temp   (ico)           &
                                                      , cpatch%wood_fliq   (ico) )

                  call is_resolvable(csite,ipa,ico)

                  !----- Initialise the leaf surface and intercellular properties. --------!
                  cpatch%lsfc_shv_open(ico)   = cmet%atm_shv
                  cpatch%lsfc_shv_closed(ico) = cmet%atm_shv
                  cpatch%lsfc_co2_open(ico)   = cmet%atm_co2
                  cpatch%lsfc_co2_closed(ico) = cmet%atm_co2
                  cpatch%lint_co2_open(ico)   = cmet%atm_co2
                  cpatch%lint_co2_closed(ico) = cmet%atm_co2
                  !------------------------------------------------------------------------!


                  !------------------------------------------------------------------------!
                  !      The intercellular specific humidity is assumed to be at           !
                  ! saturation.                                                            !
                  !------------------------------------------------------------------------!
                  cpatch%lint_shv(ico) = qslif(csite%can_prss(ipa),cpatch%leaf_temp(ico))
                  !------------------------------------------------------------------------!
               end do cohortloop1
            end do patchloop1
         end do siteloop1
      end do polyloop1


      !------------------------------------------------------------------------------------!
      !      Initialize remaining soil properties.                                         !
      !------------------------------------------------------------------------------------!
      select case(isoilstateinit)
      case (1)
         !---------------------------------------------------------------------------------!
         !      Initialize soil moisture, temperature, etc. from file specified in the     !
         ! namelist (the ED2 specific file.)                                               !
         !---------------------------------------------------------------------------------!
#if defined(RAMS_MPI)
         if (mynum /= 1) then
            call MPI_Recv(ping,1,MPI_INTEGER,recvnum,92,MPI_COMM_WORLD,MPI_STATUS_IGNORE   &
                         ,ierr)
         end if
#endif

         call read_soil_moist_temp(cgrid)

#if defined(RAMS_MPI)
         if (mynum < nnodetot) then
            call MPI_Send(ping,1,MPI_INTEGER,sendnum,92,MPI_COMM_WORLD,ierr)
         end if
#endif

      case (2)
#if defined(COUPLED)
         !----- Use the soil moisture and energy from LEAF-3 initialisation ---------------!
         call leaf2ed_soil_moist_energy(cgrid,igr)
#else
         !----- This is an offline run, there is no possibility to copy from LEAF-3. ------!
         write(unit=*,fmt='(a,1x,i5)') ' + ISOILSTATEINIT = ',isoilstateinit
         call fatal_error(' Invalid ISOILSTATEINIT in ED-2.1','ed_init_atm'                &
                         ,'ed_init_atm.F90')
#endif
      end select

      !------------------------------------------------------------------------------------!
      !      Do a simple, uniform initialization or take care of missing reanalysis        !
      ! points.                                                                            !
      !------------------------------------------------------------------------------------!
      polyloop2: do ipy = 1,cgrid%npolygons
         cpoly => cgrid%polygon(ipy)

         siteloop2: do isi = 1,cpoly%nsites
            csite => cpoly%site(isi)

            patchloop2: do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)

               if (csite%soil_tempk(1,ipa) == -100.0 .or. isoilstateinit == 0              &
                   .or. ied_init_mode == 7 ) then  ! The third is redundant, but keep it

                  groundloop2: do k = 1, nzg
                     nsoil=cpoly%ntext_soil(k,isi)

                     !----- Find the initial temperature. ---------------------------------!
                     csite%soil_tempk(k,ipa) = csite%can_temp(ipa) + stgoff(k)
                     !---------------------------------------------------------------------!

                     !------ Find the soil liquid fraction based on the temperature. ------!
                     if (csite%soil_tempk(k,ipa) > t3ple) then
                        nsoil=cpoly%ntext_soil(k,isi)
                        csite%soil_fracliq(k,ipa) = 1.0
                     elseif (csite%soil_tempk(k,ipa) < t3ple) then
                        csite%soil_fracliq(k,ipa) = 0.0
                     else
                        csite%soil_fracliq(k,ipa) = 0.5
                     end if
                     !---------------------------------------------------------------------!


                     !---------------------------------------------------------------------!
                     !    Initialise soil moisture and internal energy.                    !
                     !---------------------------------------------------------------------!
                     if (ied_init_mode /= 7) then
                        csite%soil_water(k,ipa) = ed_soil_idx2water(slmstr(k),nsoil)
                     end if
                     csite%soil_energy  (k,ipa) = cmtl2uext( soil(nsoil)%slcpd             &
                                                           , csite%soil_water(k,ipa)*wdns  &
                                                           , csite%soil_tempk(k,ipa)       &
                                                           , csite%soil_fracliq(k,ipa)     )
                     csite%soil_mstpot  (k,ipa) = matric_potential(nsoil                   &
                                                                 ,csite%soil_water(k,ipa)  )
                     !---------------------------------------------------------------------!
                  end do groundloop2

                  !----- Initial condition is with no snow/pond. --------------------------!
                  csite%nlev_sfcwater(ipa)    = 0.
                  csite%snowfac(ipa)          = 0.
                  csite%total_sfcw_depth(ipa) = 0.
                  snowloop2: do k=1,nzs
                     csite%sfcwater_energy (k,ipa) = 0.
                     csite%sfcwater_depth  (k,ipa) = 0.
                     csite%sfcwater_mass   (k,ipa) = 0.
                     csite%sfcwater_tempk  (k,ipa) = csite%soil_tempk(nzg,ipa)
                     csite%sfcwater_fracliq(k,ipa) = csite%soil_fracliq(nzg,ipa)
                  end do snowloop2
               
                  !----- Compute patch-level LAI, vegetation height, and roughness. -------!
                  call update_patch_derived_props(csite,ipa)

                  nsoil = cpoly%ntext_soil(nzg,isi)
                  nls   = csite%nlev_sfcwater(ipa)
                  nlsw1 = max(nls,1)
                  call ed_grndvap(nls,nsoil,csite%soil_water(nzg,ipa)                      &
                                 ,csite%soil_tempk(nzg,ipa),csite%soil_fracliq(nzg,ipa)    &
                                 ,csite%sfcwater_tempk(nlsw1,ipa)                          &
                                 ,csite%snowfac(ipa)                                       &
                                 ,csite%can_prss(ipa),csite%can_shv(ipa)                   &
                                 ,csite%ground_shv(ipa),csite%ground_ssh(ipa)              &
                                 ,csite%ground_temp(ipa),csite%ground_fliq(ipa)            &
                                 ,csite%ggsoil(ipa))
               else
                  !----- Compute patch-level LAI, vegetation height, and roughness. -------!
                  call update_patch_derived_props(csite,ipa)

                  nsoil = cpoly%ntext_soil(nzg,isi)
                  nls   = csite%nlev_sfcwater(ipa)
                  nlsw1 = max(nls,1)
                  call ed_grndvap(nls,nsoil,csite%soil_water(nzg,ipa)                      &
                                 ,csite%soil_tempk(nzg,ipa),csite%soil_fracliq(nzg,ipa)    &
                                 ,csite%sfcwater_tempk(nlsw1,ipa)                          &
                                 ,csite%snowfac(ipa)                                       &
                                 ,csite%can_prss(ipa),csite%can_shv(ipa)                   &
                                 ,csite%ground_shv(ipa),csite%ground_ssh(ipa)              &
                                 ,csite%ground_temp(ipa),csite%ground_fliq(ipa)            &
                                 ,csite%ggsoil(ipa))
               end if

               !----- Initialise vegetation wind and turbulence parameters. ---------------!
               call canopy_turbulence(cpoly,isi,ipa)

               !----- Compute the storage terms for CO2, energy, and water budgets. -------!
               call update_budget(csite,cpoly%lsl(isi),ipa)

            end do patchloop2

            !----- Compute basal area and AGB profiles. -----------------------------------!
            call update_site_derived_props(cpoly, 0, isi)
         end do siteloop2
      end do polyloop2

      !----- Compute some derived properties at the polygon level. ------------------------!
      call update_polygon_derived_props(cgrid)

      !----- Fuse similar patches to speed up the run. ------------------------------------!
      select case(ibigleaf)
      case (0)
         !---------------------------------------------------------------------------------!
         !    Size and age structure.  Start by fusing similar patches.                    !
         !---------------------------------------------------------------------------------!
         call fuse_patches(cgrid,igr,.true.)
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !    Loop over all polygons/sites/patches, and fuse/split/terminate cohorts as    !
         ! needed.                                                                         !
         !---------------------------------------------------------------------------------!
         polyloop3: do ipy = 1,cgrid%npolygons
            ncohorts     = 0
            npatches     = 0
            poly_lai     = 0.0
            poly_nplant  = 0.0

            cpoly => cgrid%polygon(ipy)
            poly_area_i = 1./sum(cpoly%area(:))

            siteloop3: do isi = 1,cpoly%nsites
               csite => cpoly%site(isi)
               site_area_i = 1./sum(csite%area(:))

               patchloop3: do ipa = 1,csite%npatches
                  npatches = npatches + 1
                  cpatch => csite%patch(ipa)

                  if (cpatch%ncohorts > 0) then
                     call fuse_cohorts(csite,ipa,cpoly%lsl(isi),.true.)
                     call terminate_cohorts(csite,ipa,elim_nplant,elim_lai)
                     call split_cohorts(cpatch, cpoly%green_leaf_factor(:,isi))
                  end if

                  cohortloop3: do ico = 1,cpatch%ncohorts
                     ncohorts=ncohorts+1
                     poly_lai    = poly_lai + cpatch%lai(ico) * csite%area(ipa)            &
                                            * cpoly%area(isi) * site_area_i * poly_area_i
                     poly_nplant = poly_nplant + cpatch%nplant(ico) * csite%area(ipa)      &
                                               * cpoly%area(isi) * site_area_i             &
                                               * poly_area_i
                  end do cohortloop3
               end do patchloop3
            end do siteloop3

            write(unit = *                                                                 &
                 ,fmt  = '(2(a,1x,i6,1x),2(a,1x,f9.4,1x),2(a,1x,f7.2,1x),2(a,1x,i4,1x))')  &
                     'Grid:',igr,'Poly:',ipy,'Lon:',cgrid%lon(ipy),'Lat: ',cgrid%lat(ipy)  &
                    ,'Nplants:',poly_nplant,'Avg. LAI:',poly_lai                           &
                    ,'NPatches:',npatches,'NCohorts:',ncohorts
         end do polyloop3
         !---------------------------------------------------------------------------------!


      case (1)
         !---------------------------------------------------------------------------------!
         !     Big leaf.  No need to do anything, just print the banner.                   !
         !---------------------------------------------------------------------------------!
         polyloop4: do ipy = 1,cgrid%npolygons
            ncohorts     = 0
            npatches     = 0
            poly_lai     = 0.0
            poly_nplant  = 0.0
            
            cpoly => cgrid%polygon(ipy)
            poly_area_i = 1./sum(cpoly%area(:))
            
            siteloop4: do isi = 1,cpoly%nsites
               csite => cpoly%site(isi)
               site_area_i = 1./sum(csite%area(:))
               
               !call rescale_patches(csite)
               
               patchloop4: do ipa = 1,csite%npatches
                  npatches = npatches + 1
                  cpatch => csite%patch(ipa)
               
                  cohortloop4: do ico = 1,cpatch%ncohorts
                     ncohorts=ncohorts+1
                     poly_lai    = poly_lai + cpatch%lai(ico) * csite%area(ipa)            &
                                         * cpoly%area(isi) * site_area_i * poly_area_i
                     poly_nplant = poly_nplant + cpatch%nplant(ico) * csite%area(ipa)      &
                                           * cpoly%area(isi) * site_area_i * poly_area_i
                  end do cohortloop4
               end do patchloop4
            end do siteloop4
            
            write( unit = *                                                                &
                 , fmt  = '(2(a,1x,i6,1x),2(a,1x,f9.4,1x),2(a,1x,f7.2,1x),2(a,1x,i4,1x))') &
                'Grid:',igr,'Poly:',ipy,'Lon:',cgrid%lon(ipy),'Lat: ',cgrid%lat(ipy)       &
               ,'Nplants:',poly_nplant,'Avg. LAI:',poly_lai                                &
               ,'NPatches:',npatches,'NCohorts:',ncohorts
            end do polyloop4
         !---------------------------------------------------------------------------------!
         end select
      !------------------------------------------------------------------------------------!
   end do gridloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine ed_init_atm
!==========================================================================================!
!==========================================================================================!
