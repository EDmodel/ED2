!==========================================================================================!
!==========================================================================================!
subroutine ed_init_coup_atm
  
  use ed_misc_coms,  only: ied_init_mode,runtype
  use ed_state_vars, only: edtype,polygontype,sitetype,patchtype,edgrid_g
  use soil_coms,     only: soil_rough, isoilstateinit, soil, slmstr,dslz
  use rconstants,    only: tsupercool, cliqvlme, cicevlme, t3ple
  use grid_coms,      only: nzs, nzg, ngrids
  use fuse_fiss_utils, only: fuse_patches,fuse_cohorts,terminate_cohorts,split_cohorts
  use ed_node_coms, only: nnodetot,mynum,sendnum,recvnum
  use pft_coms,only : sla
  use ed_therm_lib,only : calc_hcapveg,ed_grndvap

  implicit none

  type(edtype)     ,pointer :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch
  integer :: igr,ipy,isi,ipa,ico
  integer :: k
  integer :: nsoil
  integer :: nls
  integer :: nlsw1
  integer :: ncohorts, npatches
  real    :: site_area_i, poly_area_i
  real    :: poly_lai,poly_nplant
  real    :: elim_nplant,elim_lai
  integer :: ix,iy
  real    :: surface_temp,surface_fliq
  integer, parameter :: harvard_override = 0
  include 'mpif.h'
  integer :: ping,ierr
  ping = 6 ! Just any rubbish for MPI Send/Recv

  ! This subroutine fills the ED2 fields which depend on current 
  ! atmospheric conditions.

  do igr = 1,ngrids
     
     cgrid => edgrid_g(igr)

     ! If this is a standard ED2 restart, we will read these fields in from 
     ! a history file and therefore not worry about setting them here.
     if(ied_init_mode == 4 .or. trim(runtype) == 'HISTORY' )return
     
     ! Loop over polygons, sites and patches
     
     do ipy = 1,cgrid%npolygons
        
        cpoly => cgrid%polygon(ipy)
        
        do isi = 1,cpoly%nsites
           
           csite => cpoly%site(isi)

           do ipa = 1,csite%npatches
              
              cpatch => csite%patch(ipa)

              csite%can_temp(ipa) =   cpoly%met(isi)%atm_tmp
              csite%can_shv(ipa)  =   cpoly%met(isi)%atm_shv
              csite%can_co2(ipa)  =   cpoly%met(isi)%atm_co2

              ! Initialize stars
              csite%tstar(ipa)  = 0.
              csite%ustar(ipa)  = 0.
              csite%qstar(ipa)  = 0.
              csite%cstar(ipa)  = 0.
              
              

              csite%rshort_g(ipa) = 0.0
              csite%rlong_g(ipa) = 0.0
              
              ! Initialize soil textural class.  Soil water, energy, etc. will
              ! be initialized in the next round of loops.
              do k = 1,nzg
                 csite%ntext_soil(k,ipa) = cpoly%ntext_soil(k,isi)
              enddo
              
              csite%rough(ipa) = soil_rough
              csite%soil_tempk(1,ipa) = -100.0 ! This value functions as a flag.  Do not 
                                               ! change it here. It will be changed below.
              csite%hcapveg(ipa) = 0.
              do ico = 1,cpatch%ncohorts
                 ! Initialize vegetation properties.
                 ! For now, set heat capacity for stability.

                 cpatch%veg_temp(ico)   = cpoly%met(isi)%atm_tmp
                 cpatch%veg_fliq(ico)   = 0.0
                 cpatch%veg_water(ico)  = 0.0
                 cpatch%hcapveg(ico)    = calc_hcapveg(cpatch%bleaf(ico)           &
                                                      ,cpatch%bdead(ico)           &
                                                      ,cpatch%balive(ico)          &
                                                      ,cpatch%nplant(ico)          &
                                                      ,cpatch%hite(ico)            &
                                                      ,cpatch%pft(ico)             &
                                                      ,cpatch%phenology_status(ico))
                 cpatch%veg_energy(ico) = cpatch%hcapveg(ico)*cpatch%veg_temp(ico)
                 csite%hcapveg(ipa) = csite%hcapveg(ipa) + cpatch%hcapveg(ico)
              end do
           
           enddo
           
        enddo
        
     enddo
     ! Initialize remaining soil properties.
     if(isoilstateinit == 1)then
        ! Initialize soil moisture, temperature, etc. from file specified in 
        ! the ED_NL.

        if (mynum    /= 1) &
           call MPI_Recv(ping,1,MPI_INTEGER,recvnum,92,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

        call read_soil_moist_temp(cgrid)

        if (mynum     < nnodetot) &
           call MPI_Send(ping,1,MPI_INTEGER,sendnum,92,MPI_COMM_WORLD,ierr)
     elseif (isoilstateinit == 2) then
        ! Use the soil moisture and energy from LEAF-3 initialisation ---------------------!
        call leaf2ed_soil_moist_energy(cgrid,igr)
     end if

     ! Do a simple, uniform initialization or take care of 
     ! missing reanalysis points
     
     ! Loop over polygons, sites, patches.
     do ipy = 1,cgrid%npolygons
        
        cpoly => cgrid%polygon(ipy)
        
        do isi = 1,cpoly%nsites
    
           csite => cpoly%site(isi)

           do ipa = 1,csite%npatches
              
              cpatch => csite%patch(ipa)
              
              if(csite%soil_tempk(1,ipa) == -100.0 .or. isoilstateinit == 0)then
                 
                 csite%soil_tempk(1:nzg,ipa) = csite%can_temp(ipa)
                 
                 if(csite%can_temp(ipa) > t3ple)then
                    do k = 1, nzg
                       nsoil=csite%ntext_soil(k,ipa)
                       csite%soil_fracliq(k,ipa) = 1.0
                       csite%soil_water(k,ipa) = max( soil(nsoil)%soilcp                   &
                                                    , slmstr(k) * soil(nsoil)%slmsts)
                       csite%soil_energy(k,ipa) = soil(nsoil)%slcpd                        &
                                                * csite%soil_tempk(k,ipa)                  &
                                                + csite%soil_water(k,ipa) * cliqvlme       &
                                                * (csite%soil_tempk(k,ipa) - tsupercool)
                    end do
                 else
                    do k = 1, nzg
                       nsoil=csite%ntext_soil(k,ipa)
                       csite%soil_fracliq(k,ipa) = 0.0
                       csite%soil_water(k,ipa) = max( soil(nsoil)%soilcp                   &
                                                    , slmstr(k) * soil(nsoil)%slmsts)
                       csite%soil_energy(k,ipa) = soil(nsoil)%slcpd                        &
                                                * csite%soil_tempk(k,ipa)                  &
                                                + csite%soil_water(k,ipa)                  &
                                                * cicevlme * csite%soil_tempk(k,ipa)
                    end do
                 end if

                 !----- Initial condition is with no snow/pond. ---------------------------!
                 csite%nlev_sfcwater(ipa) = 0
                 do k=1,nzs
                    csite%sfcwater_energy (k,ipa) = 0.
                    csite%sfcwater_depth  (k,ipa) = 0.
                    csite%sfcwater_mass   (k,ipa) = 0.
                    csite%sfcwater_tempk  (k,ipa) = csite%soil_tempk(nzg,ipa)
                    csite%sfcwater_fracliq(k,ipa) = csite%sfcwater_fracliq(nzg,ipa)
                 end do

                 nls   = csite%nlev_sfcwater(ipa)
                 nlsw1 = max(nls,1)
                 
                 call ed_grndvap(nls,                    &
                      csite%ntext_soil       (nzg,ipa),  &
                      csite%soil_water       (nzg,ipa),  &
                      csite%soil_energy      (nzg,ipa),  &
                      csite%sfcwater_energy(nlsw1,ipa),  &
                      cpoly%met(isi)%rhos,  &
                      csite%can_shv(ipa),  &
                      csite%ground_shv(ipa),  &
                      csite%surface_ssh(ipa), &
                      surface_temp,&
                      surface_fliq)
              endif
              
              ! Compute patch-level LAI, vegetation height, and roughness
              call update_patch_derived_props(csite, cpoly%lsl(isi), cpoly%met(isi)%rhos, ipa)
              

           enddo
           
           ! Compute basal area and AGB profiles.
           call update_site_derived_props(cpoly, 0, isi)
           
        enddo
        
        
        
     enddo
     
     call update_polygon_derived_props(cgrid)


     call fuse_patches(cgrid,igr)

     do ipy = 1,cgrid%npolygons
        ncohorts     = 0
        npatches     = 0
        poly_lai     = 0.0
        poly_nplant  = 0.0

        cpoly => cgrid%polygon(ipy)
        poly_area_i = 1./sum(cpoly%area(:))

        do isi = 1,cpoly%nsites

           csite => cpoly%site(isi)
           site_area_i = 1./sum(csite%area(:))
           do ipa = 1,csite%npatches
              npatches = npatches + 1
              cpatch => csite%patch(ipa)
              if (cpatch%ncohorts > 0) then
                 call fuse_cohorts(csite,ipa,cpoly%green_leaf_factor(:,isi),cpoly%lsl(isi))
                 call terminate_cohorts(csite,ipa,elim_nplant,elim_lai)
                 call split_cohorts(cpatch, cpoly%green_leaf_factor(:,isi), cpoly%lsl(isi))
              end if
              do ico = 1,cpatch%ncohorts
                 ncohorts=ncohorts+1
                 poly_lai    = poly_lai + cpatch%lai(ico) * csite%area(ipa)                &
                                        * cpoly%area(isi) * site_area_i * poly_area_i
                 poly_nplant = poly_nplant + cpatch%nplant(ico) * csite%area(ipa)          &
                                           * cpoly%area(isi) * site_area_i * poly_area_i
              end do
           end do
        end do
        write(unit=*,fmt='(2(a,1x,i4,1x),2(a,1x,f9.4,1x),2(a,1x,f7.2,1x),2(a,1x,i4,1x))')  &
            'Grid:',igr,'Poly:',ipy,'Lon:',cgrid%lon(ipy),'Lat: ',cgrid%lat(ipy)           &
           ,'Nplants:',poly_nplant,'Avg. LAI:',poly_lai                                    &
           ,'NPatches:',npatches,'NCohorts:',ncohorts
     end do
  end do

  return
end subroutine ed_init_coup_atm
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_init_radiation

  use mem_radiate,only:radiate_g
  use mem_leaf,only:leaf_g
  use mem_grid,only:ngrids
  use rconstants,only:stefan
  use ed_state_vars,only:edgrid_g,edtype

  implicit none
  
  type(edtype),pointer :: cgrid
  integer :: igr,ipy
  integer :: ix,iy

  do igr=1,ngrids

     cgrid => edgrid_g(igr)
     
     ! Initilialize the upwelling longwave radiation and albedo
     ! --------------------------------------------------------
     
     ! First, do a catch all using sst and crude albedo
     radiate_g(igr)%rlongup  = 0.98 * stefan * leaf_g(igr)%seatp**4
     radiate_g(igr)%albedt   = 0.3
     
     ! Then use the air temperature to initialize a brightness temperature
     ! This is not ideal, but the air temperature just initialized the land-surface
     ! variables.  These radiation parameters are just to prevent the model from
     ! crashing.
     do ipy = 1,cgrid%npolygons
        
        ix = cgrid%ilon(ipy)
        iy = cgrid%ilat(ipy)
        
        radiate_g(igr)%rlongup(ix,iy) = 0.97 * stefan * cgrid%met(ipy)%atm_tmp**4
        radiate_g(igr)%albedt(ix,iy)  = 0.3
        
     enddo

  enddo

  return
end subroutine ed_init_radiation
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine leaf2ed_soil_moist_energy(cgrid,ifm)
   use ed_state_vars, only : edtype       & ! structure
                           , polygontype  & ! structure
                           , sitetype     & ! structure
                           , patchtype    ! ! structure
   use grid_coms    , only : nzg          & ! intent(in)
                           , nzs          ! ! intent(in)
   use ed_therm_lib , only : ed_grndvap   ! ! subroutine
   use therm_lib    , only : qwtk8        & ! subroutine
                           , qwtk         ! ! subroutine
   use rconstants   , only : wdns         & ! intent(in)
                           , tsupercool   & ! intent(in)
                           , cicevlme     & ! intent(in)
                           , cliqvlme     ! ! intent(in)
   use mem_leaf     , only : leaf_g       ! ! structure
   use leaf_coms    , only : slcpd        ! ! intent(in)
   use soil_coms    , only : soil         ! ! intent(in)
   
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(edtype)      , target     :: cgrid  ! Alias for current ED grid
   integer           , intent(in) :: ifm
   !----- Local variables -----------------------------------------------------------------!
   type(polygontype) , pointer    :: cpoly           ! Alias for current polygon
   type(sitetype)    , pointer    :: csite           ! Alias for current site
   type(patchtype)   , pointer    :: cpatch          ! Alias for current patch
   integer                        :: ntext           ! Alias for ED-2 soil texture class
   integer, dimension(nzg)        :: lsoil_text      ! LEAF-3 soil texture class
   integer                        :: ipy,isi,ipa,ico ! Counters for all structures
   integer                        :: k               ! Counter for soil layers
   integer                        :: ix,iy           ! Counter for lon/lat
   integer                        :: ksn, ksnw1      ! Alias for # of pond/snow layers
   real, dimension(nzg)           :: lsoil_temp      ! LEAF-3 soil temperature
   real, dimension(nzg)           :: lsoil_fliq      ! LEAF-3 soil liquid fraction
   real                           :: surface_temp    ! Scratch variable for ed_grndvap
   real                           :: surface_fliq    ! Scratch variable for ed_grndvap
   real                           :: fice            ! soil ice fraction
   !---------------------------------------------------------------------------------------!



   !----- Loop over land points -----------------------------------------------------------!
   polyloop: do ipy=1,cgrid%npolygons
      ix = cgrid%ilon(ipy)
      iy = cgrid%ilat(ipy)
      
      !------------------------------------------------------------------------------------!
      !    Determining initial soil temperature and liquid fraction.  The reason we find   !
      ! this for LEAF-3 instead of simply copying the soil energy and water to ED-2 is     !
      ! that depending on the way the user set up his or her RAMSIN, soil types may not    !
      ! match and this could put ED-2.1 in an inconsistent initial state.                  !
      !------------------------------------------------------------------------------------!
      do k=1,nzg
         lsoil_text(k) =nint(leaf_g(ifm)%soil_text(k,ix,iy,2))
         call qwtk(leaf_g(ifm)%soil_energy(k,ix,iy,2)                                      &
                  ,leaf_g(ifm)%soil_water(k,ix,iy,2)*wdns                                  &
                  ,slcpd(lsoil_text(k)),lsoil_temp(k),lsoil_fliq(k))
      end do

      cpoly => cgrid%polygon(ipy)

      !----- Loop over sites --------------------------------------------------------------!
      siteloop: do isi=1,cpoly%nsites
         csite => cpoly%site(isi)
         
         !----- Loop over patches ---------------------------------------------------------!
         patchloop: do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)
  
            do k=1,nzg
            
               ntext = csite%ntext_soil(k,ipa)
               !---------------------------------------------------------------------------!
               !   Soil water.  Ensuring that the initial condition is within the accept-  !
               ! able range.                                                               !
               !---------------------------------------------------------------------------!
               csite%soil_water(k,ipa) = max(soil(ntext)%soilcp                            &
                                            ,min(soil(ntext)%slmsts                        &
                                                ,leaf_g(ifm)%soil_water(k,ix,iy,2) ) )

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
               csite%soil_energy(k,ipa) = soil(ntext)%slcpd * csite%soil_tempk(k,ipa)      &
                                        + csite%soil_water(k,ipa)                          &
                                        * ( fice * cicevlme * csite%soil_tempk(k,ipa)      &
                                          + csite%soil_fracliq(k,ipa) * cliqvlme           &
                                          * (csite%soil_tempk(k,ipa) - tsupercool) )
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
            
            !----- Compute the ground properties ------------------------------------------!
            ksn   = csite%nlev_sfcwater(ipa)
            ksnw1 = max(ksn,1)
            call ed_grndvap(ksn,csite%ntext_soil(nzg,ipa),csite%soil_water(nzg,ipa)        &
                           ,csite%soil_energy(nzg,ipa),csite%sfcwater_energy(ksnw1,ipa)    &
                           ,cpoly%met(isi)%rhos,csite%can_shv(ipa),csite%ground_shv(ipa)   &
                           ,csite%surface_ssh(ipa),surface_temp,surface_fliq)
         end do patchloop
      end do siteloop
   end do polyloop
  
   return
end subroutine leaf2ed_soil_moist_energy
!==========================================================================================!
!==========================================================================================!
