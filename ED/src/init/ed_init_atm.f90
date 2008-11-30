
subroutine ed_init_atm_ar
  
  use misc_coms,     only: ied_init_mode,runtype
  use ed_state_vars, only: edtype,polygontype,sitetype,patchtype,edgrid_g
  use soil_coms,     only: soil_rough, isoilstateinit, soil, slmstr
  use consts_coms,    only: alli1000, cliq1000, cice1000, t3ple
  use grid_coms,      only: nzs, nzg, ngrids
  use fuse_fiss_utils_ar, only: fuse_patches_ar,fuse_cohorts_ar
  use ed_node_coms, only: nnodetot,mynum,sendnum,recvnum
  use pft_coms,only : sla
  use ed_therm_lib,only : update_veg_energy_ct
  
  
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
  integer :: ncohorts
  real    :: poly_lai,p_lai
  integer, parameter :: harvard_override = 0
  include 'mpif.h'
  integer :: ping,ierr
  ping = 6 ! Just any rubbish for MPI Send/Recv

  ! This subroutine fills the ED2 fields which depend on current 
  ! atmospheric conditions.

  do igr = 1,ngrids
     
     cgrid => edgrid_g(igr)
     
     ! First we need to update the meteorological fields.
     call update_met_drivers_array(cgrid)

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
              
              ! Initialize stars
              !csite%tstar(ipa)  = 0.
              !csite%ustar(ipa)  = 0.
              !csite%rstar(ipa)  = 0.
              !csite%cstar(ipa)  = 0.
              
              ! For now, choose heat/vapor capacities for stability
              csite%can_depth(ipa) = 30.0
              do k=1,nzs
                 csite%sfcwater_tempk(k,ipa) = t3ple   ! Set canopy temp triple point
                 csite%sfcwater_fracliq(k,ipa) = 1.0 ! Set to 100% liquid
              end do
              
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
              
              do ico = 1,cpatch%ncohorts

                 ! Initialize vegetation properties.
                 ! For now, set heat capacity for stability.

                 cpatch%veg_temp(ico)  = cpoly%met(isi)%atm_tmp
                 cpatch%veg_water(ico) = 0.0

                 call update_veg_energy_ct(cpatch,ico)

              enddo
           
           enddo
           
        enddo
        
     enddo
     ! Initialize remaining soil properties.
     if(isoilstateinit == 1)then
        ! Initialize soil moisture, temperature, etc. from file specified in 
        ! the ED_NL.
        if (nnodetot /= 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)
        if (mynum    /= 1) call MPI_Recv(ping,1,MPI_INTEGER,recvnum,864,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

        call read_soil_moist_temp_ar(cgrid)

        if (mynum     < nnodetot) call MPI_Send(ping,1,MPI_INTEGER,sendnum,864,MPI_COMM_WORLD,ierr)
        if (nnodetot /=        1) call MPI_Barrier(MPI_COMM_WORLD,ierr)

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
              
              if(csite%soil_tempk(1,ipa) == -100.0 .or. isoilstateinit /= 1)then
                 
                 csite%soil_tempk(1:nzg,ipa) = csite%can_temp(ipa)
                 
                 if(csite%can_temp(ipa) > t3ple)then
                    do k = 1, nzg
                       nsoil=csite%ntext_soil(k,ipa)
                       csite%soil_fracliq(k,ipa) = 1.0
                       csite%soil_water(k,ipa) = max(soil(nsoil)%soilcp,   &
                            slmstr(k) * soil(nsoil)%slmsts)
                       csite%soil_energy(k,ipa) = (csite%soil_tempk(k,ipa) - t3ple) *   &
                            (soil(nsoil)%slcpd + csite%soil_water(k,ipa) *   &
                            cliq1000) + csite%soil_water(k,ipa) * alli1000
                    enddo
                 else
                    do k = 1, nzg
                       nsoil=csite%ntext_soil(k,ipa)
                       csite%soil_fracliq(k,ipa) = 0.0
                       csite%soil_water(k,ipa) = max(soil(nsoil)%soilcp,             &
                            slmstr(k) * soil(nsoil)%slmsts)
                       csite%soil_energy(k,ipa) = (csite%soil_tempk(k,ipa) - t3ple) *   &
                            (soil(nsoil)%slcpd + csite%soil_water(k,ipa) * cice1000)
                    enddo
                 endif
              
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
                      csite%surface_ssh(ipa))
              endif
              
              ! Compute patch-level LAI, vegetation height, and roughness
              call update_patch_derived_props_ar(csite, cpoly%lsl(isi), cpoly%met(isi)%rhos, ipa)
              

           enddo
           
           ! Compute basal area and AGB profiles.
           call update_site_derived_props_ar(cpoly, 0, isi)
           
        enddo
        
        
        
     enddo
     
     call update_polygon_derived_props_ar(cgrid)

     ! Energy needs to be done after LAI and Hite are loaded
     call initialize_vegetation_energy(cgrid)


     call fuse_patches_ar(cgrid)
     do ipy = 1,cgrid%npolygons
        
        cpoly => cgrid%polygon(ipy)
        
        do isi = 1,cpoly%nsites
           
           csite => cpoly%site(isi)
           
           do ipa = 1,csite%npatches
              
              cpatch => csite%patch(ipa)

           enddo

        enddo
     enddo

     do ipy = 1,cgrid%npolygons
        
        ncohorts = 0
        poly_lai = 0.0
        
        cpoly => cgrid%polygon(ipy)
        
        do isi = 1,cpoly%nsites
           
           csite => cpoly%site(isi)
           
           do ipa = 1,csite%npatches
              
              cpatch => csite%patch(ipa)

              call fuse_cohorts_ar(csite,ipa,cpoly%green_leaf_factor(:,isi),cpoly%lsl(isi))
              
              do ico = 1,cpatch%ncohorts
                 ncohorts=ncohorts+1
                 poly_lai = poly_lai + cpatch%lai(ico) * csite%area(ipa)*cpoly%area(isi)
              enddo
              
           enddo
           
        enddo
        write(*,'(2(a,1x,i4,1x),2(a,1x,f9.4,1x),a,1x,f5.2,1x,a,1x,i4)')   &
            'Grid:',igr,'Poly:',ipy,'Lon:',cgrid%lon(ipy),'Lat: ',cgrid%lat(ipy),'Avg. LAI:',poly_lai,'NCohorts:',ncohorts

     enddo
  
  enddo

  ! Energy needs to be done after LAI and Hite are loaded
  call initialize_vegetation_energy(cgrid)

  return
end subroutine ed_init_atm_ar

!==========================================================================================!
!==========================================================================================!

subroutine update_derived_props(cgrid)
  ! Update some of the derived quantities (this may be redundant)
  use ed_state_vars, only: edtype,polygontype,sitetype
  implicit none
  type(edtype)      , target  :: cgrid
  type(polygontype) , pointer :: cpoly
  type(sitetype)    , pointer :: csite
  integer                     :: ipy, isi, ipa

   do ipy = 1,cgrid%npolygons
     
     cpoly => cgrid%polygon(ipy)
     
     do isi = 1,cpoly%nsites
        csite => cpoly%site(isi)
        do ipa = 1,csite%npatches
           call update_patch_derived_props_ar(csite, cpoly%lsl(isi), cpoly%met(isi)%rhos, ipa)
        enddo
        call update_site_derived_props_ar(cpoly, 0, isi)
     enddo
     call update_polygon_derived_props_ar(cgrid)
   enddo

   return
end subroutine update_derived_props
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine update_patch_derived_props_ar(csite, lsl, rhos, ipa)
  
  use ed_state_vars,only:sitetype,patchtype

  implicit none
  integer         , intent(in) :: ipa
  real            , intent(in) :: rhos
  integer         , intent(in) :: lsl
  type(sitetype)  , target     :: csite
  type(patchtype) , pointer    :: cpatch

  real                         :: norm_fac, ba
  integer                      :: ico
  real            , external   :: compute_water_storage_ar
  real            , external   :: compute_energy_storage_ar
  real            , external   :: compute_co2_storage_ar
  real            , external   :: ed_biomass

  ! call derived patch-level structural quantities.  These depend
  ! on the results from reproduction, which in turn depends on 
  ! structural growth results from all patches.


  ! Reset height
  csite%veg_height(ipa) = 0.0
  csite%lai(ipa)        = 0.0
  norm_fac              = 0.0
  csite%plant_ag_biomass(ipa) = 0.0
  cpatch => csite%patch(ipa)

  ! Loop over cohorts
  do ico = 1,cpatch%ncohorts
     
     ! Compute contribution to height
     ba = cpatch%nplant(ico) * cpatch%dbh(ico)**2
     norm_fac = norm_fac + ba
     csite%veg_height(ipa) = csite%veg_height(ipa) + cpatch%hite(ico) * ba
     
     ! Update LAI and AGB
     csite%lai(ipa) = csite%lai(ipa) + cpatch%lai(ico)
     csite%plant_ag_biomass(ipa)  = csite%plant_ag_biomass(ipa) +                          &
           ed_biomass(cpatch%bdead(ico),cpatch%balive(ico), cpatch%bleaf(ico)              &
                      ,cpatch%pft(ico), cpatch%hite(ico),cpatch%bstorage(ico))             &
           *cpatch%nplant(ico)           
  
  enddo
  
  if (csite%lai(ipa).lt.0.0) then
     print*,"STRANGE LAI, ncohorts:",cpatch%ncohorts,"dist type",csite%dist_type(ipa)
     do ico = 1,cpatch%ncohorts
        print*,cpatch%lai(ico),cpatch%bleaf(ico),cpatch%nplant(ico),&
             cpatch%pft(ico),cpatch%phenology_status(ico),cpatch%maintenance_costs(ico),&
             cpatch%bstorage(ico),cpatch%dbh(ico)
     enddo
     stop
  endif

  
  ! Update vegetation height
  if(norm_fac > 0.0)then
     csite%veg_height(ipa) = csite%veg_height(ipa) / norm_fac
  else
     ! this branch if there aren't any cohorts
     csite%veg_height(ipa) = 0.2
  endif
  csite%veg_rough(ipa) = 0.13 * csite%veg_height(ipa)
  

  csite%wbudget_initialstorage(ipa) = compute_water_storage_ar(csite, lsl, rhos, ipa)
  csite%ebudget_initialstorage(ipa) = compute_energy_storage_ar(csite, lsl, rhos, ipa)
  csite%co2budget_initialstorage(ipa) = compute_co2_storage_ar(csite, rhos, ipa)

  csite%cohort_count(ipa) = cpatch%ncohorts

  return
end subroutine update_patch_derived_props_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine update_site_derived_props_ar(cpoly, census_flag, isi)
  
  use ed_state_vars,only: polygontype,sitetype,patchtype

  use consts_coms,    only : pi1
  implicit none
  
  type(polygontype),target :: cpoly
  type(sitetype), pointer :: csite
  type(patchtype), pointer :: cpatch
  integer :: isi,ipa,ico
  real :: ba
  integer :: bdbh
  real, external :: ed_biomass
  integer, intent(in) :: census_flag

  cpoly%basal_area(:,:,isi) = 0.0
  cpoly%agb(:,:,isi) = 0.0
  
  csite => cpoly%site(isi)

  do ipa = 1,csite%npatches
     
     cpatch => csite%patch(ipa)

     do ico = 1,cpatch%ncohorts

        ! Update basal area, agb
        if(census_flag == 0 .or. cpatch%first_census(ico) == 1)then
           bdbh = min( int(cpatch%dbh(ico) * 0.1), 10) + 1
           ba = cpatch%nplant(ico) * cpatch%dbh(ico)**2
           cpoly%basal_area(cpatch%pft(ico), bdbh,isi) = cpoly%basal_area(cpatch%pft(ico), bdbh,isi) &
                +  csite%area(ipa) * ba * pi1 * 0.25
           cpoly%agb(cpatch%pft(ico), bdbh,isi) = cpoly%agb(cpatch%pft(ico), bdbh,isi) +  &
                ed_biomass(cpatch%bdead(ico), cpatch%balive(ico), cpatch%bleaf(ico), &
                cpatch%pft(ico), cpatch%hite(ico),cpatch%bstorage(ico)) &
                * cpatch%nplant(ico) * 10.0 * csite%area(ipa)
        endif
        
     end do
     
  end do
  
  return
end subroutine update_site_derived_props_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine ed_grndvap(nlev_sfcwater, nts, soil_water, soil_energy,    &
     sfcwater_energy, rhos, can_shv, ground_shv, surface_ssh)

  use soil_coms,   only: ed_nstyp, soil
  use grid_coms,   only: nzg
  use consts_coms, only:  pi1, grav, rvap
  use therm_lib  , only: rhovsil, qtk, qwtk8

  implicit none

  integer, intent(in) :: nlev_sfcwater ! # active levels of surface water
  integer, intent(in) :: nts           ! soil textural class (local name)
  
  real(kind=8), intent(in)  :: soil_water      ! soil water content [vol_water/vol_tot]
  real, intent(in)  :: soil_energy     ! [J/m^3]
  real, intent(in)  :: sfcwater_energy ! [J/kg]
  real, intent(in)  :: rhos            ! air density [kg/m^3]
  real, intent(in)  :: can_shv         ! canopy vapor spec hum [kg_vap/kg_air]
  real, intent(out) :: ground_shv      ! ground equilibrium spec hum [kg_vap/kg_air]
  real, intent(out) :: surface_ssh     ! surface (saturation) spec hum [kg_vap/kg_air]


  real, parameter :: gorvap = grav / rvap  ! gravity divided by vapor gas constant


  ! Local variables

  real :: slpotvn ! soil water potential [m]
  real :: alpha   ! "alpha" term in Lee and Pielke (1993)
  real :: beta    ! "beta" term in Lee and Pielke (1993)
  real :: tempk   ! surface water temp [K]
  real :: fracliq ! fraction of surface water in liquid phase

  ! surface_ssh is the saturation mixing ratio of the top soil or snow surface
  ! and is used for dew formation and snow evaporation.

  if (nlev_sfcwater > 0) then
     call qtk(sfcwater_energy,tempk,fracliq)
     surface_ssh = rhovsil(tempk) / rhos
  else
     
     ! Without snowcover, ground_shv is the effective saturation mixing
     ! ratio of soil and is used for soil evaporation.  First, compute the
     ! "alpha" term or soil "relative humidity" and the "beta" term.
     
     call qwtk8(soil_energy,soil_water*1.d3,soil(nts)%slcpd,tempk,fracliq)
     surface_ssh = rhovsil(tempk) / rhos
     
     slpotvn = soil(nts)%slpots * (soil(nts)%slmsts / soil_water) ** soil(nts)%slbs
     alpha = exp(gorvap * slpotvn / tempk)
     beta = .25 * (1. - cos (min(1.,soil_water / soil(nts)%sfldcap) * pi1)) ** 2
     ground_shv = surface_ssh * alpha * beta + (1. - beta) * can_shv
     
  endif

  return
end subroutine ed_grndvap
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine read_soil_moist_temp_ar(cgrid)

  use ed_state_vars, only: edtype, polygontype, sitetype, patchtype
  use soil_coms, only: soilstate_db, soil,slz
  use consts_coms, only: alli1000, cliq1000, cice1000, t3ple
  use grid_coms, only: nzg, ngrids
  
  implicit none

  type(edtype)      , target  :: cgrid  ! Alias for current ED grid
  type(polygontype) , pointer :: cpoly  ! Alias for current polygon
  type(sitetype)    , pointer :: csite  ! Alias for current site
  type(patchtype)   , pointer :: cpatch ! Alias for current patch
  integer :: ntext
  integer :: iland
  real :: xe_mean
  real :: ye_mean
  real :: ze_mean
  real :: glat
  real :: glon
  integer :: ipt
  integer :: ilat
  integer :: ilon
  integer :: ilatf
  integer :: ilonf
  integer :: nls
  integer :: nlsw1
  real :: soil_tempaux
  integer :: k
  real :: tmp1
  real :: tmp2
  real :: soilw1
  real :: soilw2
  logical :: l1
  integer :: ip
  integer :: ipy, isi, ipa, ico !Counters for all structures
  integer, parameter :: harvard_override = 0

! Putting these numbers as parameters, but we should think in a way to provide this info so we can make it more general.
  integer , parameter :: nlon=144, nlat=73
  real    , parameter :: dlon=2.5, dlat=2.5
  

  inquire(file=trim(soilstate_db),exist=l1)
  if(.not.l1)then
     print*,'You have ISOILSTATEINIT = 1, which means you read initial'
     print*,'soil moisture and temperature from file.  However, the '
     print*,'file you specified for SOILSTATE_DB,'
     print*
     print*,trim(soilstate_db)
     print*
     print*,'does not exist.'
     stop
  endif

  open(unit=12,file=trim(soilstate_db),form='formatted',status='old',position='rewind')
  latloop: do ilatf = 1,nlat  ! Reanalysis has 73 latitude points
     ! 1 corresponds to 90N
     lonloop: do ilonf = 1,nlon  ! Reanalysis has 144 longitude points
        ! 1 corresponds to 0E
        ! Read in reanalysis: two temperatures and moistures, corresponding to different depths
        read(unit=12,fmt=*)tmp1,tmp2,soilw1,soilw2
        ! soilw1, soilw2 are relative porosities and thus range from [0-1]
        ! tmp1, tmp2 are temperature in kelvin.

        ! Make sure it is not buggy
        if(tmp1 > 0.0 .and. tmp2 > 0.0 .and.   &
             soilw1 > 0.0 .and. soilw2 > 0.0)then

           ! Loop over land points
           polyloop: do ipy=1,cgrid%npolygons
              cpoly => cgrid%polygon(ipy)
              
              ! Land point lat, lon
              glat = cgrid%lat(ipy)
              glon = cgrid%lon(ipy)
              
              if(glon < 0.0) glon = glon + 360.0
              
              ! Find reanalysis point corresponding to this land point
              if(glat >= 0.0)then
                 ilat = nint((90.0 - glat)/dlat) + 1
              else
                 ilat = nlat - nint((90.0 - abs(glat))/dlat)
              endif
              ilon = int(glon/dlon) + 1
              
              ! If we are at the right point, fill the array
              if(ilat == ilatf .and. ilon == ilonf)then

                 ! Loop over sites and patches
                 siteloop: do isi=1,cpoly%nsites
                    csite => cpoly%site(isi)
                    
                    patchloop: do ipa=1,csite%npatches
                       cpatch => csite%patch(ipa)

                       do k=1,nzg
                          ntext = csite%ntext_soil(k,ipa)

                          if(abs(slz(k)) < 0.1)then
                             csite%soil_tempk(k,ipa) = tmp1
                             soil_tempaux = tmp1 - t3ple
                             csite%soil_water(k,ipa) = max(soil(ntext)%soilcp,   &
                                  soilw1 * soil(ntext)%slmsts)
                          else
                             csite%soil_tempk(k,ipa) = tmp2
                             soil_tempaux = tmp2 - t3ple
                             csite%soil_water(k,ipa) = max(soil(ntext)%soilcp,   &
                                  soilw2 * soil(ntext)%slmsts)
                          endif
                          if(soil_tempaux > 0.0)then
                             csite%soil_energy(k,ipa) = soil_tempaux * (soil(ntext)%slcpd   &
                                  + csite%soil_water(k,ipa) * cliq1000) +   &
                                    csite%soil_water(k,ipa) * alli1000
                             csite%soil_fracliq(k,ipa) = 1.0
                          else
                             csite%soil_energy(k,ipa) = soil_tempaux * (soil(ntext)%slcpd   &
                                  + csite%soil_water(k,ipa) * cice1000)
                             csite%soil_fracliq(k,ipa) = 0.0
                          end if
                       end do
                       if(harvard_override == 1)then
                          csite%soil_tempk(1,ipa)     = 277.6
                          csite%soil_tempk(2:4,ipa)   = 276.0
                          csite%soil_energy(1,ipa)    =   1.5293664e8
                          csite%soil_energy(2,ipa)    =   1.4789957e8
                          csite%soil_energy(3:4,ipa)  =   1.4772002e8
                          csite%soil_water(1:4,ipa)   =   0.41595
                          csite%soil_fracliq(1:4,ipa) =   1.0
                       endif
                       
                       nls   = csite%nlev_sfcwater(ipa)
                       nlsw1 = max(nls,1)
                 
                       call ed_grndvap(nls,                                &
                            csite%ntext_soil       (nzg,ipa),  &
                            csite%soil_water       (nzg,ipa),  &
                            csite%soil_energy      (nzg,ipa),  &
                            csite%sfcwater_energy(nlsw1,ipa),  &
                            cpoly%met(isi)%rhos,  &
                            csite%can_shv(ipa),  &
                            csite%ground_shv(ipa),  &
                            csite%surface_ssh(ipa))

                    end do patchloop
                 end do siteloop
                 
              end if
           end do polyloop
        end if
     end do lonloop
  end do latloop
  close(unit=12,status='keep')

  return
end subroutine read_soil_moist_temp_ar
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine update_polygon_derived_props_ar(cgrid)

  use ed_state_vars,only : edtype,polygontype

  implicit none
  integer :: ipy,isi
  type(edtype), target :: cgrid
  type(polygontype), pointer :: cpoly
  

  do ipy=1,cgrid%npolygons
     
     cgrid%agb(:,:,ipy) = 0.0
     cgrid%basal_area(:,:,ipy) = 0.0
     
     cpoly => cgrid%polygon(ipy)
     
     do isi = 1,cpoly%nsites
        cgrid%agb(:,:,ipy) = cgrid%agb(:,:,ipy) + cpoly%area(isi) * cpoly%agb(:,:,isi)
        cgrid%basal_area(:,:,ipy) = cgrid%basal_area(:,:,ipy) + cpoly%area(isi) * cpoly%basal_area(:,:,isi)
     enddo

  enddo

  return
end subroutine update_polygon_derived_props_ar
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!    This subroutine simply assigns the initial value for internal energy. The only reason !
! to do it separatedly is that we first load atmospheric-based variables, then we assign   !
! LAI and height. This should be called just at the initialization, during the run energy  !
! is what defines the temperature, not the other way.                                      !
!------------------------------------------------------------------------------------------!
subroutine initialize_vegetation_energy(cgrid)
   use ed_state_vars, only: edtype,polygontype,sitetype,patchtype
   use canopy_air_coms, only: hcapveg_ref, heathite_min
   use ed_therm_lib,only:calc_hcapveg
   use consts_coms, only: t3ple
   implicit none 
   !----- Argument ------------------------------------------------------------------------!
   type(edtype), target :: cgrid
   !----- Local variables -----------------------------------------------------------------!
   integer :: ipy,isi,ipa,ico
   type(polygontype), pointer :: cpoly
   type(sitetype)   , pointer :: csite
   type(patchtype)  , pointer :: cpatch
   real                       :: hcapveg
   !---------------------------------------------------------------------------------------!

   do ipy=1,cgrid%npolygons
      cpoly => cgrid%polygon(ipy)
      do isi=1,cpoly%nsites
         csite => cpoly%site(isi)
         do ipa=1,csite%npatches
            cpatch => csite%patch(ipa)
            do ico=1,cpatch%ncohorts
               
               hcapveg = calc_hcapveg(cpatch%bleaf(ico),cpatch%bdead(ico), &
                    cpatch%nplant(ico),cpatch%pft(ico))

               cpatch%veg_energy(ico) = hcapveg * (cpatch%veg_temp(ico)-t3ple)
            end do
         end do
      end do
   end do

   return
end subroutine initialize_vegetation_energy
!==========================================================================================!
!==========================================================================================!
