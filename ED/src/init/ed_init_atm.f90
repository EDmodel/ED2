!==========================================================================================!
!==========================================================================================!
subroutine ed_init_atm
  
  use ed_misc_coms,     only: runtype
  use ed_state_vars, only: edtype,polygontype,sitetype,patchtype,edgrid_g
  use soil_coms,     only: soil_rough, isoilstateinit, soil, slmstr, stgoff
  use consts_coms,    only: cliqvlme, cicevlme, t3ple, tsupercool, p00i, rocp,t00
  use grid_coms,      only: nzs, nzg, ngrids
  use fuse_fiss_utils, only: fuse_patches,fuse_cohorts
  use ed_node_coms, only: nnodetot,mynum,sendnum,recvnum
  use pft_coms,only : sla
  use ed_therm_lib,only : calc_hcapveg,ed_grndvap
  use therm_lib, only : thetaeiv,idealdenssh,reducedpress
  
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
  real    :: rvaux
  real    :: site_area_i, poly_area_i
  real    :: poly_lai, poly_nplant
  real    :: surface_temp, surface_fliq
  integer, parameter :: harvard_override = 0
  include 'mpif.h'
  integer :: ping,ierr
  integer :: npatches

  ping = 6 ! Just any number for MPI Send/Recv

  ! This subroutine fills the ED2 fields which depend on current 
  ! atmospheric conditions.

  do igr = 1,ngrids
     
     cgrid => edgrid_g(igr)
     
     ! First we need to update the meteorological fields.

     call update_met_drivers(cgrid)

     ! If this is a standard ED2 restart, we will read these fields in from 
     ! a history file and therefore not worry about setting them here.

     if(trim(runtype) /= 'HISTORY' ) then

     ! Loop over polygons, sites and patches
     
     do ipy = 1,cgrid%npolygons
        
        cpoly => cgrid%polygon(ipy)
        
        do isi = 1,cpoly%nsites
           
           csite => cpoly%site(isi)

           do ipa = 1,csite%npatches

              cpatch => csite%patch(ipa)

              !----------------------------------------------------------------------------!
              !      This first call is just to have the vegetation height so we can       !
              ! compute the initial canopy pressure...  It must be called again to have    !
              ! the storage right.                                                         !
              !----------------------------------------------------------------------------!
              call update_patch_derived_props(csite,cpoly%lsl(isi),cpoly%met(isi)%prss,ipa)

              csite%can_theta(ipa)    = cpoly%met(isi)%atm_theta
              csite%can_shv(ipa)      = cpoly%met(isi)%atm_shv
              csite%can_co2(ipa)      = cpoly%met(isi)%atm_co2
              csite%can_prss(ipa)     = reducedpress(cpoly%met(isi)%prss                   &
                                                    ,cpoly%met(isi)%atm_theta              &
                                                    ,cpoly%met(isi)%atm_shv                &
                                                    ,cpoly%met(isi)%geoht                  &
                                                    ,csite%can_theta(ipa)                  &
                                                    ,csite%can_shv(ipa)                    &
                                                    ,csite%can_depth(ipa))

              if (csite%can_theta(ipa) < 180.   .or. csite%can_theta(ipa) > 400. .or.      &
                  csite%can_shv(ipa)   < 1.e-8  .or. csite%can_shv(ipa) > 0.04   .or.      &
                  csite%can_prss(ipa)  < 40000. .or. csite%can_prss(ipa) > 110000.) then
                  write (unit=*,fmt='(a)') '======== Weird initial properties... ========'
                  write (unit=*,fmt='(a,f8.3)') ' LON       [ deg] = ',cgrid%lon(ipy)
                  write (unit=*,fmt='(a,f8.3)') ' LAT       [ deg] = ',cgrid%lat(ipy)
                  write (unit=*,fmt='(a,f8.3)')                                            &
                                     ' ATM_PRSS  [ hPa] = ',cpoly%met(isi)%prss  * 0.01
                  write (unit=*,fmt='(a,f8.3)')                                            &
                                     ' ATM_THETA [degC] = ',cpoly%met(isi)%atm_theta - t00
                  write (unit=*,fmt='(a,f8.3)')                                            &
                                     ' ATM_SHV   [g/kg] = ',cpoly%met(isi)%atm_shv * 1.e3
                  write (unit=*,fmt='(a,f8.3)')                                            &
                                     ' CAN_PRSS  [ hPa] = ',csite%can_prss(ipa)  * 0.01
                  write (unit=*,fmt='(a,f8.3)')                                            &
                                     ' CAN_THETA [degC] = ',csite%can_theta(ipa) - t00
                  write (unit=*,fmt='(a,f8.3)')                                            &
                                     ' CAN_SHV   [g/kg] = ',csite%can_shv(ipa)   * 1.e3
                  call fatal_error('Non-sense initial values!!!'                           &
                                  ,'ed_init_atm','ed_init_atm.f90')
              end if



              csite%can_temp(ipa)  = csite%can_theta(ipa)                                  &
                                   * (p00i *csite%can_prss(ipa)) ** rocp

              rvaux                = csite%can_shv(ipa) / (1. - csite%can_shv(ipa))
              csite%can_theiv(ipa) = thetaeiv(csite%can_theta(ipa),csite%can_prss(ipa)     &
                                             ,csite%can_temp(ipa),rvaux,rvaux,0)

              csite%can_rhos(ipa)  = idealdenssh(csite%can_prss(ipa)                       &
                                                ,csite%can_temp(ipa),csite%can_shv(ipa))

              ! Initialize stars
              csite%tstar (ipa) = 0.
              csite%ustar (ipa) = 0.
              csite%qstar (ipa) = 0.
              csite%cstar (ipa) = 0.
              
              csite%zeta  (ipa) = 0.
              csite%ribulk(ipa) = 0.
              
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

                 cpatch%veg_temp(ico)   = csite%can_temp(ipa)
                 cpatch%veg_water(ico)  = 0.0
                 cpatch%veg_fliq(ico)   = 0.0
                 cpatch%hcapveg(ico)    = calc_hcapveg(cpatch%bleaf(ico),cpatch%bdead(ico)   &
                                                      ,cpatch%balive(ico),cpatch%nplant(ico) &
                                                      ,cpatch%hite(ico),cpatch%pft(ico)      &
                                                      ,cpatch%phenology_status(ico))
                 cpatch%veg_energy(ico) = cpatch%hcapveg(ico)*cpatch%veg_temp(ico)
                 csite%hcapveg(ipa) = csite%hcapveg(ipa) + cpatch%hcapveg(ico)
              end do
           end do
        end do
     end do
     ! Initialize remaining soil properties.
     if(isoilstateinit == 1)then
        ! Initialize soil moisture, temperature, etc. from file specified in 
        ! the ED_NL.
        if (nnodetot /= 1) call MPI_Barrier(MPI_COMM_WORLD,ierr)
        if (mynum    /= 1) call MPI_Recv(ping,1,MPI_INTEGER,recvnum,110,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

        call read_soil_moist_temp(cgrid,igr)

        if (mynum     < nnodetot) call MPI_Send(ping,1,MPI_INTEGER,sendnum,110,MPI_COMM_WORLD,ierr)
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
                 csite%soil_tempk(1:nzg,ipa) = csite%can_temp(ipa) + stgoff(1:nzg)
                 do k = 1, nzg                 
                    if(csite%soil_tempk(k,ipa) > t3ple)then
                       nsoil=csite%ntext_soil(k,ipa)
                       csite%soil_fracliq(k,ipa) = 1.0
                       csite%soil_water(k,ipa)   = max(soil(nsoil)%soilcp                  &
                                                      , slmstr(k) * ( soil(nsoil)%slmsts   &
                                                                    - soil(nsoil)%soilwp)  &
                                                      + soil(nsoil)%soilwp)
                       csite%soil_energy(k,ipa)  = soil(nsoil)%slcpd                       &
                                                 * csite%soil_tempk(k,ipa)                 &
                                                 + csite%soil_water(k,ipa)  * cliqvlme     &
                                                 * (csite%soil_tempk(k,ipa) - tsupercool)
                    else
                       nsoil=csite%ntext_soil(k,ipa)
                       csite%soil_fracliq(k,ipa) = 0.0
                       csite%soil_water(k,ipa)   = max(soil(nsoil)%soilcp                  &
                                                      , slmstr(k) * ( soil(nsoil)%slmsts   &
                                                                    - soil(nsoil)%soilwp)  &
                                                      + soil(nsoil)%soilwp)
                       csite%soil_energy(k,ipa) = soil(nsoil)%slcpd                        &
                                                * csite%soil_tempk(k,ipa)                  &
                                                + csite%soil_water(k,ipa)                  &
                                                * cicevlme * csite%soil_tempk(k,ipa)
                    end if
                 end do

                 !----- Initial condition is with no snow/pond. ---------------------------!
                 csite%nlev_sfcwater(ipa)    = 0
                 csite%total_snow_depth(ipa) = 0
                 do k=1,nzs
                    csite%sfcwater_energy (k,ipa) = 0.
                    csite%sfcwater_depth  (k,ipa) = 0.
                    csite%sfcwater_mass   (k,ipa) = 0.
                    csite%sfcwater_tempk  (k,ipa) = csite%soil_tempk(nzg,ipa)
                    csite%sfcwater_fracliq(k,ipa) = csite%soil_fracliq(nzg,ipa)
                 end do

                 !----- Compute patch-level LAI, vegetation height, and roughness. --------!
                 call update_patch_derived_props(csite,cpoly%lsl(isi),cpoly%met(isi)%prss  &
                                                ,ipa)

                 !----- Compute the ground specific humidity. -----------------------------!
                 nsoil = csite%ntext_soil(k,ipa)
                 nls   = csite%nlev_sfcwater(ipa)
                 nlsw1 = max(1,nls)
                 call ed_grndvap(nls,nsoil,csite%soil_water(nzg,ipa)                       &
                                ,csite%soil_tempk(nzg,ipa),csite%soil_fracliq(nzg,ipa)     &
                                ,csite%sfcwater_tempk(nlsw1,ipa)                           &
                                ,csite%sfcwater_fracliq(nlsw1,ipa),csite%can_prss(ipa)     &
                                ,csite%can_shv(ipa),csite%ground_shv(ipa)                  &
                                ,csite%surface_ssh(ipa),surface_temp,surface_fliq)

              else
                 !----- Compute patch-level LAI, vegetation height, and roughness. --------!
                 call update_patch_derived_props(csite,cpoly%lsl(isi),cpoly%met(isi)%prss  &
                                                ,ipa)
              end if

              !----- Computing the storage terms for CO2, energy, and water budgets. ------!
              call update_budget(csite,cpoly%lsl(isi),ipa,ipa)

           end do
           
           ! Compute basal area and AGB profiles.
           call update_site_derived_props(cpoly, 0, isi)
           
        end do
        
        
        
     enddo
     
     call update_polygon_derived_props(cgrid)


     call fuse_patches(cgrid, igr)

     do ipy = 1,cgrid%npolygons
        
        ncohorts = 0
        npatches = 0
        poly_lai = 0.0
        poly_nplant = 0.0

        cpoly => cgrid%polygon(ipy)
        poly_area_i = 1./sum(cpoly%area(:))

        do isi = 1,cpoly%nsites
           
           csite => cpoly%site(isi)
           site_area_i = 1./sum(csite%area(:))
           
           do ipa = 1,csite%npatches
              npatches = npatches + 1
              cpatch => csite%patch(ipa)

              call fuse_cohorts(csite,ipa,cpoly%green_leaf_factor(:,isi),cpoly%lsl(isi))
              
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

  end if

end do

  return
end subroutine ed_init_atm
!==========================================================================================!
!==========================================================================================!
