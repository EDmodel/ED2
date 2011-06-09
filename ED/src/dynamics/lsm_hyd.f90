!! Land surface model LATERAL hydrology 
!! Subsurface hydrology based on TOPMODEL
!! Surface runoff hydrology based on site adjacency matrix

!!! IMPORTANT VARIABLE THROUGHT CODE:
!!
!!! TOPMODEL ::
!!  nzg = number of soil layers
!!  slz = array of soil layer depths
!!  cp%ntext_soil = array of site soil classes (see surfdata.f90)
!!  cp%watertable = soil water table depth (m)
!!  cp%soil_water = soil volumetric water content (m3/m3)
!!  cp%soil_energy  = soil energy (J/m3)
!!  soil(k)%slmsts = maximum volumetric soil moisture for soil type k (m2 water/m2 soil)
!!  soil(k)%soilcp = minimum volumetric soil moisture (m2/m2)
!!  soil(nsoil)%slcpd = dry soil heat capacity (units??)
!!  cs%moist_tau = characteristic redistribution timescale (seconds)
!!  cs%moist_zi = TOPMODEL equilibrium site water table depth
!!  cs%moist_W = site moisture index
!!  cpoly$zbar = polygon mean water table depth (m)

!! note: soil moisture defined volumetrically (m2/m2), surface water 
!! defined by height/mass (m/m2 = 1000 kg/m2; 1 mm precip/m2 = 1 liter = 1kg)
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! Initializes Lateral Hydrology
! defines initial state of hydrology parameters
! and redistributes soil moisture initial condition among sites
! to be in "equilibrium" and avoid long spinup
! inputs: polygon_list - list of ED2 polygons for all grid resolutions
! NOTE: skip if restarting from binary history dump
!==========================================================================================!
!==========================================================================================!
subroutine initHydrology()

  use hydrology_constants
  use grid_coms,only: ngrids
  use hydrology_coms, only: useTOPMODEL, useRUNOFF
  use ed_state_vars, only : edgrid_g
  use ed_node_coms, only: mynum
  implicit none
  
  integer :: igr,ipy
  integer :: mpolys,msites
  
  !! Initialize Subsurface lateral flow
  !! by doing a bulk redistribution of polygon water
  !! to it's initial "equilibrium" state across sites
  call initHydroSubsurface()
  mpolys=0
  msites=0
  do igr=1,ngrids
     mpolys=max(mpolys,edgrid_g(igr)%npolygons)
     do ipy=1,edgrid_g(igr)%npolygons
        msites=max(msites,edgrid_g(igr)%polygon(ipy)%nsites)
     end do
  end do
  write(unit=*,fmt='(4(a,1x,i9,1x))')  &
     'initHydrology | mynum=',mynum,'ngrids=',ngrids,'mpolys=',mpolys,'msites=',msites
  allocate(qw4out(ngrids,mpolys,msites))  !!just for debugging
  allocate(qh4out(ngrids,mpolys,msites))  !!just for debugging
  qw4out = 0.0
  qh4out = 0.0
  write(unit=*,fmt='(4(a,1x,i9,1x))')  &
     'Allocated | mynum=',mynum,'ngrids=',ngrids,'mpolys=',mpolys,'msites=',msites

  !!initialize runoff parameters
  do igr=1,ngrids
     call updateHydroParms (edgrid_g(igr))
  end do
  write(unit=*,fmt='(4(a,1x,i9,1x))')  &
     'Updated | mynum=',mynum,'ngrids=',ngrids,'mpolys=',mpolys,'msites=',msites
  deallocate(qw4out,qh4out)
  write(unit=*,fmt='(4(a,1x,i9,1x))')  &
     'Deallocated | mynum=',mynum,'ngrids=',ngrids,'mpolys=',mpolys,'msites=',msites
  return
end subroutine initHydrology
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! code for initializing site hydrology by shifting the soil water among sites 
! to initial equilibrium state
! inputs: polygon_list - list of ED2 polygons for all grid resolutions
! note: ** initialization routine does not conserve water or energy **
!==========================================================================================!
!==========================================================================================!
subroutine initHydroSubsurface()
  use hydrology_constants
  use hydrology_coms, only: useTOPMODEL, MoistRateTuning
  use soil_coms, only: soil,slz,dslz,dslzi
  use grid_coms, only: ngrids,nzg
  use ed_misc_coms, only: dtlsm
  use ed_state_vars, only: edgrid_g,edtype,polygontype,sitetype
  use consts_coms, only: wdns,cicevlme,tsupercool,t3ple,cliqvlme
 use therm_lib, only: qwtk
  implicit none

  type(edtype)      , pointer :: cgrid
  type(polygontype) , pointer :: cpoly
  type(sitetype)    , pointer :: csite
  integer :: igr, ipy, isi, ipa, k,  nsoil,slsl
  real :: zbar_site,area_poly,zmin
  !! SEE calcHydroSubsurface for variable definitions

  if(useTOPMODEL == 0) return

  !! Calculates TOPMODEL water table adjustment

  !! loop over grids and polygons
  do igr = 1,ngrids
     cgrid => edgrid_g(igr)
     !initialize variables

     do ipy=1,cgrid%npolygons
        cpoly => cgrid%polygon(ipy)
        cgrid%zbar(ipy) = 0.0
        !! FIRST calculate mean watertable (zbar)
        area_poly = 0.0        !proportion

        do isi=1,cpoly%nsites
           csite => cpoly%site(isi)
           zbar_site = 0.0  !reset site vars

           do ipa=1,csite%npatches
              !!first calculate water table depth for each patch 
              call calcWatertable(cpoly,isi,ipa)
              !!update site mean variables
              zbar_site   = zbar_site   + csite%watertable(ipa)*csite%area(ipa)
           enddo
           !update polygon mean variables
           area_poly       = area_poly       + cpoly%area(isi)
           cgrid%zbar(ipy) = cgrid%zbar(ipy) + zbar_site*cpoly%area(isi)
        end do
        !normalize by non-water site areas
        cgrid%zbar(ipy)      = cgrid%zbar(ipy)/area_poly 

        !!NEXT set new watertable depth
        do isi=1,cpoly%nsites

           csite => cpoly%site(isi)
           slsl=cpoly%lsl(isi)
           nsoil=cpoly%ntext_soil(slsl,isi)

           cpoly%moist_zi(isi) = cgrid%zbar(ipy) + (cpoly%moist_W(isi) - cgrid%wbar(ipy))/cpoly%moist_f(isi)/MoistRateTuning !TOPMODEL equilibrium water depth

           zmin = slz(slsl) + soil(nsoil)%soilcp/soil(nsoil)%slmsts*dslz(slsl) 
           !zmin = minimum height of water table assuming lowest soil layer is at minimal moisture level

           cpoly%moist_zi(isi) = max(cpoly%moist_zi(isi),zmin)
           
           do ipa=1,csite%npatches
              !!set soil moisture at specified height
              !!loop from bottom up 
              do k = cpoly%lsl(isi),nzg
                 nsoil = cpoly%ntext_soil(k,isi) !! mcd [9/30/08]
                 if(cpoly%moist_zi(isi) < slz(k+1)) then
                    csite%soil_water(k,ipa) = max(soil(nsoil)%soilcp, &
                         soil(nsoil)%slmsts*(cpoly%moist_zi(isi)-slz(k))*dslzi(k))
                    exit 
                 else
                    !! set soil to saturated
                    csite%soil_water(k,ipa) = soil(nsoil)%slmsts
                 endif
              end do
              !! reset energy
              do k = cpoly%lsl(isi),nzg
                 if(csite%soil_tempk(k,ipa) > t3ple)then
                    csite%soil_energy(k,ipa)  = soil(nsoil)%slcpd                       &
                                              * csite%soil_tempk(k,ipa)                 &
                                              + csite%soil_water(k,ipa)  * cliqvlme     &
                                              * (csite%soil_tempk(k,ipa) - tsupercool)
                 else
                    csite%soil_energy(k,ipa) = soil(nsoil)%slcpd                        &
                                             * csite%soil_tempk(k,ipa)                  &
                                             + csite%soil_water(k,ipa)                  &
                                             * cicevlme * csite%soil_tempk(k,ipa)
                 end if

              end do
              call calcWatertable(cpoly,isi,ipa)
           enddo
           cpoly%moist_tau(isi) = 100000.0*dtlsm !!set to some small rate initially
        enddo
     enddo
  enddo
end subroutine initHydroSubsurface
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine calcHydroSubsurface()
!! Main function for calculation subsurface lateral hydrologic flow
!! assumed to be called on the fast timestep dtlsm
  use ed_state_vars, only : edtype,polygontype,sitetype,edgrid_g
  use hydrology_constants
  use hydrology_coms, only: useTOPMODEL, FracLiqRunoff,MoistRateTuning
  use soil_coms, only: soil,slz,dslz
  use grid_coms, only: ngrids,nzg
  use ed_misc_coms, only: dtlsm
  use therm_lib, only: qtk,qwtk
  use consts_coms, only: wdns
  implicit none

  type(edtype)      , pointer :: cgrid                ! Alias for current ED grid-type
  type(polygontype) , pointer :: cpoly                ! Alias for current polygon
  type(sitetype)    , pointer :: csite                ! Alias for current site
  integer                     :: igr,ipy,isi,ipa      ! Grid,Polygon,Site and Patch Counters
  real                        :: zbar_site            ! site area-weighted mean water table depth (m)
  real                        :: area_land            ! polygon area that is land (proportion) 
  real                        :: area_subtract        ! polygon area that is has water removed (proportion) 
  real                        :: dz                   ! change in water table depth over time step (m)
  real                        :: dzadd                ! sum of net additions to the water table (m)
  real                        :: dw                   ! change in water content (m3/m2)
  real                        :: tempk                ! water temperature (K)
  real                        :: fracliq,fracliqtotal ! fraction of moisture that is liquid (proportion)
  real                        :: sheat                ! soil heat (units??)
  real                        :: swater               ! soil water (units??)
  real                        :: zmin                 ! minimum water table height (m)
  real                        :: soil_sat_energy      ! polygon-level mean energy in saturated zone (J/m2)
  real                        :: soil_sat_water       ! polygon-level mean water volume in saturated zone  (m3/m2)
  real                        :: soil_sat_heat        ! polygon-level mean soil heat capacity in saturated zone (J/m2)
  real                        :: energy_site          ! site-level mean energy in saturated zone (J/m2)
  real                        :: water_site           ! site-level mean water in saturated zone  (m/m2)
  real                        :: heat_site            ! site-level mean soil heat capacity in saturated zone (J/m2)
  real                        :: bf_site, bf_patch
  integer                     :: slsl,nsoil

  !!******************************************************************************!!
  !! If not using TOPMODEL, just do water table calculation then return           !!
  !!******************************************************************************!!
  if(useTOPMODEL == 0)then
     do igr = 1,ngrids
        cgrid => edgrid_g(igr)
        do ipy=1,cgrid%npolygons
           cpoly => cgrid%polygon(ipy)
           cgrid%zbar(ipy)     = 0.0
           cgrid%baseflow(ipy) = 0.0
           area_land = 0.0        !proportion
           do isi=1,cpoly%nsites
              csite => cpoly%site(isi)
              zbar_site   = 0.0  !reset site vars
              do ipa=1,csite%npatches
                 !!first calculate water table depth for each patch 
                 call calcWatertable(cpoly,isi,ipa)
                 !!update site mean variables
                 zbar_site = zbar_site + csite%watertable(ipa)*csite%area(ipa)
              end do
              !update polygon mean variables
              area_land       = area_land + cpoly%area(isi)
              cgrid%zbar(ipy) = cgrid%zbar(ipy)+ zbar_site * cpoly%area(isi)
           enddo
           !normalize variables by non-water site areas
           cgrid%zbar(ipy) = cgrid%zbar(ipy)/area_land 
        enddo
     enddo
     return
  endif

  ! Calculates TOPMODEL water table adjustment


  do igr = 1,ngrids
     cgrid => edgrid_g(igr)
     do ipy=1,cgrid%npolygons
        cpoly => cgrid%polygon(ipy)

        !!*************************************************************************!!
        !! Initialize variables                                                    !!
        !!*************************************************************************!!
        cgrid%zbar(ipy)     = 0.0 !! Polygon-average water table depth (m)
        cgrid%baseflow(ipy) = 0.0 !! Polygon total baseflow
        soil_sat_energy = 0.0  !mean energy in saturated zone (J/m2)
        soil_sat_water = 0.0   !mean water volume in saturated zone  (m3/m2)
        soil_sat_heat = 0.0    !mean soil heat capacity in saturated zone (J/m2)
        area_land = 0.0        !proportion of polygon that is land


        !!*************************************************************************!!
        !! Calculate patch-level watertable depth, polygon-level mean watertable   !!
        !! depth (ZBAR) and find out the amount of water, energy, and heat capacity!!
        !! is in the saturated zone                                                !! 
        !!*************************************************************************!!
        do isi=1,cpoly%nsites
           csite => cpoly%site(isi)

           zbar_site   = 0.0  !reset site vars
           energy_site = 0.0
           heat_site   = 0.0
           water_site  = 0.0

           do ipa=1,csite%npatches

              !!first calculate water table depth for each patch 
              call calcWatertable(cpoly,isi,ipa)

              !!update site-level mean variables
              zbar_site   = zbar_site   + csite%watertable(ipa)      * csite%area(ipa)
              energy_site = energy_site + csite%soil_sat_energy(ipa) * csite%area(ipa)
              heat_site   = heat_site   + csite%soil_sat_heat(ipa)   * csite%area(ipa)
              water_site  = water_site  + csite%soil_sat_water(ipa)  * csite%area(ipa)

              !!zero variables to be used soon
              csite%moist_dz(ipa) = 0.0
           end do
           !update polygon mean variables
           area_land       = area_land       +                  cpoly%area(isi)
           cgrid%zbar(ipy) = cgrid%zbar(ipy) + zbar_site      * cpoly%area(isi)
           soil_sat_energy = soil_sat_energy + energy_site    * cpoly%area(isi)
           soil_sat_heat   = soil_sat_heat   + heat_site      * cpoly%area(isi)
           soil_sat_water  = soil_sat_water  + water_site     * cpoly%area(isi)

        end do

        !!*************************************************************************!!
        !! normalize variables by non-water site areas                             !!
        !!*************************************************************************!!
        cgrid%zbar(ipy) = cgrid%zbar(ipy)/area_land 
        soil_sat_energy = soil_sat_energy/area_land
        soil_sat_heat   = soil_sat_heat/area_land
        soil_sat_water  = soil_sat_water/area_land

        !!*************************************************************************!!
        !! Compute mean temperature and liquid fraction of soil water              !!
        !! soil water converted from m3 -> kg                                      !!
        !!*************************************************************************!!
        call qwtk(soil_sat_energy,soil_sat_water*wdns,soil_sat_heat,tempk,fracliqtotal)


        !!*************************************************************************!!
        !!  If there is no saturated water or if saturated soil is more than       !!
        !! threshold frozen, skip over soil hydrology.                             !!
        !!*************************************************************************!!
        if(soil_sat_water > 1.e-6 .and. fracliqtotal > FracLiqRunoff) then 
           cgrid%swliq(ipy) = 1.0


           !!**********************************************************************!!
           !! Zero variables that track the amount of water removed from patches   !!
           !!**********************************************************************!!
           dzadd = 0.0            !total height to add to patches gaining water
           soil_sat_energy = 0.0  !total energy in removed water
           soil_sat_water = 0.0   !total removed water
           area_subtract = 0.0    !total area of sites where water is removed

           !!**********************************************************************!!
           !! First we will remove water from patches                              !!
           !!**********************************************************************!!
           do isi=1,cpoly%nsites

              csite => cpoly%site(isi)
              slsl=cpoly%lsl(isi)
              nsoil=cpoly%ntext_soil(slsl,isi)

              !!*******************************************************************!!
              !! Calculate new site-level equilibrium watertable depth (MOIST_ZI)  !!
              !!*******************************************************************!!
              cpoly%moist_zi(isi) = cgrid%zbar(ipy) + (cpoly%moist_W(isi) - &
                   cgrid%wbar(ipy))/cpoly%moist_f(isi)/MoistRateTuning !TOPMODEL water depth
              zmin = slz(slsl) + soil(nsoil)%soilcp/soil(nsoil)%slmsts*dslz(slsl)
              cpoly%moist_zi(isi)  = max(cpoly%moist_zi(isi),zmin)

              !!*******************************************************************!!
              !! Calculate rate constant for lateral redistribution time scale     !!
              !!*******************************************************************!!
              cpoly%moist_tau(isi) = soil(cpoly%ntext_soil(nzg-1,isi))%slmsts / & 
                   (MoistRateTuning*cpoly%moist_f(isi)*cgrid%Te(ipy)* &
                   exp(MoistRateTuning*cpoly%moist_f(isi)*min(0.0,cgrid%zbar(ipy))) & !characteristic redistribution timescale
                   *exp(-cgrid%wbar(ipy))*fracliqtotal) !!added a linear liquid fraction adjustment
              !! added a min(0,zbar) to get sensible behaviour when watertable perched -> flux is overland, not subsurface
              
              !!*******************************************************************!!
              !! Calculate rate constant for baseflow out of the bottom of the soil!!
              !!*******************************************************************!!
              cpoly%baseflow(isi) = soil(cpoly%ntext_soil(nzg-1,isi))%slmsts/ &
                   (MoistRateTuning*cpoly%moist_f(isi)*cpoly%moist_tau(isi))*1000.0 !kg/m2/s
              bf_site = 0.0
              

              do ipa=1,csite%npatches

                 !!****************************************************************!!
                 !! remove baseflow water and heat                                 !!
                 !!****************************************************************!!
                 bf_patch=cpoly%baseflow(isi)
                 call updateWatertableBaseflow(cpoly,isi,ipa,bf_patch)
                 bf_site = bf_site + bf_patch*csite%area(ipa)

                 
                 !!****************************************************************!!
                 !! calculate change in water table depth over time step           !!
                 !!****************************************************************!!
                 dz = (cpoly%moist_zi(isi)-csite%watertable(ipa))*dtlsm/cpoly%moist_tau(isi)
                 csite%moist_dz(ipa) = dz !! assign dz to patch

                 !!****************************************************************!!
                 !! if water flows out of the patch, attempt to reduce watertable  !!
                 !! by DZ (meters) and keep track of water and energy removed      !!
                 !!****************************************************************!!
                 if(dz < 0.0) then  
                    call updateWatertableSubtract(cpoly,isi,ipa,dz,sheat,swater)
                    soil_sat_water  = soil_sat_water  + swater  * csite%area(ipa) * cpoly%area(isi)
                    soil_sat_energy = soil_sat_energy + sheat   * csite%area(ipa) * cpoly%area(isi)
                    area_subtract   = area_subtract   +           csite%area(ipa) * cpoly%area(isi)
                 else                    
                    dzadd           = dzadd           + dz      * csite%area(ipa) * cpoly%area(isi)
                 end if

              end do

              !! Save the site-level mean baseflow
              cpoly%baseflow(isi) = bf_site

           end do

           !!**********************************************************************!!
           !! calculate mean heat content of water added to patches                !!
           !! (J/m3 water i.e. per volumetic water fraction)                       !!
           !!  assumes redistributed water is well mixed                           !!
           !!**********************************************************************!!
           sheat = soil_sat_energy / max(1.e-15,soil_sat_water)
           cgrid%sheat(ipy) = sheat
           if(soil_sat_water .lt. 1.e-15) soil_sat_water = 0.0

           !!**********************************************************************!!
           !! temperature and liquid fraction of water in flux                     !!
           !!**********************************************************************!!
           call qtk(sheat/1000.0,tempk,fracliq)

           !!**********************************************************************!!
           !! Next, add water to other patches                                     !!
           !!**********************************************************************!!
           do isi=1,cpoly%nsites
              csite => cpoly%site(isi)
              do ipa=1,csite%npatches
                 !change in water table depth
                 if(csite%moist_dz(ipa) > 0.0 .and. dzadd > tiny(1.0)) then

                    !!*************************************************************!!
                    !! calculate water (m3/m2) to add to each patch as a proportion!!
                    !! of the water that was removed from the upslope patches to   !!
                    !! ensure conservation of mass                                 !!
                    !!*************************************************************!!
                    dw = csite%moist_dz(ipa)/dzadd*soil_sat_water
                    if(dw /= dw) then
                       call fatal_error('NaN in water table','calcHydroSubsurface','lsm_hyd.f90')
                    endif
                    csite%moist_dz(ipa) = dw/0.435 !! approximation for output
                    
                    !!*************************************************************!!
                    !! Add water and energy to patch                               !!
                    !!*************************************************************!!
                    call updateWatertableAdd(cpoly,isi,ipa,dw,sheat)

                 end if
              end do

              !!*************************************************************!!
              !! increment polygon-level baseflow                            !!
              !!*************************************************************!! 
              cgrid%baseflow(ipy) = cgrid%baseflow(ipy) + cpoly%area(isi)*cpoly%baseflow(isi)

           end do
        else
           !!*************************************************************!!
           !! if soil water too frozen, set flag                          !!
           !!*************************************************************!!
           cgrid%swliq(ipy) = 0.0
        endif

        !!*************************************************************!!
        !! normalize fluxes for land area                              !!                                   
        !!*************************************************************!!
        cgrid%baseflow(ipy) = cgrid%baseflow(ipy)/area_land

     end do  !! loop over polygons
  end do     !! loop over grids
  return
end subroutine calcHydroSubsurface
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!MLO. Mike, I left this subroutine commented because it doesn't seem it's been called      !
!     anywhere in the code.                                                                !
!------------------------------------------------------------------------------------------!
!!!!!! subroutine calcHydroSubsurfaceNode(polygon_list)
!!!!!! !! Calculate fast-time step component of subsurface hydrology on the nodes
!!!!!! !! variables calculated with across-patch information are held constant between updates
!!!!!! !! and assumed to change slowly (equilibrium water depth, rate constant, heat content of shifted water)
!!!!!! !! inputs: polygon_list - list of ED2 polygons for all grid resolutions
!!!!!! !! note: used only when passing by site or patch
!!!!!! !!       if whole polygons are passed, call calcHydroSubsurface on each node
!!!!!!   use site_type
!!!!!!   use hydrology_constants
!!!!!!   use hydrology_coms, only: useTOPMODEL
!!!!!!   use soil_coms, only: soil
!!!!!!   use grid_coms, only: ngrids
!!!!!!   use ed_misc_coms, only: dtlsm
!!!!!!   implicit none
!!!!!!   include "mpif.h"
!!!!!! !  include 'rcommons.h'
!!!!!! !  include 'rpara.h'
!!!!!!   type(plist), dimension(ngrids) :: polygon_list
!!!!!!   type(polygon), pointer :: cpoly
!!!!!!   type(site), pointer :: cs
!!!!!!   type(patch), pointer :: cp
!!!!!!   integer :: ifm
!!!!!!   real :: sheat  ! soil heat in flux (dummy)
!!!!!!   real :: swater ! soil water in flux (dummy)
!!!!!!   real :: dz,dw     ! change in water table depth over time step (m)

!!!!!!   if(useTOPMODEL == 0) return

!!!!!!   !! Calculates TOPMODEL water table adjustment
!!!!!!   !! Assumes that ZI ("equilibrium" water table height) 
!!!!!!   !! tau (temporal rate), soil moisture heat content (sheat), 
!!!!!!   !! and whether the soil is frozen (swliq) stay unchanged between passes

!!!!!!   do ifm = 1,ngrids
!!!!!!      cpoly => polygon_list(ifm)%first_polygon
!!!!!!      do while(associated(cpoly))
!!!!!!         if(cpoly%swliq .eq.1) then
!!!!!!            cs => cpoly%first_site
!!!!!!            do while(associated(cs))
!!!!!!               cp => cs%oldest_patch
!!!!!!               do while(associated(cp))
!!!!!!                  !!first calculate water table depth for each patch 
!!!!!!                  call calcWatertable(cp)
!!!!!!                  dz = (cs%moist_zi-cp%watertable)*dtlsm/cs%moist_tau !RAMS rate correction
!!!!!!                  !change in water table depth
!!!!!!                  if(dz .gt.0.0) then
!!!!!!                     !!dw = cp%moist_dz/dzadd*soil_sat_water
!!!!!!                     dw = cp%moist_dz*0.45 !! quick fix that makes same assumptions as leaf2
!!!!!!                     call updateWatertableAdd(cp,dw,cpoly%sheat)
!!!!!!                  else
!!!!!!                     call updateWatertableSubtract(cp,dz,sheat,swater)
!!!!!!                  endif
!!!!!!                  cp => cp%younger
!!!!!!               enddo
!!!!!!               cs => cs%next_site
!!!!!!            enddo
!!!!!!         endif
!!!!!!         cpoly => cpoly%next_polygon
!!!!!!      enddo
!!!!!!   enddo
!!!!!! end subroutine calcHydroSubsurfaceNode
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine calcWatertable(cpoly,isi,ipa)
  ! Estimates  water table depth (meters)
  ! For stability, constrains the rate at water table estimates can change over time
  ! based on the variable Moist_dWT
  ! inputs: cpoly   - current polygon
  !         isi     - current site ID
  !         ipa     - current patch ID
  ! outputs: sets patch variables: 
  !         watertable      - depth of saturated zone (m)
  !         soil_sat_energy - energy in saturated zone 
  !         soil_sat_water  - water volume in saturated zone
  !         soil_sat_heat   - heat capacity of saturated zone
  use ed_state_vars, only: polygontype,sitetype
  use hydrology_constants
  use hydrology_coms, only: MoistSatThresh
  use soil_coms, only: soil,slz,dslz
  use grid_coms, only: nzg
  implicit none
  type(sitetype)   , pointer :: csite
  type(polygontype), target  :: cpoly
  integer, intent(in) :: isi,ipa
  
  integer :: k            ! soil layer counter
  integer :: nsoil        ! soil class
  real :: fracw           ! fraction of saturation (m/m)

  !!******************************************************************************!!
  !! Compute water table depth starting at the bottom of the lowest soil layer and!!
  !! summing over all saturated levels, defined here as levels that have more     !!
  !! than 95% of the full moisture capacity), but accounting for any deficit      !!
  !! in those levels.  Sum the total soil energy ( soil_energy* dz), water        !!
  !! content (soil_water * dz), and soil heat capacity (slcpd * dz ) in order     !!
  !! to compute a mean temperature of transported water.                          !!
  !!******************************************************************************!!

  csite => cpoly%site(isi)

  !WTold = csite%watertable(ipa)

  csite%watertable(ipa) = slz(cpoly%lsl(isi))   !initialize water table depth to deepest soil layer
  csite%soil_sat_energy(ipa) = 0.0
  csite%soil_sat_water(ipa)  = 0.0
  csite%soil_sat_heat(ipa)   = 0.0

  layerloop: do k = cpoly%lsl(isi),nzg
     nsoil = cpoly%ntext_soil(k,isi)  !look up soil type (switched to using SITE level soils [mcd 9/30/08]
     fracw = csite%soil_water(k,ipa) / soil(nsoil)%slmsts !calculate fraction of moisture capacity
     csite%watertable(ipa)      = csite%watertable(ipa)      + fracw *dslz(k)
     csite%soil_sat_energy(ipa) = csite%soil_sat_energy(ipa) + csite%soil_energy(k,ipa)*dslz(k)
     csite%soil_sat_water(ipa)  = csite%soil_sat_water(ipa)  + csite%soil_water(k,ipa)*dslz(k)
     csite%soil_sat_heat(ipa)   = csite%soil_sat_heat(ipa)   + soil(nsoil)%slcpd*dslz(k)
     if (fracw < MoistSatThresh) exit layerloop!change from original version to go one layer up
  end do layerloop
  !!store watertable depth integer
  csite%ksat(ipa) = real(max(cpoly%lsl(isi),k-1))

  if(k > nzg) then
     !!soil completely saturated, add surface water too
     if(csite%nlev_sfcwater(ipa) >= 1) then
        do k=1,csite%nlev_sfcwater(ipa)
           csite%soil_sat_energy(ipa) = csite%soil_sat_energy(ipa) + csite%sfcwater_energy(k,ipa) * &
                csite%sfcwater_mass(k,ipa)
        enddo
        k = csite%nlev_sfcwater(ipa)
        csite%soil_sat_water(ipa) = csite%soil_sat_water(ipa) + csite%sfcwater_depth(k,ipa)
        csite%watertable(ipa) = csite%watertable(ipa) + csite%sfcwater_depth(k,ipa)
     endif
  endif

   !! check that rate of change in water table depth is reasonable

   !! check that not first time and not surface water (allowed to change fast)
   !  if(WTold .lt. 0.0 .and. cpoly%moist_tau(isi) > 0.0) then 
        !! compute allowed rate of change
   !     dWT = abs(Moist_dWT*dtlsm/cpoly%moist_tau(isi))
   !     if( abs(WTold-csite%watertable(ipa)) > dWT) then
           !! need to limit rate of change to dz
   !        if(csite%watertable(ipa) > WTold) then
   !           csite%watertable(ipa) = WTold+dWT
   !        else
   !           csite%watertable(ipa) = WTold-dWT
   !        endif
   !     endif
   !  endif
   return
end subroutine calcWatertable
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine updateWatertableAdd(cpoly,isi,ipa,dw,sheat)
  !! Add subsurface water to a patch
  !! inputs: cpoly - polygon
  !!         isi   - site ID
  !!         ipa   - patch ID
  !!         dw    - change in water content (m3/m2)
  !!         sheat - heat content of added water (J/m3)
  use ed_state_vars, only: sitetype,polygontype
  use soil_coms, only: soil,slz,dslz,dslzi
  use grid_coms, only: nzg
  use therm_lib, only: qwtk
  use consts_coms, only: wdns
  implicit none
  type(polygontype) , target        :: cpoly
  integer           , intent(in)    :: isi,ipa
  real              , intent(inout) :: dw
  real              , intent(in)    :: sheat       !heat content of water to add (J/m3)

  type(sitetype)    , pointer       :: csite
  real                              :: tempk       !water temperature (K)
  real                              :: fracliq     !water fraction liquid (proportion)
  integer                           :: i,k         !counters
  integer                           :: nsoil       !soil type
  real                              :: wcap        !available capacity for adding moisture (m)
  real                              :: fracw       !proportion of saturation for a soil layer
  real                              :: dw_layer    !volumetric water content added to a layer (m3/m3)
  logical                           :: done

  !!******************************************************************************!!
  !! Return immediately if added water will be lost in round-off                  !!
  !!******************************************************************************!!
  if(dw .lt. tiny(1.0)) return

  done = .false.
  csite => cpoly%site(isi)

  !!******************************************************************************!!
  !! loop from lowest soil layer to the surface, adding moisture and energy       !!
  !!******************************************************************************!!
   layerloop: do k = cpoly%lsl(isi),nzg
      nsoil = cpoly%ntext_soil(k,isi)  !look up soil type
      fracw = csite%soil_water(k,ipa) / soil(nsoil)%slmsts  !calculate fraction of moisture capacity
      if(fracw < 1.0) then

         wcap = soil(nsoil)%slmsts-csite%soil_water(k,ipa)*dslz(k) !!m3/m2
         if(dw > wcap) then

            !!********************************************************************!!
            !! layer can be saturated, add as much as we can and subtract from dw !!                   
            !!********************************************************************!!
            dw_layer = wcap
            csite%soil_water(k,ipa) = soil(nsoil)%slmsts
            dw = dw - wcap

         else

            !!********************************************************************!!
            !! layer can accept all of the dw, add everything and exit the loop   !!                   
            !!********************************************************************!!
            dw_layer = dw
            csite%soil_water(k,ipa) = csite%soil_water(k,ipa)+dw*dslzi(k)
            dw = 0.0
            done = .true.

         endif

         !!************************************************************************!!
         !! update layer's soil heat based on the amount of water added            !!
         !!************************************************************************!!
         csite%soil_energy(k,ipa)  = csite%soil_energy(k,ipa) + dw_layer*sheat
         if(csite%soil_energy(k,ipa) /= csite%soil_energy(k,ipa)) then
            print*,"dw_layer",dw_layer
            print*,"sheat",sheat
            print*,"k= ",k," ipa= ",ipa
            print*,"wcap",wcap," dw", dw
            call fatal_error('Failed soil_energy sanity check in lsm_hyd' &
                 ,'updateWatertableAdd','lsm_hyd.f90')
         end if

         !!***********************************************************************!!
         !!  update soil temperature and liquid fraction                          !!
         !!***********************************************************************!!
         call qwtk(csite%soil_energy(k,ipa),csite%soil_water(k,ipa)*wdns            &
                   ,soil(nsoil)%slcpd,tempk,fracliq)
         csite%soil_tempk(k,ipa) = tempk
         csite%soil_fracliq(k,ipa) = fracliq

         if(done) exit layerloop
      end if
   end do layerloop

   if(.not. done) then
      !!***************************************************************************!!
      !! reached top layer and still more water add,  put in surface water         !!
      !! note: sfcwater_mass units = kg/m2                                         !!
      !!       sfcwater_energy = J/(kg H2O)                                        !!
      !!***************************************************************************!!

      if(csite%nlev_sfcwater(ipa) == 0 .or. csite%sfcwater_mass(1,ipa) < tiny(1.0)) then
         !!currently no surface water, create as liquid water
         csite%sfcwater_depth(1,ipa) = dw
         csite%sfcwater_mass(1,ipa) = dw*wdns
         csite%sfcwater_energy(1,ipa) = 0.001 * sheat
         csite%nlev_sfcwater(ipa) = 1
      else
         !!fill from bottom 
         !!assumes surface water is liquid
         !!need to work on case where may be snow or snow/liquid mix

         ! convert J/kg to J/m2
         csite%sfcwater_energy(1,ipa) = csite%sfcwater_energy(1,ipa) * csite%sfcwater_mass(1,ipa)
         ! still in J/m2
         csite%sfcwater_energy(1,ipa) = csite%sfcwater_energy(1,ipa) + sheat*dw
         ! get new mass
         csite%sfcwater_mass(1,ipa) = csite%sfcwater_mass(1,ipa) + dw*1000.0
         ! convert back to J/kg
         if(csite%sfcwater_mass(1,ipa) > 1.0e-10)then
            csite%sfcwater_energy(1,ipa) = csite%sfcwater_energy(1,ipa) / csite%sfcwater_mass(1,ipa)
         else
            csite%sfcwater_energy(1,ipa) = 0.0
         end if
         if(csite%sfcwater_energy(1,ipa) /= csite%sfcwater_energy(1,ipa))then
            call fatal_error('Failed sfcwater_energy sanity check in lsm_hyd' &
                 ,'updateWaterTableAdd','lsm_hyd.f90')
         end if
         do i= 1,csite%nlev_sfcwater(ipa)
            csite%sfcwater_depth(i,ipa) = csite%sfcwater_depth(i,ipa) + dw
         end do
      end if
   end if

   return
end subroutine updateWatertableAdd
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine updateWatertableSubtract(cpoly,isi,ipa,dz,sheat,swater)
!! Remove subsurface water from a patch
!! inputs: 
!!         cpoly  - current polygon
!!         isi    - current site ID
!!         ipa    - current patch ID
!!         dz     - change in water table height (meters, negative)
!! outputs:
!!         sheat  - heat content of removed water (J/m2?)
!!         swater - removed water (m3/m2)

   !!need to subtract water from patch
   use hydrology_constants
   use hydrology_coms, only: MoistSatThresh
   use ed_state_vars,  only: polygontype, sitetype 
   use consts_coms, only : cliqvlme,wdns,tsupercool
   use soil_coms, only: soil,slz,dslz,dslzi
   use grid_coms, only: nzg
   use therm_lib, only : qwtk

   implicit none
   type(polygontype), target      :: cpoly
   integer          , intent(in) :: ipa,isi
   real             , intent(in)   :: dz      ! change in water table height (always negative)
   real             , intent(out) :: sheat
   real             , intent(out) :: swater

   type(sitetype)   , pointer    :: csite

   integer :: k       ! counters
   integer :: nsoil   ! soil type
   real    :: dzl     ! local copy of dz, safe to change 
   real    :: wcap    ! available capacity for adding moisture (m)
   real    :: fracw   ! proportion of saturation for a soil layer
   real    :: dh      ! change in heat (J/m3)
   real    :: dw      ! change in water (m3/m3)
   real    :: fracliq ! water fraction liquid (proportion)
   real    :: tempk   ! water temperature (K)
   real    :: wtmp    ! water content of layer before adjustment
   logical :: done    ! status variable to end loop early

   if( dz > 0.0) then
      write(unit=*,fmt=*) 'dz > 0 in updateWatertableSubtract',dz
      call fatal_error('invalid dz!','updateWatertableSubtract','lsm_hyd.f90') 
   end if

   csite => cpoly%site(isi)

   done = .false.
   dzl = dz        !make local copy
   sheat = 0.0
   swater = 0.0

   !!start by finding top of saturated layer
   k = int(csite%ksat(ipa))

   nsoil = cpoly%ntext_soil(k,isi)  !look up soil type
   fracw = csite%soil_water(k,ipa) / soil(nsoil)%slmsts  !calculate fraction of moisture capacity
   if(fracw > MoistSatThresh) then  !if layer above threshold, move up one
      k = k+1
      nsoil = cpoly%ntext_soil(k,isi)
      fracw = csite%soil_water(k,ipa)
   end if

   if(k > nzg) then
      !! soil completely saturated, reset to top layer
      k = nzg
      nsoil = cpoly%ntext_soil(k,isi)  !look up soil type
      fracw = csite%soil_water(k,ipa) / soil(nsoil)%slmsts  !calculate fraction of moisture capacity
   end if

   !!work down from saturation layer, removing water & heat
   do while (k >= cpoly%lsl(isi))
      call qwtk(csite%soil_energy(k,ipa),csite%soil_water(k,ipa)*1.e3,soil(nsoil)%slcpd,tempk,fracliq)
      !capacity for layer to loose moisture (unit = meters)
      wcap = (fracw - soil(nsoil)%soilcp/soil(nsoil)%slmsts) * dslz(k)*fracliq 

      if (wcap > tiny(1.0)) then 
         wtmp = real(csite%soil_water(k,ipa))
         !!      if(wcap > dzl) then 
         if(-dzl > wcap) then 
            !!empty soil layer completely 
            dw = -1.0*(csite%soil_water(k,ipa) - soil(nsoil)%soilcp)
            dzl = dzl + wcap
         else
            !!reduce water by dz
            dw = dzl*soil(nsoil)%slmsts*dslzi(k)
            dzl = 0.0
            done = .true.
         end if
         
         !!update soil water
         csite%soil_water(k,ipa) = csite%soil_water(k,ipa) + dw
         !! error message for positive dw or soil below capacity
         if(dw .gt. 1.0E-9 .or. (csite%soil_water(k,ipa)+1.0e-8) < soil(nsoil)%soilcp) then 
            print*," "
            print*,"dw = ",dw," in updateWatertableSubtract"
            print*,dzl,wcap,k,csite%watertable(ipa),done,fracw,nsoil
            print*,wtmp,csite%soil_water(k,ipa),soil(nsoil)%soilcp,soil(nsoil)%slmsts
            call fatal_error ('Bad dw!','updateWatertableSubtract','lsm_hyd.f90')
         endif
         
         !!update soil heat
         dh = dw*cliqvlme*(tempk-tsupercool)
         csite%soil_energy(k,ipa) = csite%soil_energy(k,ipa) + dh
         sheat = sheat - dh*dslz(k)   !cumulative sum as return value
         swater = swater - dw*dslz(k)
         if(csite%soil_energy(k,ipa) /= csite%soil_energy(k,ipa)) then
            call fatal_error('Failed soil_energy sanity check in lsm_hyd' &
                 ,'updateWatertableSubtract','lsm_hyd.f90')
         end if
         
         
         !!update soil temperature
         call qwtk(csite%soil_energy(k,ipa),csite%soil_water(k,ipa)*wdns               &
              ,soil(nsoil)%slcpd,tempk,fracliq)
         csite%soil_tempk(k,ipa) = tempk
         csite%soil_fracliq(k,ipa) = fracliq
      end if

      !!iterate
      if(done .or. k == cpoly%lsl(isi)) then
         k = 0
      else
         k = k-1
         nsoil = cpoly%ntext_soil(k,isi)  !look up soil type
         fracw = csite%soil_water(k,ipa) / soil(nsoil)%slmsts
      endif
   enddo
   return
end subroutine updateWatertableSubtract
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine updateWatertableBaseflow(cpoly,isi,ipa,baseflow)
   use ed_state_vars, only: polygontype, sitetype
   use soil_coms, only: soil,slz,dslz,dslzi,slcons1
   use ed_misc_coms, only: dtlsm
   use consts_coms, only: cliqvlme, tsupercool
   use therm_lib, only : qwtk
   implicit none

   real, parameter :: freezeCoef = 7.0        !! should probably move to the com
   type(polygontype) , target        :: cpoly
   integer           , intent(in)    :: isi, ipa
   real              , intent(inout) :: baseflow !read in as mm/s
  
   type(sitetype)    , pointer       :: csite
   real                              :: bf   !baseflow in water content (meters/timestep)
   real                              :: potn_fd,wflux_fd ! potential & flux for free-drainage
   integer                           :: nsoil,slsl
   real                              :: tempk, fracliq, freezeCor


   csite => cpoly%site(isi)
   slsl  = cpoly%lsl(isi)

   !! determine freeze
   call qwtk(csite%soil_energy(slsl,ipa),csite%soil_water(slsl,ipa)*1.e3,soil(nsoil)%slcpd,tempk,fracliq)
   freezeCor = fracliq
   if(freezeCor .lt. 1.0) freezeCor = 10.0**(-freezeCoef*(1.0-freezeCor))

   !! calc max free-drainage as cap to baseflow
   !! assumes layer below is permenantly at minimal water capacity
   slsl  = cpoly%lsl(isi)
   nsoil = cpoly%ntext_soil(slsl,isi)
   potn_fd = -dslzi(slsl)+soil(nsoil)%slpots* &
        ((soil(nsoil)%slmsts/soil(nsoil)%soilcp)**soil(nsoil)%slbs - &
        (soil(nsoil)%slmsts/csite%soil_water(slsl,ipa))**soil(nsoil)%slbs)
   wflux_fd = slcons1(slsl,nsoil)  &
                 * (csite%soil_water(1,ipa)/ soil(nsoil)%slmsts)**(2. * soil(nsoil)%slbs + 3.)  &
                 * potn_fd * freezeCor
   bf = min(-wflux_fd,baseflow)*dtlsm*0.001

   if(bf < 0.0) then
      write (unit=*,fmt=*) "bf in wrong direction",bf
      call fatal_error('bad bf!','updateWatertableBaseflow','lsm_hyd.f90')
   end if

   !! first, remove baseflow water and heat
   csite%soil_water(slsl,ipa) = csite%soil_water(slsl,ipa)-bf*dslzi(slsl)
   if(csite%soil_water(slsl,ipa) < soil(nsoil)%soilcp) then
      write(unit=*,fmt=*) 'Bottom dry:',csite%soil_water(1,ipa),bf,soil(nsoil)%soilcp
      bf = bf + (csite%soil_water(slsl,ipa)-soil(nsoil)%soilcp)*dslz(slsl)
      csite%soil_water(slsl,ipa) = soil(nsoil)%soilcp
   end if
   csite%soil_energy(slsl,ipa) = csite%soil_energy(slsl,ipa)-bf*cliqvlme*(csite%soil_tempk(slsl,ipa)-tsupercool)

   if(csite%soil_energy(slsl,ipa) /= csite%soil_energy(slsl,ipa)) then
                    call fatal_error('Failed soil_energy sanity check in lsm_hyd' &
                                    ,'updateWatertableBaseflow','lsm_hyd.f90')
                 end if

   baseflow = bf*1000.0/dtlsm !! reassign for return (m/step->mm/sec)
   return
end subroutine updateWatertableBaseflow
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
! MLO. Mike, I added the node-number in the file name, just to avoid calling mpi in the    !
!      testing subroutine...
!==========================================================================================!
!==========================================================================================!
subroutine writeHydro()
  !! code for writing out data for lateral hydrology testing and validation
  !! not intended as a standard output, key hydologic variables should be migrated to other files
  !! inputs:
  !!        
  !! output:
  !!        writes to <output_filepath>-HYD<year>
  !!        after every HydroOutputPeriod calls
  use ed_state_vars       , only : edgrid_g,edtype, polygontype,sitetype
  use ed_misc_coms           , only : ffilout, current_time
  use hydrology_constants
  use hydrology_coms      , only : HydroOutputPeriod, useRUNOFF, useTOPMODEL
  use grid_coms            , only : nzg,ngrids
  use soil_coms           , only : dslzi
  use ed_node_coms        , only : mynum
  use consts_coms         , only : t00
  use therm_lib           , only : qtk

  implicit none
  type(edtype)      , pointer          :: cgrid                              !Alias for "Current ED grid structure"
  type(polygontype) , pointer          :: cpoly                              !Alias for "Current polygon"
  type(sitetype)    , pointer          :: csite                              !Alias for "Current site"

  character(len=512)                   :: fname                              !output file name
  character(len=512) , save            :: last_fname=' '                     !file name last time function was called
  integer            , save            :: count = 1                          !function calls since last output
  logical            , save            :: first = .true.                     !indicator for first time called
  integer                              :: igr, ipy, isi, ipa                 !Grid, Polygon, Site and Patch counters.
  integer                              :: k                                  !soil layer counter
  real                                 :: WTbar                              !site-mean water table depth (meters)
  real                                 :: DZbar                              !site-mean change in water table height (meters)
  real                                 :: SSE,SSW,SSH,ksat
  real, save, allocatable,dimension(:) :: WPbar                              !site-mean volumetric soil moisture by soil layer (m3/m3) 
  real, save, allocatable,dimension(:) :: WPHbar                             !site-mean soil heat by soil layer
  real, save, allocatable,dimension(:) :: WPTbar                             !site-mean soil temperature by soil layer
  real                                 :: tempk, fracliq, runoff_t,runoff_fl !subsurface moisture flux temperature (K) 
  real                                 :: area_land
  character(len=255)                   :: format_str
  real                                 :: tpw,tsw !! total site and patch water

  return  
!!! ALL variables in this section of code that are NOT for debugging should be migrated to the analysis files
!!! MCD

  if(useTOPMODEL == 0) return

  !! allocate dynamic array once
  if(first) then
     allocate(WPbar(nzg))
     allocate(WPHbar(nzg))
     allocate(WPTbar(nzg))
     first = .false.
  endif

  !!if time to output
  if(count == HydroOutputPeriod) then
     count = 1        ! reset counter

     gridloop: do igr=1,ngrids

        cgrid => edgrid_g(igr)
        !! open file 
        write(fname,'(a,a,i4.4,2(a,i2.2))')trim(ffilout),'-HYD',current_time%year,'-grid_',igr,'-node_',mynum  !generate file name
        if(.not.(trim(fname).eq.trim(last_fname))) then
           !open file new
           open(unit=12,file=trim(fname),form='formatted',status='replace')
           last_fname = trim(fname)
        else
           open(unit=12,position='APPEND',file=trim(fname),form='formatted',status='old')
        end if

        !iterate over polygons and sites
        polyloop: do ipy=1,cgrid%npolygons
           cpoly => cgrid%polygon(ipy)
           if(useRUNOFF == 0) cgrid%runoff(ipy) = 0.0
           area_land = 0.0

           !calc temperature and liquid fraction of water in flux
           call qtk(cgrid%sheat(ipy)*0.001,tempk,fracliq)

           siteloop: do isi=1,cpoly%nsites

               csite => cpoly%site(isi)

               WTbar  = 0.0
               DZbar  = 0.0
               WPbar  = 0.0
               WPHbar = 0.0
               WPTbar = 0.0
               SSE    = 0.0
               SSW    = 0.0
               SSH    = 0.0
               ksat   = 0.0
               runoff_t = t00
               tsw    = 0.0
               if(useRUNOFF == 0) cgrid%runoff(ipy) = cgrid%runoff(ipy) + cpoly%area(isi)*cpoly%runoff(isi)
               area_land = area_land + cpoly%area(isi)
               if(isi == cpoly%nsites .and. useRUNOFF == 0) cgrid%runoff(ipy) = cgrid%runoff(ipy)/area_land
               if(cpoly%runoff(isi) > 0.0) call qtk(cpoly%avg_runoff_heat(isi)/cpoly%runoff(isi),runoff_t,runoff_fl)

               do ipa=1,csite%npatches

                  tpw = 0.0
                  do k = cpoly%lsl(isi),nzg
                     tpw = tpw + csite%soil_water(k,ipa)*dslzi(k)
                  enddo
                  tsw = tsw + tpw*csite%area(ipa)

                  !! calculate site-level means
                  WTbar = WTbar + csite%watertable(ipa)               * csite%area(ipa)
                  DZbar = DZbar + csite%moist_dz(ipa)                 * csite%area(ipa)
                  SSE   = SSE   + csite%soil_sat_energy(ipa)          * csite%area(ipa)
                  SSW   = SSW   + csite%soil_sat_water(ipa)           * csite%area(ipa)
                  SSH   = SSH   + csite%soil_sat_heat(ipa)            * csite%area(ipa)
                  ksat  = ksat  + csite%ksat(ipa)                     * csite%area(ipa)
                  do k =cpoly%lsl(isi),nzg
                     WPbar(k)  = WPbar(k)  + csite%soil_water(k,ipa)  * csite%area(ipa)
                     WPTbar(k) = WPTbar(k) + csite%soil_tempk(k,ipa)  * csite%area(ipa)
                  end do
               end do
             
              write(format_str,fmt='(a,i2.2,a)') &
                '(2(i2.2,1x),es16.9,1x,3(i2,1x),',20+3*nzg,'(es16.9,1x))'          
              !! output site data
              write(unit=12,fmt=trim(format_str)) current_time%month,current_time%date, current_time%time, &
                   ipy,isi, &  !!index variables
                   cgrid%swliq(ipy),cgrid%sheat(ipy),tempk,cgrid%zbar(ipy),cpoly%moist_zi(isi), &
                   cpoly%moist_tau(isi),WTbar,DZbar, &
                   WPbar,WPTbar,SSE,SSW,SSH,ksat, & !! subsurface variables WPHbar
                   cpoly%runoff(isi),cpoly%avg_runoff_heat(isi),runoff_t, &
                   qw4out(igr,ipy,isi),qh4out(igr,ipy,isi), & !! surface variables
                   cpoly%baseflow(isi),cgrid%baseflow(ipy),cgrid%runoff(ipy), & !!discharge variables
                   tsw
           end do siteloop
        end do polyloop
        close(unit=12,status='keep')
     end do gridloop
  else
     !! if not output, increment counter
     count = count + 1
  endif
  return
end subroutine writeHydro
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine calcHydroSurface()
  !! Main function for calculating surface runoff
  !! Calculates overland flow (runoff rate kg/s/m2)
  !! Runoff direction is based on site adjacency matrix
  !! Runoff rate is based on Manning's equation and is a function of:
  !!    water depth, slope, surface roughness 
  !! surface roughness is calculated in updateHydroParms
  !! Inputs:
  !! Outputs:
  !!         sets cgrid%runoff, cgrid%runoff_heat, cpoly%runoff, cpoly%runoff_heat
  use ed_state_vars, only : edgrid_g,edtype,polygontype,sitetype
  use hydrology_constants
  use hydrology_coms, only: useRUNOFF,FracLiqRunoff,runoff_Vmax
  use grid_coms, only: ngrids, nzg
  use ed_misc_coms, only: dtlsm
  use grid_coms,only:ngrids
  use consts_coms, only : cliqvlme,cliq,t3ple,qicet3,qliqt3,tsupercool
  use soil_coms, only: water_stab_thresh
  use therm_lib, only: qtk
  implicit none
  type(edtype)       , pointer :: cgrid
  type(polygontype)  , pointer :: cpoly
  type(sitetype)     , pointer :: csite
  
  integer :: igr,ipy,isi,ines,ipa  ! Grid, polygon, site, neighbouring site and patch counters
  integer :: top_surf_water        ! index of top surface water layer
  integer :: curr_site             ! site counter
  integer :: i
 ! integer :: nsites            ! number of sites in the polygon
  real, allocatable :: qwout(:) ! flux of water out of a site (kg/s/m2?)
  real, allocatable :: qhout(:) ! flux of heat out of a site ( units /s/m2? ******************)
  real :: tempk                 ! water temperature (K)
  real :: fracliq               ! water fraction liquid (proportion)
  real :: surf_water_depth      ! depth of liquid surface water (m)
  real :: surf_water_heat       ! energy of liquid surface water (J/m3 ??)
  real :: swd_i                 ! surface water depth in layer i (m)
  real :: flow_vel,flow_denom   ! overland flow velocity (m/s)
  real :: patch_water_out       ! patch water flux out (kg/m2/s)
  real :: patch_heat_out        ! patch heat flux out (J/m2/s??? ***************)
  real :: site_water_in         ! site water flux in (kg/m2/s)
  real :: site_heat_in          ! site heat flux in (J/m2/s)


  !!! Mike,what does this eleven stand for?
  real,dimension(11) :: hts     ! heights of surface water layers (m)




  if(useRUNOFF == 0) return

  !! Loop over grids and polygons
  do igr = 1,ngrids
     cgrid => edgrid_g(igr)
     do ipy=1,cgrid%npolygons
        cpoly => cgrid%polygon(ipy)
        !! FIRST PASS, CALCULATE GROSS RUNOFF

        !! initialize variables
        allocate(qwout(cpoly%nsites))
        allocate(qhout(cpoly%nsites))
        qwout = 0.0
        qhout = 0.0

        do isi=1,cpoly%nsites
           csite => cpoly%site(isi)

           !site specific parms
           do ipa=1,csite%npatches
              !first check that there IS surface water
              top_surf_water = csite%nlev_sfcwater(ipa)
              if(top_surf_water > 0) then
                 !second check that some is liquid
                 call qtk(csite%sfcwater_energy(top_surf_water,ipa),tempk,fracliq)
                 if(fracliq > FracLiqRunoff) then
                    !! calculate water depth
                    swd_i = csite%sfcwater_mass(top_surf_water,ipa)*0.001*fracliq !convert liquid fraction from kg to meters
                    surf_water_depth = swd_i
                    surf_water_heat = swd_i*cliqvlme*(tempk-tsupercool)
                    do i=(top_surf_water-1),1,-1
                       call qtk(csite%sfcwater_energy(i,ipa),tempk,fracliq)
                       swd_i = csite%sfcwater_mass(top_surf_water,ipa)*0.001*fracliq
                       surf_water_depth = surf_water_depth + swd_i
                       surf_water_heat = surf_water_heat + swd_i*cliqvlme*(tempk-tsupercool)
                    end do
                    !!sanity check
                    if(surf_water_depth > 1.0) then  !!if there's more than 1 m of standing water
                      write (unit=*,fmt=*) 'FLOOD WARNING !!!!!!!!!!'
                      write(unit=*,fmt='(4(a,1x))') 'Level','Sfc.Water','SfcWater.Mass','Sfc.Water.Energy'
                      do i=1,top_surf_water
                         write(unit=*,fmt='(2(i9,1x),2(es14.7,1x))') &
                            i,csite%nlev_sfcwater(ipa),csite%sfcwater_mass(i,ipa),csite%sfcwater_energy(i,ipa)
                      end do
                      call fatal_error('I''m drowning!!!','calcHydroSurface','lsm_hydro.f90')
                    endif

                    !! calculate flow velocity (m/s)
                    flow_denom = 1.0+csite%runoff_a(2,ipa)*surf_water_depth**(4/3) & 
                         - csite%runoff_a(3,ipa)*surf_water_depth**(7/3)
                    if(flow_denom > 1.0) then
                       flow_vel = surf_water_depth**(2/3)*csite%runoff_a(1,ipa)/sqrt(flow_denom)
                    else
                       flow_vel = surf_water_depth**(2/3)*csite%runoff_a(1,ipa) !! asymptotic vel w/o vegetation
                    endif
                    flow_vel = min(flow_vel,runoff_vmax) !! clamp runoff velocity to maximum sensible value

                    !!patch level flux out (kg/m2/s)
                    patch_water_out = flow_vel*surf_water_depth*1000.0
                    patch_heat_out  = flow_vel*surf_water_heat  !energy/m2/s

                    !!update vars
                    csite%runoff(ipa) = patch_water_out
                    csite%avg_runoff_heat(ipa) = patch_heat_out
                    qwout(isi) = qwout(isi) + patch_water_out*csite%area(ipa)
                    qhout(isi) = qhout(isi) + patch_heat_out*csite%area(ipa)
                 endif
              endif
           end do
           !! for output
           qw4out(igr,ipy,isi) = qwout(isi)
           qh4out(igr,ipy,isi) = qhout(isi)
        end do
        


        !! SECOND PASS: CALCULATE RUN-IN
        do isi=1,cpoly%nsites
           !!calc NET site-level run-in rate (adjusts patch net water, prev calc was gross patch runoff)
           site_water_in = 0.0
           site_heat_in = 0.0

           do ines=1,isi-1
              site_water_in = site_water_in + qwout(ines)*cgrid%site_adjacency(ines,isi,ipy)
              site_heat_in  = site_heat_in  + qhout(ines)*cgrid%site_adjacency(ines,isi,ipy)
           end do
           do ines=isi+1,cpoly%nsites
              site_water_in = site_water_in + qwout(ines)*cgrid%site_adjacency(ines,isi,ipy)
              site_heat_in  = site_heat_in  + qhout(ines)*cgrid%site_adjacency(ines,isi,ipy)
           end do

           !!calc NET site-level run-off rate (just for output)
           cpoly%runoff(isi)          = 0.0
           cpoly%avg_runoff_heat(isi) = 0.0
           do ines=1,isi-1
                 cpoly%runoff(isi)          = cpoly%runoff(isi) + qwout(isi)*cgrid%site_adjacency(isi,ines,ipy)
                 cpoly%avg_runoff_heat(isi) = cpoly%avg_runoff_heat(isi) + qhout(isi)*cgrid%site_adjacency(isi,ines,ipy)
           end do
           do ines=isi+1,cpoly%nsites+1
                 cpoly%runoff(isi)          = cpoly%runoff(isi) + qwout(isi)*cgrid%site_adjacency(isi,ines,ipy)
                 cpoly%avg_runoff_heat(isi) = cpoly%avg_runoff_heat(isi) + qhout(isi)*cgrid%site_adjacency(isi,ines,ipy)
           end do

           !! calc patch layer run in
           do ipa=1,csite%npatches
              !! adjust runoff rates for run-on
              !! negative values indicate a net gain of water from overland flow
              csite%runoff(ipa)          = (1.-cgrid%site_adjacency(isi,isi,ipy))*csite%runoff(ipa)          - site_water_in 
              csite%avg_runoff_heat(ipa) = (1.-cgrid%site_adjacency(isi,isi,ipy))*csite%avg_runoff_heat(ipa) - site_heat_in
              ! changed runoff code to assume recycled runoff occurrs at a patch level, not a intra-site level
              ! less realistic, but testing whether can avoid putting intersite runoff in the integrator

              !!compute contribution of each surface water layer -> hts(i)/surf_water_depth
              surf_water_depth = 0.0
              top_surf_water = csite%nlev_sfcwater(ipa)
              do i=1,top_surf_water
                 call qtk(csite%sfcwater_energy(i,ipa),tempk,fracliq)
                 hts(i) = csite%sfcwater_mass(i,ipa)*0.001*fracliq
                 surf_water_depth = surf_water_depth + hts(i)
              end do

              !! calc actual plot-level runoff/run-on
              do i=1,top_surf_water
                 ! Convert to J/m2
                 csite%sfcwater_energy(i,ipa) = csite%sfcwater_energy(i,ipa) * csite%sfcwater_mass(i,ipa)
                 csite%sfcwater_mass(i,ipa) = csite%sfcwater_mass(i,ipa) - csite%runoff(ipa)*dtlsm*hts(i)/surf_water_depth
                 ! Subtract and convert back to J/kg
                 if(csite%sfcwater_mass(i,ipa) > 1.0e-10)then
                    csite%sfcwater_energy(i,ipa) = (csite%sfcwater_energy(i,ipa) -   &
                         csite%avg_runoff_heat(ipa)*dtlsm*hts(i)/surf_water_depth) /  &
                         csite%sfcwater_mass(i,ipa)
                 else
                    csite%sfcwater_energy(i,ipa) = 0.0
                 end if

                 !!runoff sanity checks
                 if(csite%sfcwater_mass(i,ipa) > 100.0) then !!layer > 10cm
                    print*,"PREFLOOD"
                    print*,"patch water",csite%sfcwater_mass(i,ipa)
                    print*,"patch energy [J/kg]",csite%sfcwater_energy(i,ipa)
                    print*,"site water",qwout
                    print*,"site heat",qhout
                    print*,"curr_site",curr_site
                    print*,"runoff",csite%runoff(ipa)
                    print*,"heights",hts
                    print*,"surf_water_depth",surf_water_depth
                 end if
                 !! if surface water =< 0, there is no more surface water
                 if(csite%sfcwater_mass(i,ipa) <= 0.0) then
                    csite%sfcwater_mass(i,ipa)   = 0.0
                    csite%sfcwater_energy(i,ipa) = 0.0
                    csite%sfcwater_depth(i,ipa)  = 0.0
                    csite%nlev_sfcwater(ipa)     = 0
                 else
                    !!if surface water is too small, keep in equilibrium with soil
                    if(csite%sfcwater_mass(i,ipa) < water_stab_thresh   &
                         .and. i == 1) then
                       csite%sfcwater_energy(i,ipa) = cliq * (csite%soil_tempk(nzg,ipa)-tsupercool)
                    end if
                 end if
                 !! if sfcwater_depth < amount of water, set to water depth
                 if(csite%sfcwater_depth(i,ipa) < csite%sfcwater_mass(i,ipa)*0.001) then
                    csite%sfcwater_depth(i,ipa) = csite%sfcwater_mass(i,ipa)*0.001
                 end if
                 !! if sfcwater_energy < 0, set to freezing
                 if(csite%sfcwater_energy(i,ipa) < qicet3) then
                    csite%sfcwater_energy(i,ipa) = qliqt3
                 end if
                 !!check for NaN
                 if(csite%sfcwater_energy(i,ipa) /= csite%sfcwater_energy(i,ipa))then
                    call fatal_error('Failed sfcwater_energy sanity check in lsm_hyd' &
                                    ,'calcHydroSurface','lsm_hyd.f90')
                 end if
                 if(csite%soil_energy(nzg,ipa) /= csite%soil_energy(nzg,ipa)) then
                    call fatal_error('Failed soil_energy sanity check in lsm_hyd' &
                                    ,'calcHydroSurface','lsm_hyd.f90')
                 end if
              end do
           end do !patch
        end do !site
        deallocate(qwout)
        deallocate(qhout)

        !!Finally, calculate polygon-level runoff rate (runoff leaving terrestrial sites)
        cgrid%runoff(ipy) = 0.0
        do i=1,cpoly%nsites
              cgrid%runoff(ipy) = cgrid%runoff(ipy) + qwout(i)*cgrid%site_adjacency(i,(cpoly%nsites+1),ipy)
        end do
     end do !polygon

  end do !! grid
  return
end subroutine calcHydroSurface
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!MLO. Mike, I left the following two subroutines commented for now because the parallel POI! 
!     is still not ready. I think for simplicity the POI will still have a master node and !
!     slaves with different roles, which means that the subroutine would be really similar !
!     to these two, just with the new structure. But it would never reach this point in    !
!     this version because of the SOI parallel. And for regional runs and POI serial, all  !
!     sites will be in the same polygon so the previous should work fine there             !
!------------------------------------------------------------------------------------------!
!!!!!! subroutine calcHydroSurfaceMaster (polygon_list)
!!!!!! !! tallies overland water and heat fluxes among sites
!!!!!! !! on the head node when distributing by site (not patch or polygon)
!!!!!! !! STILL IN DEVELOPMENT - assumes only one polygon
!!!!!! !! 
!!!!!! !! input:  polygon_list
!!!!!! !! output: none
!!!!!! !! communicates with calcHydroSurfaceNode
!!!!!!   use site_type
!!!!!!   use hydrology_constants
!!!!!!   use hydrology_coms, only: useRUNOFF
!!!!!!   use ed_para_coms,only: mainnum,nmachs,machnum
!!!!!!   use grid_coms,only:ngrids
!!!!!!   implicit none
!!!!!!   include "mpif.h"
!!!!!!   type(plist),target,dimension(ngrids)::polygon_list
!!!!!!   type(polygon),pointer :: cpoly
!!!!!!   type(site),pointer :: cs
!!!!!!   integer :: nsites,cnt,tcnt,i,ierr,nm,tag=1212
!!!!!!   integer, allocatable :: sitenums(:),sitetmp(:)
!!!!!!   real, allocatable :: qwout(:),qhout(:),qwtmp(:),qhtmp(:)

!!!!!!   if(useRUNOFF .eq. 0) return

!!!!!!   cpoly => polygon_list(1)%first_polygon
!!!!!!   nsites = cpoly%nsites
!!!!!!   allocate(sitenums(1:nsites))
!!!!!!   allocate(sitetmp(1:nsites))
!!!!!!   allocate(qwout(1:nsites))
!!!!!!   allocate(qhout(1:nsites))
!!!!!!   allocate(qwtmp(1:nsites))
!!!!!!   allocate(qhtmp(1:nsites))
!!!!!!
!!!!!!   !first, broadcast number of sites
!!!!!!   call MPI_Bcast(nsites,1,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
!!!!!!
!!!!!!   !!Next, recieve fluxes from the nodes
!!!!!!   cnt=0
!!!!!!   do nm=1,nmachs
!!!!!!
!!!!!!      !!get fluxes
!!!!!!      call MPI_Recv(tcnt,1,MPI_INTEGER,machnum(nm),tag,MPI_COMM_WORLD,ierr) !number of sites I will be getting
!!!!!!      call MPI_Recv(sitetmp(1:tcnt),tcnt,MPI_INTEGER,machnum(nm),tag,MPI_COMM_WORLD,ierr) !!site numbers
!!!!!!      call MPI_Recv(qwtmp(1:tcnt),tcnt,MPI_REAL,machnum(nm),MPI_COMM_WORLD,ierr) !!water flux
!!!!!!      call MPI_Recv(qhtmp(1:tcnt),tcnt,MPI_REAL,machnum(nm),tag,MPI_COMM_WORLD,ierr) !!heat flux
!!!!!!   
!!!!!!      !!figure out where each goes
!!!!!!      do i = 1,tcnt
!!!!!!         cs => cpoly%first_site
!!!!!!         do while(associated(cs))
!!!!!!            if(cs%sitenum .eq. sitetmp(i)) exit
!!!!!!            cs => cs%next_site
!!!!!!         enddo
!!!!!!         sitenums(cnt+i) = cs%sitenum
!!!!!!         qwout(cnt+i) = qwtmp(i)
!!!!!!         qhout(cnt+i) = qhtmp(i)
!!!!!!      enddo
!!!!!!      cnt = cnt + tcnt
!!!!!!   enddo

!!!!!!   !!for output
!!!!!!   qw4out = qwout
!!!!!!   qh4out = qhout

!!!!!!   !!broadcast of fluxes among sites
!!!!!!   CALL MPI_Bcast(sitenums,nsites,MPI_INTEGER,mainnum,MPI_COMM_WORLD,ierr)
!!!!!!   CALL MPI_Bcast(qwout,nsites,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)
!!!!!!   CALL MPI_Bcast(qhout,nsites,MPI_REAL,mainnum,MPI_COMM_WORLD,ierr)

!!!!!!   !!Calculate polygon-level runoff rate (runoff leaving terrestrial sites)
!!!!!!   cpoly%runoff = 0.0
!!!!!!   do i=1,cpoly%nsites
!!!!!!      cpoly%runoff = cpoly%runoff + qwout(i)*cpoly%site_adjacency(i,(cpoly%nsites+1))
!!!!!!   enddo

!!!!!!   !! clean up
!!!!!!   deallocate(sitenums)
!!!!!!   deallocate(sitetmp)
!!!!!!   deallocate(qwout)
!!!!!!   deallocate(qhout)
!!!!!!   deallocate(qwtmp)
!!!!!!   deallocate(qhtmp)

!!!!!! end subroutine calcHydroSurfaceMaster


!!!!!! !! Calculates overland flow (runoff rate kg/s/m2)
!!!!!! !! on nodes when distributed by site
!!!!!! !! STILL IN DEVELOPMENT - only supports a single polygon
!!!!!! !! inputs: polygon_list
!!!!!! !! outputs: sets cp%runoff, cp%avg_runoff_heat
!!!!!! subroutine calcHydroSurfaceNode (polygon_list)
!!!!!!   use site_type
!!!!!!   use hydrology_constants
!!!!!!   use hydrology_coms, only: useRUNOFF, FracLiqRunoff,runoff_Vmax
!!!!!!   use soil_coms, only: nzg, water_stab_thresh
!!!!!!   use grid_coms,only: ngrids
!!!!!!   use ed_misc_coms, only: dtlsm
!!!!!!   use ed_node_coms, only: master_num
!!!!!!   use const_coms, only : g, 
!!!!!!   use therm_lib, only: qtk
!!!!!!   implicit none
!!!!!!   include "mpif.h"
!!!!!!   type(plist),dimension(ngrids) :: polygon_list
!!!!!!   type(polygon), pointer ::cpoly
!!!!!!   type(site), pointer ::cs
!!!!!!   type(patch), pointer :: cpatch
!!!!!! !  type(cohort), pointer :: cc
!!!!!!   integer :: tag=1212
!!!!!!   integer :: ierr
!!!!!!   integer :: top_surf_water     ! index of top surface water layer
!!!!!!   integer :: curr_site                ! site counter
!!!!!!   integer :: i
!!!!!!   integer :: nsites             ! number of sites in the polygon
!!!!!!   real, allocatable :: qwout(:) ! flux of water out of a site (kg/s/m2?)
!!!!!!   real, allocatable :: qhout(:) ! flux of heat out of a site ( units /s/m2? ******************)
!!!!!!   real :: tempk                 ! water temperature (K)
!!!!!!   real :: fracliq               ! water fraction liquid (proportion)
!!!!!!   real :: surf_water_depth      ! depth of liquid surface water (m)
!!!!!!   real :: surf_water_heat       ! energy of liquid surface water (J/m3 ??)
!!!!!!   real :: swd_i                 ! surface water depth in layer i (m)
!!!!!!   real :: flow_vel, flow_denom  ! overland flow velocity (m/s)
!!!!!!   real :: patch_water_out       ! patch water flux out (kg/m2/s)
!!!!!!   real :: patch_heat_out        ! patch heat flux out (J/m2/s??? ***************)
!!!!!!   real :: site_water_in         ! site water flux in (kg/m2/s)
!!!!!!   real :: site_heat_in          ! site heat flux in (J/m2/s)
!!!!!! !  real :: patch_water_in        ! patch water flux in (kg/m2/s)
!!!!!! !  real :: patch_heat_in         ! patch heat flux in (J/m2/s)
!!!!!!   real,dimension(11) :: hts     ! heights of surface water layers (m)
!!!!!!   


!!!!!!   if(useRUNOFF .eq. 0) return

!!!!!!   !! check to make sure we can handle the current set-up
!!!!!!   if(ngrids .gt. 1) then
!!!!!!      print*,"multiple grids not currently supported in calcHydroSurfaceNode"
!!!!!!      stop
!!!!!!   endif
!!!!!!   if(polygon_list(1)%npoly .gt. 1) then
!!!!!!      print*,"multiple polygons not currently supported in calcHydroSurfaceNode"
!!!!!!      stop
!!!!!!   endif
!!!!!!   if(polygon_list(1)%passmode .eq. 2) then
!!!!!!      if(polygon_list(1)%first_polygon%nsites .gt. 2) then
!!!!!!         print*,"PASSMODE 2 not supported for multi-site runs in calcHydroSurfaceNode"
!!!!!!         stop
!!!!!!      endif
!!!!!!      !!if single terrestrial site then do old runoff calculation
!!!!!!      print*,"PASSMODE 2 not supported in calcHydroSurfaceNode"
!!!!!!      return
!!!!!!   endif

!!!!!!   !!first, get the total number of sites from the master
!!!!!!   call MPI_Bcast(nsites,1,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr)

!!!!!!   cpoly => polygon_list(1)%first_polygon
!!!!!!   allocate(qwout(1:nsites))
!!!!!!   allocate(qhout(1:nsites))
!!!!!!   qwout = 0.0
!!!!!!   qhout = 0.0
!!!!!!   
!!!!!!   !! FIRST PASS, CALCULATE RUNOFF
!!!!!!   curr_site = 1
!!!!!!   cs => cpoly%first_site
!!!!!!   do while(associated(cs))

!!!!!!      !site specific parms
!!!!!!      cpoly%sitenums(curr_site) = cs%sitenum
!!!!!!      do while(associated(cpatch))
!!!!!!         !first check that there IS surface water
!!!!!!         top_surf_water = cpatch%nlev_sfcwater
!!!!!!         if(top_surf_water .gt. 0) then
!!!!!!            !second check that some is liquid
!!!!!!            call qtk(cpatch%sfcwater_energy(top_surf_water),tempk,fracliq)
!!!!!!            if(fracliq .gt. FracLiqRunoff) then
!!!!!!               !! water depth
!!!!!!               swd_i = cpatch%sfcwater_mass(top_surf_water)*0.001*fracliq
!!!!!!               surf_water_depth = swd_i
!!!!!!               surf_water_heat = swd_i*cliqvlme*(tempk-tsupercool)
!!!!!!               do i=(top_surf_water-1),1,-1
!!!!!!                  call qtk(cpatch%sfcwater_energy(i),tempk,fracliq)
!!!!!!                  swd_i = cpatch%sfcwater_mass(top_surf_water)*0.001*fracliq
!!!!!!                  surf_water_depth = surf_water_depth + swd_i
!!!!!!                  surf_water_heat = surf_water_heat + swd_i*cliqvlme*(tempk-tsupercool)
!!!!!!               enddo
!!!!!!               
!!!!!!               !! calculate flow velocity (m/s)                 
!!!!!!               flow_denom = 1+cpatch%runoff_a(2)*surf_water_depth**(4/3) & 
!!!!!!                    - cpatch%runoff_a(3)*surf_water_depth**(7/3)
!!!!!!               if(flow_denom .gt. 1.0) then
!!!!!!                  flow_vel = surf_water_depth**(2/3)*cpatch%runoff_a(1)/sqrt(flow_denom)
!!!!!!               else
!!!!!!                  flow_vel = surf_water_depth**(2/3)*cpatch%runoff_a(1) !! asymptotic vel w/o vegetation
!!!!!!               endif
!!!!!!               flow_vel = min(flow_vel,runoff_vmax) !! clamp runoff velocity to maximum sensible value
!!!!!!               
!!!!!!               !!patch level flux out (kg/m2/s)
!!!!!!               patch_water_out = flow_vel*surf_water_depth*1000
!!!!!!               patch_heat_out = flow_vel*surf_water_heat  !energy/m2/s
!!!!!!               !!update vars
!!!!!!               cpatch%runoff = patch_water_out
!!!!!!               cpatch%avg_runoff_heat = patch_heat_out
!!!!!!               qwout(curr_site) = qwout(curr_site) + patch_water_out*cpatch%area
!!!!!!               qhout(curr_site) = qhout(curr_site) + patch_heat_out*cpatch%area
!!!!!!            endif
!!!!!!         endif
!!!!!!         cp => cpatch%younger
!!!!!!      enddo
!!!!!!      cs => cs%hydro_next
!!!!!!      curr_site = curr_site + 1
!!!!!!   enddo
!!!!!!   curr_site = curr_site - 1

!!!!!!   !!send my fluxes to the master
!!!!!!   call MPI_Send(curr_site,1,MPI_INTEGER,master_num,tag,MPI_COMM_WORLD,ierr) !number of sites I have
!!!!!!   call MPI_Send(cpoly%sitenums(1:curr_site),curr_site,MPI_INTEGER,master_num,tag,MPI_COMM_WORLD,ierr) !!my site numbers
!!!!!!   call MPI_Send(qwout(1:curr_site),curr_site,MPI_REAL,master_num,tag,MPI_COMM_WORLD,ierr) !!my water flux
!!!!!!   call MPI_Send(qhout(1:curr_site),curr_site,MPI_REAL,master_num,tag,MPI_COMM_WORLD,ierr) !!my heat flux
!!!!!!   

!!!!!!   !!recieve broadcast of fluxes among sites
!!!!!!   CALL MPI_Bcast(cpoly%sitenums,nsites,MPI_INTEGER,master_num,MPI_COMM_WORLD,ierr) 
!!!!!!   CALL MPI_Bcast(qwout,nsites,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)
!!!!!!   CALL MPI_Bcast(qhout,nsites,MPI_REAL,master_num,MPI_COMM_WORLD,ierr)


!!!!!!   !! SECOND PASS: CALCULATE RUN-IN
!!!!!!   cs => cpoly%first_site
!!!!!!   do while(associated(cs))

!!!!!!      !!determine my site
!!!!!!      do curr_site=1,nsites
!!!!!!         if(cpoly%sitenums(curr_site) == cs%sitenum) exit
!!!!!!      enddo
!!!!!!      !!calc site-level run-in
!!!!!!      site_water_in = 0
!!!!!!      site_heat_in = 0
!!!!!!      do i=1,nsites
!!!!!!         if(i .ne. curr_site) then
!!!!!!            site_water_in = site_water_in + qwout(i)*cpoly%site_adjacency(i,curr_site)
!!!!!!            site_heat_in = site_heat_in + qhout(i)*cpoly%site_adjacency(i,curr_site)
!!!!!!         endif
!!!!!!      enddo

!!!!!!      !!calc site-level run-off
!!!!!!      cs%runoff = 0.0
!!!!!!      cs%avg_runoff_heat = 0.0
!!!!!!      do i=1,(nsites+1)
!!!!!!         if(i .ne. curr_site) then
!!!!!!            cs%runoff = cs%runoff + qwout(curr_site)*cpoly%site_adjacency(curr_site,i)
!!!!!!            cs%avg_runoff_heat = cs%avg_runoff_heat + qhout(curr_site)*cpoly%site_adjacency(curr_site,i)
!!!!!!         endif
!!!!!!      enddo
!!!!!!  
!!!!!!      !! calc patch layer run in
!!!!!!      cp => cs%oldest_patch
!!!!!!      do while(associated(cp))
!!!!!!      
!!!!!!         !!calc runoff rates
!!!!!!         cpatch%runoff = (1-cpoly%site_adjacency(curr_site,curr_site))*cpatch%runoff - site_water_in 
!!!!!!         cpatch%avg_runoff_heat = (1-cpoly%site_adjacency(curr_site,curr_site))*cpatch%avg_runoff_heat - site_heat_in

!!!!!!         surf_water_depth = 0.0
!!!!!!         top_surf_water = cpatch%nlev_sfcwater
!!!!!!         !!compute fraction to distribute to each layer
!!!!!!         do i=1,top_surf_water
!!!!!!            call qtk(cpatch%sfcwater_energy(i),tempk,fracliq)
!!!!!!            hts(i) = cpatch%sfcwater_mass(i)*0.001*fracliq
!!!!!!            surf_water_depth =  surf_water_depth + hts(i)
!!!!!!         enddo

!!!!!!         !! calc actual runoff/run-on
!!!!!!         do i=1,top_surf_water

!!!!!!            ! Convert to J/m2
!!!!!!            cpatch%sfcwater_energy(i) = cpatch%sfcwater_energy(i) * cpatch%sfcwater_mass(i)
!!!!!!            cpatch%sfcwater_mass(i) = cpatch%sfcwater_mass(i)-cpatch%runoff*dtlsm*hts(i)/surf_water_depth
!!!!!!            ! Subtract and convert back to J/kg.
!!!!!!            if(cpatch%sfcwater_mass(i) > 1.0e-10)then
!!!!!!               cpatch%sfcwater_energy(i) = (cpatch%sfcwater_energy(i) -  &
!!!!!!                    cpatch%avg_runoff_heat*dtlsm*hts(i)/surf_water_depth) /   &
!!!!!!                    cpatch%sfcwater_mass(i)
!!!!!!            else
!!!!!!               cpatch%sfcwater_energy(i) = 0.0
!!!!!!            endif

!!!!!!            !!sanity checks
!!!!!!            if(cpatch%sfcwater_mass(i).le.0) then
!!!!!!               cpatch%sfcwater_mass(i) = 0.0
!!!!!!               cpatch%sfcwater_energy(i) = 0.0
!!!!!!               cpatch%sfcwater_depth(i) = 0.0
!!!!!!               cpatch%nlev_sfcwater=0
!!!!!!            else
!!!!!!               !!if surface water is too small, keep in equilibrium with soil
!!!!!!               if(cpatch%sfcwater_mass(i) .lt. water_stab_thresh .and. i .eq. 1) then
!!!!!!                  cpatch%sfcwater_energy(i) = cliq*(cpatch%soil_tempk(nzg)-tsupercool)
!!!!!!               endif
!!!!!!            endif
!!!!!!            if(cpatch%sfcwater_depth(i) .lt. cpatch%sfcwater_mass(i)*0.001) then
!!!!!!               cpatch%sfcwater_depth(i) = cpatch%sfcwater_mass(i)*0.001
!!!!!!            endif
!!!!!!            if(cpatch%sfcwater_energy(i) .lt. 0.0) then
!!!!!!               cpatch%sfcwater_energy(i) = qliqt3
!!!!!!            endif
!!!!!!            !!check for NaN
!!!!!!            if(cpatch%sfcwater_energy(i) /= cpatch%sfcwater_energy(i))then
!!!!!!               print*,"failed sfcwater_energy sanity check in lsm_hyd"
!!!!!!               stop
!!!!!!            endif
!!!!!!            if(cpatch%soil_energy(nzg) /= cpatch%soil_energy(nzg))then
!!!!!!               print*,"failed soil_energy sanity check in lsm_hyd"
!!!!!!               stop
!!!!!!            endif
!!!!!!         enddo
!!!!!!         cpatch => cpatch%younger
!!!!!!      enddo
!!!!!!      cs => cs%next_site
!!!!!!   enddo
!!!!!!   
!!!!!!   deallocate(qwout)
!!!!!!   deallocate(qhout)

!!!!!! end subroutine calcHydroSurfaceNode
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine updateHydroParms (cgrid)
  !! Calculates hydrology parameters
  !! Updates runoff roughness for vegetation growth
  !! based on USGS Water-Supply Paper 2339
  !! inputs: 
  !! output: updates cgrid%runoff_a(:) - vector of runoff parameters
  !! only needs to be updated after changes in cohort density/size 
  !! (e.g. annually and/or after disturbance)
  !! needs tables for n3,n4  
  !! needs CWD implementation for n3
  use ed_state_vars, only : edtype, polygontype,sitetype,patchtype
  use hydrology_constants
  use hydrology_coms, only: useRUNOFF, grassLAImax
  use consts_coms, only: grav
  use pft_coms, only: include_pft_ag
  implicit none
  type(edtype)      , target  :: cgrid           ! Alias for current grid
  type(polygontype) , pointer :: cpoly           ! Alias for current polygon
  type(sitetype)    , pointer :: csite           ! Alias for current site
  type(patchtype)   , pointer :: cpatch          ! Alias for current patch
  integer                     :: ipy,isi,ipa,ico ! Polygon, site, patch and cohort counter
  real                        :: nb              ! soil roughness
  real                        :: n3              ! obstruction correction (e.g. CWD)
  real                        :: n4              ! non-woody vegetation correction
  real                        :: n0              ! total non-tree roughness
  real                        :: beta            ! slope (degrees)
  real                        :: sigma           ! patch 
!  integer :: ifm ! grid counter

  if(useRUNOFF == 0) return
  write (unit=*,fmt='(a)') 'What am I doing here????'
  !! Surface Runoff Parameterization
  !! Parameters for Manning's Formula
  !! based on USGS Water-Supply Paper 2339
  !! used in the calculation v = h^(2/3)*a0/sqrt(1+a1*h^(4/3)-a2*h^(7/3)) (m/s)
  !! only needs to be update after changes in cohort density/size 
  !! (e.g. annually and/or after disturbance)

  do ipy=1,cgrid%npolygons
     cpoly => cgrid%polygon(ipy)
     do isi=1,cpoly%nsites
        csite => cpoly%site(isi)
        !site specific parms
        beta = cpoly%slope(isi) !slope (degrees)
        nb   = 0.025 ! soil roughness
        do ipa=1,csite%npatches
           cpatch => csite%patch(ipa)
           ! obstruction correction, approx 0.1* propr. x-section area obstructed by CWD
           n3    = 0.0 
           ! non-woody veg correction, approx 0.1*(grass biomass)/(max grass biomass)
           n4    = 0.0   
           sigma = 0.0
           do ico=1,cpatch%ncohorts
              if(include_pft_ag(cpatch%pft(ico))) then
                 !! update non-woody correction
                 n4 = n4 + 0.01*cpatch%LAI(ico)/GrassLAImax
              else
                 !! update woody cross-sectional area (ft)
                 sigma = sigma + cpatch%nplant(ico)*cpatch%DBH(ico)*m2f
              end if
           end do ! cohort
           n0 = nb + n3 + n4 
!           cp%runoff_a(1) = 1.486*m2f**(-1/3)*sqrt(tan(beta))/n0
           csite%runoff_a(1,ipa) = m2f**(-1/3)*sqrt(tan(beta))/n0
           csite%runoff_a(2,ipa) = c1 * m2f**(4/3)*sigma*1.49*1.49/(2.*grav*m2f*n0*n0)
           csite%runoff_a(3,ipa) = c2 * m2f**(7/3)*sigma*1.49*1.49/(2.*grav*m2f*n0*n0)
        end do !patch
     end do !site 
  end do !polygon
!  write(unit=*,fmt=*) "END updateHydroParms"
  return
end subroutine updateHydroParms
!==========================================================================================!
!==========================================================================================!

