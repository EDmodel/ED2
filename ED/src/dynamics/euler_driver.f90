subroutine euler_timestep_ar(cgrid)

  use ed_state_vars,only: edtype,polygontype, &
       sitetype,patchtype
  use misc_coms, only: dtlsm
  use soil_coms, only: soil_rough
  use consts_coms, only: cp,mmdryi

  implicit none

  type(edtype),target       :: cgrid
  type(polygontype),pointer :: cpoly
  type(sitetype),pointer    :: csite
  type(patchtype),pointer   :: cpatch
  integer :: ipy,isi,ipa
  real :: thetaatm
  real :: thetacan
  real, parameter :: snowrough=0.001
  real :: hgtmin


  do ipy = 1,cgrid%npolygons

     cpoly => cgrid%polygon(ipy)

     do isi = 1,cpoly%nsites

        csite => cpoly%site(isi)
        
        thetaatm = cpoly%met(isi)%atm_tmp / cpoly%met(isi)%exner

        do ipa = 1,csite%npatches

           cpatch => csite%patch(ipa)

           thetacan = csite%can_temp(ipa) / cpoly%met(isi)%exner
           csite%rough(ipa) = max(soil_rough, csite%veg_rough(ipa)) *   &
                (1.0 - csite%snowfac(ipa)) + snowrough
           hgtmin = max(csite%rough(ipa) + 1.0e-3,   &
                max(50.0, cpoly%met(isi)%geoht))
           if(csite%can_temp(ipa) < cpoly%met(isi)%atm_tmp)then
              cpoly%met(isi)%vels = cpoly%met(isi)%vels_stab
           else
              cpoly%met(isi)%vels = cpoly%met(isi)%vels_unstab
           endif

           call ed_stars(thetaatm, cpoly%met(isi)%atm_shv,   &
                cpoly%met(isi)%atm_co2, thetacan, hgtmin, cpoly%met(isi)%vels, &
                csite%rough(ipa), csite%ustar(ipa), csite%rstar(ipa),   &
                csite%tstar(ipa), csite%cstar(ipa), csite%can_shv(ipa), &
                csite%can_co2(ipa))

           csite%co2budget_loss2atm(ipa) = - cpoly%met(isi)%rhos *   &
                csite%ustar(ipa) * csite%cstar(ipa) * mmdryi * dtlsm
           csite%ebudget_loss2atm(ipa) = - cp * cpoly%met(isi)%rhos * &
                csite%ustar(ipa) * csite%tstar(ipa) * cpoly%met(isi)%exner * dtlsm
           csite%wbudget_loss2atm(ipa) = - cpoly%met(isi)%rhos *  &
                csite%ustar(ipa) * csite%rstar(ipa) * dtlsm

           ! This is like the LEAF-3 implemented in OLAM.
           call leaf3_land_ar(csite,ipa,csite%nlev_sfcwater(ipa),   &
                csite%ntext_soil(:,ipa), csite%soil_water(:,ipa),   &
                csite%soil_energy(:,ipa), csite%sfcwater_mass(:,ipa),   &
                csite%sfcwater_energy(:,ipa), csite%sfcwater_depth(:,ipa),   &
                csite%rshort_s(:,ipa), csite%rshort_g(ipa), csite%rlong_s(ipa),  &
                csite%rlong_g(ipa), csite%veg_height(ipa), csite%veg_rough(ipa),  &
                csite%lai(ipa), csite%can_depth(ipa), cpoly%met(isi)%rhos,   &
                cpoly%met(isi)%vels, cpoly%met(isi)%prss,   &
                cpoly%met(isi)%pcpg * dtlsm,   &
                cpoly%met(isi)%qpcpg * dtlsm,  &
                cpoly%met(isi)%dpcpg * dtlsm,   &
                csite%ebudget_loss2atm(ipa),   &
                csite%wbudget_loss2atm(ipa),  &
                csite%co2budget_loss2atm(ipa),   &
                csite%ustar(ipa), csite%snowfac(ipa), csite%surface_ssh(ipa),  &
                csite%ground_shv(ipa), csite%can_temp(ipa), csite%can_shv(ipa),  &
                csite%can_co2(ipa), csite%rough(ipa), cpoly%lsl(isi),  &
                cpoly%leaf_aging_factor(:,isi),   &
                cpoly%green_leaf_factor(:,isi))

        enddo
     enddo
  enddo

  return
end subroutine euler_timestep_ar

!*****************************************************************************

subroutine leaf3_land_ar(csite,ipa, nlev_sfcwater,   &
     ntext_soil, soil_water,   &
     soil_energy, sfcwater_mass,   &
     sfcwater_energy, sfcwater_depth,   &
     rshort_s, rshort_g, rlong_s,  &
     rlong_g, veg_height, veg_rough,  &
     veg_tai, can_depth, rhos,   &
     vels, prss, pcpg,   &
     qpcpg, dpcpg,   &
     sxfer_t, sxfer_r, sxfer_c,   &
     ustar, snowfac, surface_ssh,  &
     ground_shv, can_temp, can_shv,  &
     can_co2, rough, lsl, leaf_aging_factor, green_leaf_factor)

  
  use ed_state_vars, only: sitetype
  use grid_coms, only: nzg, nzs
  use soil_coms, only: soil_rough, snow_rough, soil
  use misc_coms, only: dtlsm
  use max_dims, only: n_pft
  use therm_lib, only: qtk,qwtk8
  use ed_therm_lib,only:ed_grndvap
  use consts_coms, only: wdns

  implicit none

  type(sitetype), target :: csite
  integer, intent(in) :: ipa
  integer, intent(inout) :: nlev_sfcwater
  integer, dimension(nzg), intent(in) :: ntext_soil
  real(kind=8), dimension(nzg), intent(inout) :: soil_water
  real, dimension(nzg), intent(inout) :: soil_energy
  real, dimension(nzs), intent(inout) :: sfcwater_mass
  real, dimension(nzs), intent(inout) :: sfcwater_energy
  real, dimension(nzs), intent(inout) :: sfcwater_depth
  real, intent(inout), dimension(nzs) :: rshort_s
  real, intent(in) :: rshort_g
  real, intent(in) :: rlong_s
  real, intent(in) :: rlong_g
  real, intent(in) :: veg_height
  real, intent(inout) :: veg_rough
  real, intent(inout) :: veg_tai
  real, intent(in) :: can_depth
  real, intent(in) :: rhos
  real, intent(in) :: vels
  real, intent(in) :: prss
  real, intent(in) :: pcpg
  real, intent(in) :: qpcpg
  real, intent(in) :: dpcpg
  real, intent(in) :: sxfer_t
  real, intent(in) :: sxfer_r
  real, intent(in) :: sxfer_c
  real, intent(in) :: ustar
  real, intent(in) :: snowfac
  real, intent(inout) :: surface_ssh
  real, intent(inout) :: ground_shv
  real, intent(inout) :: can_temp
  real, intent(inout) :: can_shv
  real, intent(inout) :: can_co2
  real, intent(inout) :: rough
  integer, intent(in) :: lsl
  real, intent(in), dimension(n_pft) :: leaf_aging_factor, green_leaf_factor

  real, dimension(nzg) :: soil_tempk
  real, dimension(nzg) :: soil_fracliq
  real :: soil_rfactor    (nzg) ! soil thermal resistivity [K m^2/W]
  real, dimension(nzs) :: sfcwater_tempk
  real, dimension(nzs) :: sfcwater_fracliq
  real :: hxferg        (nzg+1) ! heat xfer between soil layers [J/m^2]

  real :: wxfer       (nzg+1) ! soil water xfer [m]
  real :: qwxfer      (nzg+1) ! soil energy xfer from water xfer [J/m^2] 
  real :: psiplusz    (nzg)   ! soil water potential (including grav) [m]
  real :: hydraul_cond(nzg)   ! soil hydraulic conductivity [m/s]

  real :: surface_temp,surface_fliq

  integer :: k     ! vertical index over soil layers
  integer :: nlsw1 ! maximum of (1,nlev_sfcwater)
  
  integer :: ktrans ! vertical index of soil layer supplying transpiration
  
  real :: transp  ! transpiration xfer this LEAF timestep [kg/m^2]
  real, dimension(nzg) :: ed_transp ! transpired water from each soil level; ED2 cells only [kg/m^2]
  real :: hxfergc ! heat xfer from ground (soil) to can_air this step [J/m^2]
  real :: wxfergc ! vapor xfer from ground (soil) to can_air this step [kg_vap/m^2  ]
  real :: hxfersc ! heat xfer from sfcwater to can_air this step [J/m^2]
  real :: wxfersc ! vapor xfer from sfcwater to can_air this step [kg_vap/m^2]
  real :: wshed   ! water shed from veg this timestep [kg/m^2]
  real :: qwshed  ! water energy shed from veg this timestep [J/m^2]
  real :: rdi     ! (soil or surface water)-to-can_air conductance [m/s]
  

  !=================================================================
  ! TEMPORARY RUNOFF VARIABLES
  real :: runoff
  real :: qrunoff
  !=================================================================

  soil_tempk(lsl:nzg) = csite%soil_tempk(lsl:nzg,ipa)
  soil_fracliq(lsl:nzg) = csite%soil_fracliq(lsl:nzg,ipa)
  sfcwater_tempk(1:nlev_sfcwater) =   &
       csite%sfcwater_tempk(1:nlev_sfcwater,ipa)
  sfcwater_fracliq(1:nlev_sfcwater) =   &
       csite%sfcwater_fracliq(1:nlev_sfcwater,ipa)

  ! Evaluate turbulent exchanges of heat and moisture between vegetation 
  ! and canopy air and also between soil or snow surface and canopy air.  
  ! Evaluate transfer of precipitation moisture and heat to vegetation and 
  ! shed from vegetation to surface.  Update vegetation and canopy air 
  ! temperatures resulting from these plus radiative fluxes.

  nlsw1 = max(csite%nlev_sfcwater(ipa),1)

  call canopy_ar(                                        &
       nlev_sfcwater,         ntext_soil,            &
       ktrans,                &
       soil_water,            soil_fracliq,          &
       soil_tempk    (nzg),   sfcwater_mass (nlsw1), &
       sfcwater_tempk(nlsw1), veg_height,            &
       veg_rough,             veg_tai,               &
       can_depth,             rhos,                  &
       vels,                  prss,                  &
       pcpg,                  qpcpg,                 &
       sxfer_t,               sxfer_r,               &
       sxfer_c,               ustar,                 &
       snowfac,                                   &
       surface_ssh,           ground_shv,            &
       can_temp,              can_shv,               &
       can_co2,                                      &
       wshed,                 qwshed,                &
       transp,                  &
       hxfergc,               wxfergc,               &
       hxfersc,               wxfersc,               &
       rdi,                   &
       lsl,                   ed_transp,             &
       csite,                                        &
       ipa,                                          &
       leaf_aging_factor, green_leaf_factor)


  ! CALL SFCWATER:
  !  1. Compute soil and sfcwater heat conductivities
  !  2. Compute surface heat flux for top layer of soil or sfcwater
  !  3. Evaluate internal and bottom fluxes for sfcwater
  !  4. Update sfcwater layer energies due to heat flux and solar radiation
  !  5. Evaluate melting and percolation of liquid through sfcwater layers
  
  call sfcwater(  &
       nlev_sfcwater,    ntext_soil,       &
       soil_rfactor,     soil_water,       &
       soil_energy,      sfcwater_mass,    &
       sfcwater_energy,  sfcwater_depth,   &
       soil_tempk,       soil_fracliq ,    &
       sfcwater_tempk,   sfcwater_fracliq, &
       rshort_s,         hxfersc,          &
       wxfersc,          rlong_g,          &
       rlong_s,          pcpg,             &
       qpcpg,            dpcpg,            &
       wshed,            qwshed,           &
       lsl            )
  
  call soil_euler_ar(                       &
       nlev_sfcwater, &
       ntext_soil,   ktrans,        &
       soil_tempk,   soil_fracliq,  &
       soil_rfactor, hxfergc,       &
       wxfergc,      rshort_g,      &
       rlong_g,      transp,        &
       soil_energy,  soil_water,    &
       hxferg,       wxfer,         &
       qwxfer,       psiplusz,      &
       hydraul_cond, lsl,           &
       ed_transp,    csite, ipa)
  
  ! Compute ground vap mxrat for next timestep; put into ground_ssh.
  
  nlsw1 = max(nlev_sfcwater,1)
  
  call ed_grndvap(nlev_sfcwater,&
       ntext_soil       (nzg),  &
       soil_water       (nzg),  &
       soil_energy      (nzg),  &
       sfcwater_energy(nlsw1),  &
       rhos,  &
       can_shv,  &
       ground_shv,  &
       surface_ssh, &
       surface_temp,surface_fliq)
  
  do k = lsl,nzg
     call qwtk8(soil_energy(k),soil_water(k)*dble(wdns),  &
          soil(ntext_soil(k))%slcpd,soil_tempk(k),soil_fracliq(k))
  enddo
  
  ! Diagnose surface water temperature and liquid fraction
  
  do k = 1,nlev_sfcwater
     call qtk(sfcwater_energy(k),sfcwater_tempk(k),sfcwater_fracliq(k))
  enddo
  
  do k = lsl,nzg
     csite%soil_tempk(k,ipa) = soil_tempk(k)
     csite%soil_fracliq(k,ipa) = soil_fracliq(k)
  enddo
  
  do k = 1,nlev_sfcwater
     csite%sfcwater_tempk(k,ipa) = sfcwater_tempk(k)
     csite%sfcwater_fracliq(k,ipa) = sfcwater_fracliq(k)
  enddo
  
!-----------------------------------------------------------------------------
! TEMPORARY UNTIL FULL LEAF-HYDRO MODEL IS COMPLETED WITH STREAM/RIVER RUNOFF:
! Simple representation of runoff

  call remove_runoff(nlev_sfcwater,    &
       sfcwater_fracliq, sfcwater_mass,    &
       sfcwater_tempk,   sfcwater_energy,  &
       sfcwater_depth,   runoff,           &
       qrunoff                             )
!-----------------------------------------------------------------------------

  csite%wbudget_loss2runoff(ipa) = runoff
  csite%ebudget_loss2runoff(ipa) = qrunoff

  return
end subroutine leaf3_land_ar

!***************************************************************************

subroutine canopy_ar(nlev_sfcwater, ntext_soil, ktrans,   &
                  soil_water, soil_fracliq, soil_tempk,                   &
                  sfcwater_mass, sfcwater_tempk,                          &
                  veg_height, veg_rough, veg_tai,                &
                  can_depth,                                     &
                  rhos, vels, prss, pcpg, qpcpg,                          &
                  sxfer_t, sxfer_r, sxfer_c, ustar,     &
                  snowfac, surface_ssh, ground_shv,                   &
                  can_temp, can_shv, can_co2,                 &
                  wshed, qwshed, transp,                      &
                  hxfergc, wxfergc, hxfersc, wxfersc,            &
                  rdi,                         &
                  lsl, ed_transp, csite, ipa,   &
                  leaf_aging_factor, green_leaf_factor)


use soil_coms,   only: soil_rough, dslz
use grid_coms, only: nzg
use misc_coms, only: dtlsm
use consts_coms, only: cp, vonk, alvl, cliq, cice, alli, rvap

use ed_state_vars,only:sitetype,patchtype
use max_dims, only: n_pft

implicit none

integer, intent(in)  :: nlev_sfcwater   ! # active levels of surface water
integer, intent(in)  :: ntext_soil(nzg) ! soil textural class
integer, intent(out) :: ktrans          ! k index of soil layer supplying transp

real(kind=8), intent(in) :: soil_water(nzg)   ! soil water content [vol_water/vol_tot]
real, intent(in) :: soil_fracliq(nzg) ! fraction of soil moisture in liquid phase
real, intent(in) :: soil_tempk        ! soil temp [K]
real, intent(in) :: sfcwater_mass     ! surface water mass [kg/m^2]
real, intent(in) :: sfcwater_tempk    ! surface water temp [K]
real, intent(in) :: veg_height  ! veg height [m]
real, intent(in) :: veg_rough   ! veg roughess height [m]
real, intent(in) :: veg_tai     ! veg total area index
real, intent(in) :: can_depth   ! canopy depth for heat and vap capacity [m]
real, intent(in) :: rhos        ! atmospheric air density [kg/m^3]
real, intent(in) :: vels        ! surface wind speed [m/s]
real, intent(in) :: prss        ! air pressure [Pa]
real, intent(in) :: pcpg        ! new precip amount this leaf timestep [kg/m^2]
real, intent(in) :: qpcpg       ! new precip energy this leaf timestep [J/m^2]
real, intent(in) :: sxfer_t     ! surface heat xfer this step [kg_air K/m^2]
real, intent(in) :: sxfer_r     ! surface vapor xfer this step [kg_vap/m^2]
real, intent(in) :: sxfer_c     ! surface CO2 xfer this step [umol/m^2]
real, intent(in) :: ustar       ! friction velocity [m/s]
real, intent(in) :: snowfac     ! fractional veg burial by snowcover
real, intent(in) :: surface_ssh ! surface sat spec hum [kg_vap/kg_air]
real, intent(in) :: ground_shv  ! soil vapor spec hum [kg_vap/kg_air]

real, intent(inout) :: can_temp    ! canopy air temp [K]
real, intent(inout) :: can_shv     ! canopy air vapor spec hum [kg_vap/kg_air]
real, intent(inout) :: can_co2     ! canopy air CO2 mixing ratio [umol/mol]

real, intent(out) :: hxfergc ! soil-to-can_air heat xfer this step [J/m^2]
real, intent(out) :: wxfergc ! soil-to-can_air vap xfer this step [kg_vap/m^2]
real, intent(out) :: hxfersc ! sfc_water-to-can_air heat xfer this step [J/m^2]
real, intent(out) :: wxfersc ! sfc_water-to-can_air vap xfer this step [kg_vap/m^2]
real, intent(out) :: wshed   ! water shed from veg this LEAF timestep [kg/m^2]
real, intent(out) :: qwshed  ! water energy shed from veg this LEAF timestep [J/m^2]
real, intent(out) :: transp  ! transpiration xfer this LEAF timestep [kg_vap/m^2]
real, intent(out) :: rdi     ! (soil or surface water)-to-can_air conductance [m/s]

real, dimension(nzg) :: ed_transp ! transpired water from each soil level; ED2 cells only [kg/m^2]

type(sitetype),target :: csite
integer :: ipa

real, intent(in), dimension(n_pft) :: leaf_aging_factor, green_leaf_factor

! Local parameters

real, parameter :: exar = 2.5  ! for computing rasveg
real, parameter :: covr = 2.16 ! scaling tai value for computing wtveg
real, parameter :: c1 = 116.6  ! for computing rb

!     Note: c1=261.5*sqrt((1.-exp(-2.*exar))/(2.*exar)) 
!     from Lee's dissertation, Eq. 3.36.  The factor of 261.5 is
!     100 * ln((h-d)/zo) / vonk   where d = .63 * h and zo = .13 * h.
!     The factor of 100 is 1/L in Eq. 3.37.  Thus, c1 * ustar is the
!     total expression inside the radical in Eq. 3.37.
!bob      parameter(exar=3.5,covr=2.16,c1=98.8)

! Intercept, slope parameters for stomatal conductance factors

real, parameter :: brad = 196.   , srad = .047    ! for s/w radiative flux
real, parameter :: btlo = 281.5  , stlo = .26     ! for low canopy temperature
real, parameter :: bthi = 310.1  , sthi = -.124   ! for high canopy temperature
real, parameter :: bvpd = 4850.  , svpd = -.0051  ! for vapor pressure deficit
real, parameter :: bswp = -1.07e6, sswp = 7.42e-6 ! for soil water potential

! Local variables

real :: factv       ! for computing rasveg
real :: aux         ! for computing rasveg
real :: rasveg      ! full-veg value of rd [s/m]
real :: wtveg       ! weighting of rasveg in computing rd
real :: zognd       ! soil roughness height [m]
real :: zoveg       ! vegetation roughness height [m]
real :: zdisp       ! vegetation displacement height remaining after snowcover [m]
real :: zveg        ! vegetation height remaining after snowcover [m]
real :: canair
real :: canhcap

integer, intent(in) :: lsl

! Canopy air (vapor) and heat capacity

canair = rhos * can_depth
canhcap = cp * canair

! Initialize wshed = qwshed = 0

wshed  = 0.
qwshed = 0.

! If vegetation is sufficiently abundant and not covered by snow,
! COMPUTE CANOPY XFER WITH VEGETATION INFLUENCE

! Compute ground-canopy resistance rd.  Assume zognd not affected by snow.
! Assume (zoveg,zdisp) decrease linearly with snow depth, attaining
! the values (zognd,0) when veg covered.

zognd = soil_rough
zoveg = veg_rough * (1. - snowfac) + zognd * snowfac
zveg  = veg_height * (1. - snowfac)
zdisp = zveg * .63
!bob      rasgnd = log(zts / zognd) * log((zdisp + zoveg) / zognd)
!bob     +      / (vonk * vonk * vels)

! Aerodynamic resistance (rd) between surface and canopy air are weighted
! between areas without and with vegetation.

!bob   factv  = log(zts / zoveg) / (vonk * vonk * vels)
factv  = 1. / (vonk * ustar)
aux    = exp(exar * (1. - (zdisp + zoveg) / zveg))
rasveg = factv * zveg / (exar * (zveg - zdisp)) * (exp(exar) - aux)
wtveg  = max(0.,min(1., 1.1 * veg_tai / covr))
rdi    = ustar / (5. * (1. - wtveg) + ustar * rasveg * wtveg)
   
! Check if any surface water layers currently exist

if (nlev_sfcwater >= 1) then

! If surface water is present, compute heat and vapor xfer between 
! surface water and canopy air.

   hxfergc = 0.
   hxfersc = cp * rhos * dtlsm * rdi * (sfcwater_tempk - can_temp)
   
   wxfergc = 0.
   wxfersc = min(sfcwater_mass,                                 &
        rhos * dtlsm * rdi * (surface_ssh - can_shv) )
   
else

! If surface water is absent, compute heat xfer between soil and canopy air

   hxfergc = cp * rhos * dtlsm * rdi * (soil_tempk - can_temp)
   hxfersc = 0.

! Compare saturation vapor specific humidity at soil temperature against 
! canopy air vapor specific humidity

   if (surface_ssh < can_shv) then

! If saturation vapor specific humidity at soil temperature is less than 
! canopy air vapor specific humidity, compute dew formation contribution 
! to surface water

      wxfergc = 0.
      wxfersc = rhos * dtlsm * rdi * (surface_ssh - can_shv)
      
   elseif (ground_shv > can_shv) then
      
! Else, if equilibrium soil vapor specific humidity is greater than 
! canopy air vapor specific humidity, compute evaporation from soil

      wxfergc = max(0.,rhos * dtlsm * rdi * (ground_shv - can_shv))
      wxfersc = 0.
      
   else

! If neither of the above is true, both vapor fluxes are zero.

      wxfergc = 0.
      wxfersc = 0.
      
   endif

endif

ktrans = 0 
transp = 0.
ed_transp(:) = 0.

call canopy_update_euler_ar(csite,ipa, vels, rhos, prss, pcpg, qpcpg,   &
     wshed, qwshed, canair, canhcap, dtlsm, hxfergc, sxfer_t,  &
     wxfergc, hxfersc, wxfersc, sxfer_r, ed_transp, ntext_soil,  &
     soil_water, soil_fracliq, lsl, leaf_aging_factor, green_leaf_factor,  &
     sxfer_c)


return
end subroutine canopy_ar

!***************************************************************************
subroutine sfcwater(nlev_sfcwater,ntext_soil,                       &
                    soil_rfactor, soil_water, soil_energy,                 &
                    sfcwater_mass, sfcwater_energy, sfcwater_depth,        &
                    soil_tempk, soil_fracliq, sfcwater_tempk,              &
                    sfcwater_fracliq, rshort_s,                            &
                    hxfersc, wxfersc, rlong_g, rlong_s,                    &
                    pcpg, qpcpg, dpcpg, wshed, qwshed,                     &
                    lsl                               )

use soil_coms, only: slz, dslz, dslzi, dslzo2, soil
use grid_coms, only: nzg, nzs                       
use misc_coms, only: dtlsm
use consts_coms, only: alvi, cice, cliq, alli, t3ple, qicet3, qliqt3, tsupercool, wdns, wdnsi
use therm_lib, only: qwtk

implicit none

integer, intent(inout) :: nlev_sfcwater      ! # of active sfc water levels
integer, intent(in)    :: ntext_soil   (nzg) ! soil textural class

real, intent(out)   :: soil_rfactor    (nzg) ! soil thermal resistance [K m^2/W]
real(kind=8), intent(inout) :: soil_water      (nzg) ! soil water content [water_vol/total_vol]
real, intent(inout) :: soil_energy     (nzg) ! soil internal energy [J/m^3]
real, intent(inout) :: sfcwater_mass   (nzs) ! surface water mass [kg/m^2]
real, intent(inout) :: sfcwater_energy (nzs) ! surface water energy [J/kg]
real, intent(inout) :: sfcwater_depth  (nzs) ! surface water depth [m]
real, intent(inout) :: soil_tempk      (nzg) ! soil temperature [K]
real, intent(inout) :: soil_fracliq    (nzg) ! fraction of soil water in liq phase
real, intent(inout) :: sfcwater_tempk  (nzs) ! surface water temp [K]
real, intent(inout) :: sfcwater_fracliq(nzs) ! fraction of sfc water in liq phase
real, intent(inout) :: rshort_s        (nzs) ! s/w net rad flux to sfc water [W/m^2]
real, intent(in)    :: hxfersc ! sfc_water-to-can_air heat xfer this step [J/m^2]
real, intent(in)    :: wxfersc ! sfc_water-to-can_air vap xfer this step [kg_vap/m^2]
real, intent(in)    :: rlong_g ! l/w net rad flux to soil [W/m^2]
real, intent(in)    :: rlong_s ! l/w net rad flux to sfc water [W/m^2]
real, intent(in)    :: pcpg    ! new pcp amount this leaf timestep [kg/m^2]
real, intent(in)    :: qpcpg   ! new pcp energy this leaf timestep [J/m^2]
real, intent(in)    :: dpcpg   ! new pcp depth this leaf timestep [m]
real, intent(in)    :: wshed   ! water shed from veg this timestep [kg/m^2]
real, intent(in)    :: qwshed  ! water energy shed from veg this timestep [J/m^2]
integer, intent(in) :: lsl

! Local variables

real :: hxfers        (nzs+1) ! sfcwater heat xfer [J/m2] 
real :: sfcwater_rfactor(nzs) ! sfcwater thermal resistivity [K m^2/W]
real :: mass_new        (nzs) ! mass of new set of sfcwater layers [kg/m^2]
real :: energy_new      (nzs) ! energy of new set of sfcwater layers [J/kg]
real :: depth_new       (nzs) ! depth of new set of sfcwater layers [m]
real :: rshort_snew     (nzs) ! s/w rad flux to new set of sfcwater layers [W/m^2]
real :: energy_per_m2   (nzs) ! sfcwater energy per square meter [J/m^2]

integer :: k         ! vertical index over sfcwater layers
integer :: kold      ! vertical index of adjacent lower sfcwater layer
integer :: icelayers ! # of sfcwater layers that contain some ice
integer :: maxlayers ! max allowed # of sfcwater layers (1 if none contain ice)
integer :: nlev_new  ! new number of sfcwater layers after adjustment

real :: hxfergs   ! energy transfer from soil to sfcwater this step [J/m^2]
real :: rfac      ! bounded sfcwater thermal resistivity at k=1 [K m^2/W]
real :: snden     ! sfcwater density [kg/m^3]
real :: wfree     ! free liquid in sfcwater layer that can percolate out [kg/m^2]
real :: qwfree    ! energy carried by wfree [J/m^2]
real :: dwfree    ! depth carried by wfree [m]
real :: fracstep  ! ratio of leaf timestep to snow density exponential decay time
real :: totsnow   ! sum of mass over sfcwater layers [kg/m^2]
real :: wt        ! sfcwater(1) + soil(nzg) water masses (impl balance) [kg/m^2]
real :: qwt       ! sfcwater(1) + soil(nzg) energies (impl balance) [J/m^2]
real :: soilhcap  ! soil(nzg) heat capacity [J/(m^2 K)]
real :: soilcap   ! capacity of top soil layer to accept surface water [kg/m^2]
real :: wtnew     ! weight for new sfcwater layer when adjusting layer thickness
real :: wtold     ! weight for old sfcwater layer when adjusting layer thickness
real :: dwtold    ! change in wtold for partial mass transfer from old layer
real :: wdiff     ! change in sfcwater mass when adjusting layer thickness [kg/m^2]
real :: soilcond  ! soil thermal conductivity [W/(K m)]
real :: waterfrac ! soil water fraction in soil layer [vol_water/vol_total]
real :: tempk     ! Kelvin temperature [K]
real :: fracliq   ! fraction of water in liquid phase returned from qwtk
real :: flmin     ! lower bound on sfcwater_fracliq(1) in balance with soil
real :: flmax     ! upper bound on sfcwater_fracliq(1) in balance with soil
real :: specvol   ! specific volume of sfcwater involved in vapor xfer [m^3/kg]

! Local parameters

real, parameter :: sndenmax = wdns     ! max sfcwater density [kg/m^3]
real, parameter :: snowmin = 11.       ! min sfcwater layer mass with multiple layers [kg/m^2] 
real, parameter :: snowmin_expl = 10.  ! min sfcwater mass for explicit heat xfer [kg/m^2]
real, parameter :: rfac_snowmin = .01 ! min sfcwater rfactor [K m^2/W]

real, save, dimension(10,10) :: thick  ! snowlayer thickness scaling factor

data thick(1:10, 1)/  1., .00, .00, .00, .00, .00, .00, .00, .00, .00/
data thick(1:10, 2)/ .50, .50, .00, .00, .00, .00, .00, .00, .00, .00/
data thick(1:10, 3)/ .25, .50, .25, .00, .00, .00, .00, .00, .00, .00/
data thick(1:10, 4)/ .17, .33, .33, .17, .00, .00, .00, .00, .00, .00/
data thick(1:10, 5)/ .10, .20, .40, .20, .10, .00, .00, .00, .00, .00/
data thick(1:10, 6)/ .07, .14, .29, .29, .14, .07, .00, .00, .00, .00/
data thick(1:10, 7)/ .05, .09, .18, .36, .18, .09, .05, .00, .00, .00/
data thick(1:10, 8)/ .03, .07, .13, .27, .27, .13, .07, .03, .00, .00/
data thick(1:10, 9)/ .02, .04, .09, .18, .34, .18, .09, .04, .02, .00/
data thick(1:10,10)/ .02, .03, .06, .13, .26, .26, .13, .06, .03, .02/

! Compute soil heat resistance times HALF layer depth (soil_rfactor).

do k = lsl,nzg
   waterfrac = real(soil_water(k)) / soil(ntext_soil(k))%slmsts
   soilcond =        soil(ntext_soil(k))%soilcond0  &
      + waterfrac * (soil(ntext_soil(k))%soilcond1  &
      + waterfrac *  soil(ntext_soil(k))%soilcond2  )
   soil_rfactor(k) = dslzo2(k) / soilcond
enddo

! Check whether surface water was present at the beginning of this leaf step

if (nlev_sfcwater > 0) then

! Surface water was present.

! Compute snow heat resistance times HALF layer depth (sfcwater_rfactor).
! Sfcwater_tempk(k) should be correctly balanced value at this point, so 
! sfcwater_rfactor(k) should have correct value.  Formula applies to snow,
! so limit temperature to no greater than T3ple.

   do k = 1,nlev_sfcwater
      snden = sfcwater_mass(k) / sfcwater_depth(k)
      tempk = min(t3ple,sfcwater_tempk(k))

      sfcwater_rfactor(k) = .5 * sfcwater_depth(k)  &
         / (1.093e-3 * exp(.028 * tempk) *          &
         (.030 + snden * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))
   enddo

! Zero out sfcwater internal heat transfer array at bottom and top surfaces.
! Energy transfer at bottom and top are applied separately.

   hxfers(1) = 0.
   hxfers(nlev_sfcwater+1) = 0. 

! Compute internal sfcwater energy xfer if at least two layers exist [J/m2]

   if (nlev_sfcwater >= 2) then	
      do k = 2,nlev_sfcwater
         hxfers(k) = dtlsm * (sfcwater_tempk(k-1) - sfcwater_tempk(k))  &
                   / (sfcwater_rfactor(k-1) + sfcwater_rfactor(k))      
      enddo
   endif

! Diagnose sfcwater energy per square meter, and add contributions from 
! internal transfers of sensible heat and vapor (latent heat) and from
! shortwave radiative transfer.  This excludes energy transfer from internal
! gravitational draining of free water mass.

   do k = 1,nlev_sfcwater
      energy_per_m2(k) = sfcwater_energy(k) * sfcwater_mass(k)  &
         + hxfers(k) - hxfers(k+1) + dtlsm * rshort_s(k)
   enddo

! Evaluate conductive heat transfer between top soil layer and bottom sfcwater
! layer.  Impose minimum value on sfcwater_rfactor(1).  If minimum applies, 
! energy will be implicitly re-balanced later.  Apply heat transfer to bottom
! sfcwater layer and top soil layer.

   rfac = max(rfac_snowmin,sfcwater_rfactor(1))

   hxfergs = dtlsm * (soil_tempk(nzg) - sfcwater_tempk(1))   &
           / (soil_rfactor(nzg) + rfac)

   energy_per_m2(1) = energy_per_m2(1) + hxfergs

   soil_energy(nzg) = soil_energy(nzg) - hxfergs * dslzi(nzg)

! Apply longwave radiative transfer and sfcwater-to-can_air sensible heat
! transfer to top sfcwater layer

   energy_per_m2(nlev_sfcwater) = energy_per_m2(nlev_sfcwater)  &
      + dtlsm * rlong_s - hxfersc

endif

! If no sfcwater layers exist and there are net positive contributions to 
! sfcwater from precipitation, shedding of water from vegetation, and 
! vapor flux with canopy air, create a new sfcwater layer and initialize 
! prognostic sfcwater fields to zero.

if (nlev_sfcwater == 0 .and. wshed - wxfersc > 1.e-9) then

   nlev_sfcwater = 1

   sfcwater_mass  (1) = 0.
   sfcwater_energy(1) = 0.
   energy_per_m2  (1) = 0.
   sfcwater_depth (1) = 0.
endif

! Return if no sfcwater layers now exist

if (nlev_sfcwater < 1) return

! Sfcwater layers do exist

! Apply mass, energy, and depth contributions to top sfcwater layer from
! precipitation, shedding of water from vegetation, and vapor flux with 
! canopy air.  Get value for specific volume of sfcwater involved in vapor xfer.

specvol = wdnsi
if (wxfersc > 0.) specvol =  &
   sfcwater_depth(nlev_sfcwater) / sfcwater_mass(nlev_sfcwater)

sfcwater_mass(nlev_sfcwater) = sfcwater_mass(nlev_sfcwater)  &
   + wshed                                     &
   - wxfersc

energy_per_m2(nlev_sfcwater) = energy_per_m2(nlev_sfcwater)  &
   + qwshed                                  &
   - wxfersc * alvi 

sfcwater_depth(nlev_sfcwater) = sfcwater_depth(nlev_sfcwater)  &
   + wshed * wdnsi                              &
   - wxfersc * specvol

! Prepare to transfer water downward through snow layers by percolation.
! Fracliq is the fraction of liquid in the snowcover or surface water.
! wfree [kg/m^2] is the quantity of that liquid that is free (not attached to
! snowcover) and therefore available to drain into the layer below.
! soilcap [kg/m^2] is the amount of water that can fit in the unfilled pore
! space of the top soil layer.  wfree in the lowest sfcwater layer is limited
! by this value.

! First, prepare to sum sfcwater mass over all existing layers

totsnow = 0.0

! Loop downward through all existing sfcwater layers beginning with top layer

do k = nlev_sfcwater,1,-1

! Diagnose sfcwater density.

   snden = sfcwater_mass(k) / sfcwater_depth(k)

! Assume that as snow ages on ground, its density difference with a limiting
! maximum value (currently set to 400 kg/m^3) decays exponentially (with a
! decay time currently set to about 3 weeks).  If sfcwater density is less
! than this limiting value, apply the density increase for this timestep.

! This formulation and decay constants are very crude approximations to a few
! widely variable observations of snowpack density and are intended only as a 
! rough representation of the tendency for snowcover to compress with time.  
! A better formulation that accounts for several environmental factors would
! be desirable here.

   if (snden < 400.) then
      fracstep = .5e-6 * dtlsm  ! .5e-6 is inverse decay time scale
      snden = snden * (1. - fracstep) + 400. * fracstep
      sfcwater_depth(k) = sfcwater_mass(k) / snden
   endif

! Diagnose sfcwater temperature and liquid fraction now that new mass and
! energy contributions have been applied.  Use qwtk instead of qtk in case
! sfcwater_mass(k) is too small; "dryhcap" = 100 is very small value    

   call qwtk(energy_per_m2(k),sfcwater_mass(k),100.,  &
             sfcwater_tempk(k),sfcwater_fracliq(k))

! If this is bottom layer, diagnose sfcwater_rfactor.  Since sfcwater_tempk(k)
! may not be a stable computation at this point, assume that tempk = 0, which
! gives the minimum rfactor for a given density.

   if (k == 1) then
      tempk = t3ple

      sfcwater_rfactor(k) = .5 * sfcwater_depth(k)  &
         / (1.093e-3 * exp(.028 * tempk) *          &
         (.030 + snden * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))

   endif

! If this is bottom layer and sfcwater rfactor is less than minimum stable 
! value, bring bottom sfcwater and top soil layer into thermal equilibrium 
! by exchanging heat between them.

   if (k == 1 .and.   &
       (sfcwater_mass(1) < snowmin_expl .or.  &
        sfcwater_rfactor(1) < rfac_snowmin)) then

! Combined sfcwater and soil water mass per square meter

      wt = sfcwater_mass(1) + sngl(soil_water(nzg)) * 1.e3 * dslz(nzg)

! Combined sfcwater and soil energy per square meter.

      qwt = energy_per_m2(1) + soil_energy(nzg) * dslz(nzg)

! Soil heat capacity per square meter

      soilhcap = soil(ntext_soil(nzg))%slcpd * dslz(nzg)

! Diagnose equilibrium temperature and fractional liquid/ice water phases

      call qwtk(qwt,wt,soilhcap,tempk,fracliq)

! Diagnose new energy value for sfcwater based on qwt value.

      if (qwt < wt*qicet3) then
      
! Case of equilibrium temperature below 0 deg C.  Sfcwater fracliq = 0.

         sfcwater_fracliq(1) = 0.
         sfcwater_tempk(1) = tempk
         energy_per_m2(1) = sfcwater_mass(1) * cice * tempk 

      elseif (qwt > wt * qliqt3) then
      
! Case of equilibrium temperature above 0 deg C.  Sfcwater fracliq = 1.

         sfcwater_fracliq(1) = 1.
         sfcwater_tempk(1) = tempk
         energy_per_m2(1) = sfcwater_mass(1) * cliq * (tempk - tsupercool)

      else
      
! Equilibrium temperature is 0 deg C.  Determine separate values for
! sfcwater_fracliq(1) and soil_fracliq(nzg) using constraint that the sum
! of (mass * fracliq) over both components is (wt * fracliq).

! Lower bound on sfcwater_fracliq(1): case with soil_water(nzg) all liquid

         flmin = (fracliq * wt - sngl(soil_water(nzg)) * wdns * dslz(nzg))  &
               / sfcwater_mass(1)         

! Upper bound on sfcwater_fracliq(1): case with soil_water(nzg) all ice

         flmax = fracliq * wt / sfcwater_mass(1)         

! New sfcwater_fracliq(1) value becomes closest value within bounds to old value.

         sfcwater_fracliq(1) = max(0.,flmin,min(1.0,flmax,sfcwater_fracliq(1)))
         sfcwater_tempk(1) = t3ple
         energy_per_m2(1) = sfcwater_mass(1) * sfcwater_fracliq(1) * qliqt3

      endif

! New energy value for soil is combined energy minus new sfcwater energy

      soil_energy(nzg) = (qwt - energy_per_m2(1)) * dslzi(nzg)

   else

! Current sfcwater layer is either not the bottom one or is thick enough
! to not require implicit thermal balance with top soil layer.

! Diagnose sfcwater temperature and liquid fraction.  Use qwtk instead of qtk
! in case sfcwater_mass(k) is too small; "dryhcap" = 100 is neglible value    
  
      call qwtk(energy_per_m2(k),sfcwater_mass(k),100.,  &
                sfcwater_tempk(k),sfcwater_fracliq(k))

   endif

! If liquid exists in current sfcwater layer, any low-density ice structure
! tends to collapse.  Increase density accordingly using simple linear relation.

   if (snden < 1.e3 * sfcwater_fracliq(k)) then
      snden = 1.e3 * sfcwater_fracliq(k)
      sfcwater_depth(k) = sfcwater_mass(k) / snden
   endif

! Assume that excess of sfcwater_fracliq over 10% is free to drain out of layer

   if (sfcwater_fracliq(k) > .10) then
      wfree = sfcwater_mass(k) * (sfcwater_fracliq(k) - .10) / .90

! Check whether this is lowest sfcwater layer

      if (k == 1) then

! This is lowest sfcwater layer.  Reduce wfree if necessary to maximum 
! amount that can percolate into top soil layer

         soilcap = 1.e3 * max (0.,dslz(nzg)  &
                 * (soil(ntext_soil(nzg))%slmsts - sngl(soil_water(nzg))))
         if (wfree > soilcap) wfree = soilcap

      endif

! Evaluate energy and depth transferred in wfree (which is in liquid phase)

      qwfree = wfree * cliq * (sfcwater_tempk(k) - tsupercool)
      dwfree = wfree * wdnsi

! Check if essentially all of sfcwater_mass(k) will drain from layer

      if (wfree > .999 * sfcwater_mass(k)) then
      
! All sfcwater_mass(k) drains from layer.  Set layer quantities to zero to
! avoid truncation error.

         sfcwater_mass(k)  = 0.
         energy_per_m2(k)  = 0.
         sfcwater_depth(k) = 0.
          
      else

! Not all sfcwater_mass(k) drains from layer.  Drain mass, energy, and depth 
! of free water out of current layer

         sfcwater_mass(k)  = sfcwater_mass(k)  - wfree
         energy_per_m2(k)  = energy_per_m2(k)  - qwfree
         sfcwater_depth(k) = sfcwater_depth(k) - dwfree

      endif

! Check whether this is lowest sfcwater layer

      if (k == 1) then

! This is lowest sfcwater layer.  Drained water goes to top soil layer

         soil_water(nzg) = soil_water(nzg) + 1.d-3 * dble(wfree) * dble(dslzi(nzg))
         soil_energy(nzg) = soil_energy(nzg) + qwfree * dslzi(nzg)

      else

! This is not lowest sfcwater layer.  Drained water goes to next lower
! sfcwater layer

         sfcwater_mass(k-1)  = sfcwater_mass(k-1)  + wfree
         energy_per_m2(k-1)  = energy_per_m2(k-1)  + qwfree
         sfcwater_depth(k-1) = sfcwater_depth(k-1) + dwfree

      endif

   endif

! Add remaining sfcwater mass in current layer to mass sum

   totsnow = totsnow + sfcwater_mass(k)

enddo

! Check whether any sfcwater mass remains (after possibly having all drained into soil)

if (totsnow < 1.e-9) then

! Total sfcwater mass is very close to zero.  Set sfcwater layer count to zero,
! set sfcwater arrays to exactly zero, and return

   nlev_sfcwater = 0

   sfcwater_mass(1:nzs)   = 0.
   sfcwater_energy(1:nzs) = 0.
   sfcwater_depth(1:nzs)  = 0.

   return
   
endif

! Sfcwater mass remains.  Re-distribute mass, energy, and depth among layers to
! maintain prescribed distribution of mass

! Find maximum number of layers for which thinnest layer (top and bottom) would 
! be thicker than snowmin.

maxlayers = 1

do while (maxlayers < nzs .and. maxlayers < 10 .and.  &
   thick(1,maxlayers+1) * totsnow > snowmin)

   maxlayers = maxlayers + 1
enddo

! Count up all existing snow layers that are not totally liquid

icelayers = 0

do k = 1,nlev_sfcwater
   if (sfcwater_mass(k) > 1.e-9 .and.  &
       energy_per_m2(k) < sfcwater_mass(k) * alli) then

      icelayers = icelayers + 1
   endif
enddo

! Determine new number of layers.  This number may be at most one greater than
! nlev_sfcwater and one greater than icelayers, but no greater than maxlayers.

nlev_new = min (nlev_sfcwater+1, icelayers+1, maxlayers)

! Set index for first old layer; set transfer weights for first old layer

kold = 1
wtold = 1. ! fraction of "remaining" mass in old layer

! Loop over new set of layers

do k = 1,nlev_new

! To begin current new layer, set fraction of "unfilled" mass in new layer to 1.

   wtnew = 1.

! Set mass of current new layer (already determined)

   mass_new(k) = totsnow * thick(k,nlev_new)

! Initialize energy, depth, and s/w flux of current new layer to zero

   energy_new (k) = 0.
   depth_new  (k) = 0.
   rshort_snew(k) = 0.

10 continue

! Compute "unfilled" mass in current new layer minus remaining mass in
! current old layer

   wdiff = wtnew * mass_new(k) - wtold * sfcwater_mass(kold)

! Check sign of wdiff

   if (wdiff > 0.) then

! If "unfilled" mass in current new layer exceeds remaining mass in current old
! layer, transfer all of remaining energy, depth, and s/w flux from old layer

      energy_new (k) = energy_new (k) + wtold * energy_per_m2 (kold)
      depth_new  (k) = depth_new  (k) + wtold * sfcwater_depth(kold)
      rshort_snew(k) = rshort_snew(k) + wtold * rshort_s      (kold)

! Reduce fraction of "unfilled" mass in current new layer, jump to next old
! layer, and return wtold to 1.0 to indicate no mass in old layer yet removed.

      wtnew = wtnew - wtold * sfcwater_mass(kold) / mass_new(k)
      kold = kold + 1
      wtold = 1.

! If old-layer counter does not exceed top old layer, repeat transfer operation 
! for current old layer

      if (kold <= nlev_sfcwater) go to 10

   else

! If "unfilled" mass in current new layer is less than remaining mass in current
! old layer, transfer only the portion of remaining energy, depth, and s/w flux
! from old layer that fit into new layer.

! Change in wtold

      dwtold = wtnew * mass_new(k) / sfcwater_mass(kold)

      wtold = wtold - dwtold

! Energy, depth, and s/w flux transfer to current new layer

      energy_new (k) = energy_new (k) + dwtold * energy_per_m2 (kold)
      depth_new  (k) = depth_new  (k) + dwtold * sfcwater_depth(kold)
      rshort_snew(k) = rshort_snew(k) + dwtold * rshort_s      (kold)

   endif

enddo

! Now that mass, energy, depth, and s/w flux have been transferred to new layers,
! copy back to original arrays

do k = 1,nlev_new
   sfcwater_mass(k)   = mass_new(k)
   sfcwater_energy(k) = energy_new(k) / sfcwater_mass(k)
   sfcwater_depth(k)  = depth_new(k)
   rshort_s(k)        = rshort_snew(k)

! Replace sfcwater energy limiter from earlier code versions with energy
! check and message.  If message ever gets printed, investigate reasons.

   if (sfcwater_energy(k) > 6.e5 .or. sfcwater_energy(k) < -2.e5) then
      print*, ' '
      print*, 'Sfcwater energy is outside allowable range. '
      print*, 'k,sfcwater_energy = ',k,sfcwater_energy(k)
      print*, 'p1',energy_new(k),sfcwater_mass(k),nlev_sfcwater,nlev_new
      print*, 'p2',kold, energy_per_m2(kold)
      stop 'stop sfcwater energy'
   endif

enddo

nlev_sfcwater = nlev_new

return
end subroutine sfcwater
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine soil_euler_ar(nlev_sfcwater, ntext_soil, ktrans,    &
                soil_tempk, soil_fracliq, soil_rfactor,                  &
                hxfergc, wxfergc, rshort_g, rlong_g, transp,             &
                soil_energy, soil_water, hxferg, wxfer, qwxfer,          &
                psiplusz, hydraul_cond, lsl, ed_transp, csite,ipa         )

use soil_coms, only: dslz, dslzi, slzt, dslzidt, dslztidt, soil, slcons1 
use grid_coms, only: nzg
use consts_coms, only: cliqvlme, allivlme, alvi,t3ple, tsupercool

use ed_state_vars,only:sitetype,patchtype

use misc_coms, only: dtlsm

implicit none

integer, intent(in) :: nlev_sfcwater    ! # active levels of surface water
integer, intent(in) :: ntext_soil (nzg) ! soil textural class

integer, intent(in) :: ktrans           ! k index of soil layer supplying transpiration

real, intent(in) :: soil_tempk    (nzg) ! soil temperature (K)
real, intent(in) :: soil_fracliq  (nzg) ! fraction of soil water that is liquid
real, intent(in) :: soil_rfactor  (nzg) ! soil thermal resistance
real, intent(in) :: hxfergc             ! heat xfer from ground to canopy [kg/m^2]
real, intent(in) :: wxfergc             ! water xfer from ground to canopy [kg/m^2]
real, intent(in) :: rshort_g            ! s/w radiative flux abs by ground [W/m^2]
real, intent(in) :: rlong_g             ! l/w radiative flux abs by ground [W/m^2]
real, intent(in) :: transp              ! transpiration loss [kg/m^2]

real, intent(inout) :: soil_energy(nzg) ! [J/m^3]
real(kind=8), intent(inout) :: soil_water (nzg) ! soil water content [vol_water/vol_tot]

real, intent(out) :: hxferg      (nzg+1) ! soil internal heat xfer (J/m^2]
real, intent(out) :: wxfer       (nzg+1) ! soil water xfer [m]
real, intent(out) :: qwxfer      (nzg+1) ! soil energy xfer from water xfer [J/m^2] 
real, intent(out) :: psiplusz    (nzg)   ! soil water potential (including grav) [m]
real, intent(out) :: hydraul_cond(nzg)   ! soil hydraulic conductivity [m/s]

type(sitetype),target :: csite
integer :: ipa

integer, intent(in) :: lsl
! Local variables

real :: half_soilair(nzg) ! half of available airspace in soil layer [m]
real :: soil_liq    (nzg) ! liq water in soil layer available for transport [m]

integer :: k     ! vertical index over soil layers
integer :: nts   ! soil textural class

real :: watermid ! soil water content midway between layers [vol_water/vol_tot]
real :: wloss    ! soil water loss from transpiration [vol_water/vol_tot]
real :: qwloss   ! soil energy loss from transpiration [J/vol_tot]
real, dimension(nzg) :: ed_transp

! Remove transpiration water from ktrans soil layer
! Units of wloss are [vol_water/vol_tot], of transp are [kg/m^2].

do k = lsl, nzg
   wloss = ed_transp(k) * dslzi(k) * 1.e-3
   qwloss = wloss * cliqvlme * (soil_tempk(k) - tsupercool)
   soil_water(k) = soil_water(k) - dble(wloss)
   soil_energy(k) = soil_energy(k) - qwloss
   csite%mean_latflux(ipa) =   &
        csite%mean_latflux(ipa) + qwloss * dslz(k) / dtlsm

enddo

! Soil bottom, top, and internal sensible heat xfers [J/m2]

hxferg(lsl) = 0.
hxferg(nzg+1) = 0.

do k = lsl+1,nzg
   hxferg(k) = dtlsm * (soil_tempk(k-1) - soil_tempk(k))   &
             / (soil_rfactor(k-1) + soil_rfactor(k))      
enddo

! Update soil Q values [J/m3] from internal sensible heat and upward water 
! vapor (latent heat) xfers, and from longwave and shortwave fluxes at the
! top of the soil.  This excludes effects of dew/frost formation, 
! precipitation, shedding, and percolation, which were already applied 
! to the top soil layer in subroutine sfcwater.  Update top soil moisture 
! from evaporation only if sfcwater was absent.

do k = lsl,nzg
   soil_energy(k) = soil_energy(k) + dslzi(k) * (hxferg(k) - hxferg(k+1))
enddo

soil_energy(nzg) = soil_energy(nzg)  &
   + dslzi(nzg) * (dtlsm * (rshort_g + rlong_g) - hxfergc - wxfergc * alvi)

if (nlev_sfcwater == 0) then
   soil_water(nzg) = soil_water(nzg) - 1.d-3 * dble(wxfergc) * dble(dslzi(nzg))
endif

! [12/07/04] Revisit the computation of water xfer between soil layers in
! relation to the extreme nonlinearity of hydraulic conductivity with respect
! to soil moisture.  What is the best value for hydraulic conductivity (or
! resistivity) at the interface between the two layers?  The answer is
! definitely not the average of the resistivity values of the layers
! because a very dry layer would shut down xfer of water into it.  The average 
! conductivity would, on the other hand, over-estimate water xfer between layers 
! when one is wet and the other dry.  A good compromise seems to be to average
! the fractional moisture content between the two layers and to apply this
! average value in computing hydraulic resistance for the bottom half of the
! upper layer and the top half of the lower layer.  Then, add these resistances
! to obtain the total hydraulic resistance between the two layers.

! Compute gravitational potential plus moisture potential z + psi [m],
! liquid water content [m], and half the remaining water capacity [m].

do k = lsl,nzg
   nts = ntext_soil(k)

   psiplusz(k) = soil(nts)%slpots * (soil(nts)%slmsts / sngl(soil_water(k))) ** soil(nts)%slbs  &
               + slzt(k) 

   soil_liq(k) = dslz(k)  &
      * sngl(min(soil_water(k) - dble(soil(nts)%soilcp) , soil_water(k) * dble(soil_fracliq(k)) ) )

   half_soilair(k) = (soil(nts)%slmsts - sngl(soil_water(k))) * dslz(k) * .5
enddo

! Find amount of water transferred between soil layers (wxfer) [m]
! modulated by the liquid water fraction

wxfer(lsl)      = 0.
wxfer(nzg+1)  = 0.
qwxfer(lsl)     = 0.
qwxfer(nzg+1) = 0.

do k = lsl+1,nzg
   nts = ntext_soil(k)

   watermid = .5 * sngl(soil_water(k) + soil_water(k-1))
   
   hydraul_cond(k) = slcons1(k,nts)  &
      * (watermid / soil(nts)%slmsts) ** (2. * soil(nts)%slbs + 3.)
      
   wxfer(k) = dslztidt(k) * hydraul_cond(k) * (psiplusz(k-1) - psiplusz(k)) &
      * .5 * (soil_fracliq(k) + soil_fracliq(k-1))

! Limit water transfers to prevent over-saturation and over-depletion
! Compute q transfers between soil layers (qwxfer) [J/m2]

   if (wxfer(k) > 0.) then
      wxfer(k) = min(wxfer(k),soil_liq(k-1),half_soilair(k))
   else
      wxfer(k) = - min(-wxfer(k),soil_liq(k),half_soilair(k-1))
   endif

   qwxfer(k) = wxfer(k) * cliqvlme * (soil_tempk(k) - tsupercool)

enddo

! Update soil moisture (impose minimum value of soilcp) and q value.

do k = lsl,nzg
   nts = ntext_soil(k)
   soil_water(k) = max(dble(soil(nts)%soilcp),soil_water(k)  &
      - dble(dslzi(k)) * (dble(wxfer(k+1)) - dble(wxfer(k))))
   soil_energy(k) = soil_energy(k) - dslzi(k) * (qwxfer(k+1) - qwxfer(k))
enddo

! Compute soil respiration
call soil_respiration_ar(csite,ipa)

return
end subroutine soil_euler_ar

!*****************************************************************************

!*****************************************************************************

subroutine remove_runoff(ksn, sfcwater_fracliq, sfcwater_mass,   &
     sfcwater_tempk, sfcwater_energy, sfcwater_depth, runoff, qrunoff)

  use soil_coms, only: runoff_time
  use grid_coms, only: nzs
  use consts_coms, only: alli, cliq,t3ple, wdnsi, tsupercool
  use misc_coms, only: dtlsm
  use therm_lib, only: qtk

  implicit none

  integer, intent(in) :: ksn

  real, intent(inout)    :: sfcwater_fracliq(nzs)
  real, intent(inout)    :: sfcwater_tempk  (nzs)
  real, intent(inout) :: sfcwater_mass   (nzs)
  real, intent(inout) :: sfcwater_energy (nzs)
  real, intent(inout) :: sfcwater_depth  (nzs)
  real, intent(out) :: runoff
  real, intent(out) :: qrunoff

  ! Get rid of runoff

  runoff = 0.0
  qrunoff = 0.0

  if(ksn >= 1)then
     if(sfcwater_fracliq(ksn) > 0.1)then
        call qtk(sfcwater_energy(ksn), sfcwater_tempk(ksn),   &
             sfcwater_fracliq(ksn))
        runoff = sfcwater_mass(ksn) * (sfcwater_fracliq(ksn) - 0.1) /   &
             (0.9 * runoff_time) * dtlsm
        qrunoff = runoff * cliq * (sfcwater_tempk(ksn) - tsupercool)
        
        sfcwater_energy(ksn) = (sfcwater_energy(ksn) *   &
             sfcwater_mass(ksn) - qrunoff ) / (sfcwater_mass(ksn) - runoff)
        sfcwater_mass(ksn) = sfcwater_mass(ksn) - runoff
        sfcwater_depth(ksn) = sfcwater_depth(ksn) - wdnsi * runoff
     endif
  endif
  
  return
end subroutine remove_runoff




