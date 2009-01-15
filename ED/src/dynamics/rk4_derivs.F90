subroutine leaf_derivs_ar(initp, dinitp, csite, ipa,isi,ipy, rhos, prss, &
      pcpg, qpcpg, atm_tmp, exner, geoht, vels, atm_shv, atm_co2, lsl)
  
  use ed_state_vars,only:sitetype,rk4patchtype
  use consts_coms, only : cp
  use grid_coms, only: nzg

  implicit none

  type(sitetype),target :: csite
  integer :: ipa,isi,ipy

  integer, intent(in) :: lsl
  type(rk4patchtype), target :: initp
  type(rk4patchtype), target :: dinitp

  real, intent(in) :: rhos
  real, intent(in) :: prss
  real, intent(in) :: pcpg
  real, intent(in) :: qpcpg
  real, intent(in) :: atm_tmp
  real, intent(in) :: exner
  real, intent(in) :: geoht
  real, intent(in) :: vels
  real, intent(in) :: atm_shv
  real, intent(in) :: atm_co2

#if USE_INTERF
  interface
    subroutine ed_stars(tha,rva,chia,thv,zpm,um,rough,ustar,rstar,tstar,  &
         cstar,can_shv,can_co2)
      implicit none
      real, intent(in) :: rough
      real, intent(out) :: ustar
      real, intent(out) :: rstar
      real, intent(out) :: tstar
      real, intent(out) :: cstar
      real, intent(in) :: can_shv
      real, intent(in) :: can_co2
      real ::tha,rva,chia,thv,zpm,um
    end subroutine ed_stars

    subroutine leaftw_derivs_ar(initp, dinitp, csite, ipa,isi,ipy, rhos, prss, pcpg,   &
         qpcpg, atm_tmp, exner, geoht, lsl)
      
      use ed_state_vars,only:rk4patchtype,sitetype
      
      implicit none
      type(sitetype),target :: csite
      integer :: ipa,isi,ipy
      integer, intent(in) :: lsl
      type (rk4patchtype) ,target :: initp
      type (rk4patchtype) ,target :: dinitp
      real, intent(in) :: rhos
      real, intent(in) :: prss
      real, intent(in) :: pcpg
      real, intent(in) :: qpcpg
      real, intent(in) :: atm_tmp
      real, intent(in) :: exner
      real, intent(in) :: geoht
    end subroutine leaftw_derivs_ar
  end interface
#endif

  dinitp%ebudget_latent = 0.0
  ! Compute friction velocities
  call ed_stars(atm_tmp, atm_shv, atm_co2, initp%can_temp, geoht, vels, &
       initp%rough, initp%ustar, initp%rstar, initp%tstar, initp%cstar,  &
       initp%can_shv, initp%can_co2)
  initp%tstar = initp%tstar * cp / exner
  call leaftw_derivs_ar(initp, dinitp, csite,ipa,isi,ipy, rhos, prss, pcpg, qpcpg, &
       atm_tmp, exner, geoht, lsl)
  dinitp%nlev_sfcwater = initp%nlev_sfcwater

  return
end subroutine leaf_derivs_ar

!==================================================================

subroutine leaftw_derivs_ar(initp, dinitp, csite,ipa,isi,ipy, rhos, prss, pcpg, qpcpg,   &
     atm_tmp, exner, geoht, lsl)

  use max_dims, only : nzgmax,nzsmax
  use consts_coms, only : alvl, cliq1000, cpi, alvi, alli1000, t3ple
  use grid_coms, only: nzg,nzs
  use soil_coms, only : soil, slz, dslz, dslzi, water_stab_thresh,  &
       infiltration_method, dslzti, slcons1, slzt, min_sfcwater_mass,ss,isoilbc
  use misc_coms, only: current_time
  use canopy_radiation_coms, only: lai_min

  use ed_state_vars,only:sitetype,patchtype,rk4patchtype
  use ed_therm_lib,only: ed_grndvap
  use therm_lib, only : qtk, qwtk, qwtk8

  implicit none

  type(sitetype),target :: csite
  integer :: ipa,isi,ipy
  integer, intent(in) :: lsl
  type(rk4patchtype), target :: initp,dinitp
  
  real, intent(in) :: pcpg
  real, intent(in) :: qpcpg
  real, intent(in) :: exner
  real, intent(in) :: atm_tmp
  real, intent(in) :: rhos
  real, intent(in) :: prss
  real, intent(in) :: geoht
  integer :: k,ksn,nsoil
  real :: wgpfrac,soilcond,snden,hflxgc,wflxgc,dewgnd,wshed
  real, dimension(nzgmax+nzsmax)   :: rfactor
  real, dimension(nzgmax+nzsmax+1) :: hfluxgsc
  real, dimension(nzgmax)             :: soil_liq,psiplusz,soilair99
  real, dimension(nzg+nzs+1) :: w_flux,qw_flux
  real, dimension(nzs+1) :: d_flux
  real :: qwshed
  
  !----------
  ! This is the free surface water transfer inverse time.  I made this up.  It 
  ! requires testing.
  !  real, parameter :: taui = 1.0/30.0 !!!!  Standard
  real, parameter :: taui = 1.0/300.0
  !----------
  ! This is the exponent in the frozen soil hydraulic conductivity correction
  real, parameter :: freezeCoef = 7.0

  real :: wgpmid,wloss
  real :: dqwt

  real :: fracliq,tempk
  integer :: k1,k2
  real :: qwloss
  real :: surface_water !! pool of availible liquid water on the soil surface (kg/m2?)
  real :: infilt !! surface infiltration rate
  real :: snowdens !! snow density (kg/m2)
  real :: freezeCor !! correction to conductivity for partially frozen soil

  logical, parameter :: debug = .false.

#if USE_INTERF
  interface
     subroutine canopy_derivs_two_ar(initp, dinitp, csite,ipa,isi,ipy, hflxgc, wflxgc,   &
         dewgndflx, wshed_tot, qwshed_tot, rhos, prss, pcpg, qpcpg, &
         exner, geoht, atm_tmp, lsl)

       use ed_state_vars,only: rk4patchtype,sitetype,patchtype

       implicit none
       type(sitetype),target :: csite
       integer :: ipa,isi,ipy
       integer, intent(in) :: lsl
       real, intent(in) :: rhos
       real, intent(in) :: atm_tmp
       real, intent(in) :: prss
       real, intent(in) :: pcpg
       real, intent(in) :: qpcpg
       real, intent(in) :: exner
       real, intent(in) :: geoht
       type (rk4patchtype) ,target :: initp
       type (rk4patchtype) ,target :: dinitp
       real, intent(out) :: hflxgc,wflxgc,dewgndflx,wshed_tot,qwshed_tot
     end subroutine canopy_derivs_two_ar
  end interface
#endif

  ! Variables of convenience
  ksn = initp%nlev_sfcwater
  
  w_flux = 0.0
  qw_flux = 0.0
  d_flux = 0.0

  ! Initialize derivatives to zero
  ! initp contains the current values of the state variables of the patch
  ! dinitp contains the derivative of the state variable in the patch
  dinitp%soil_energy(:) = 0.0
  dinitp%soil_water(:) = 0.0d+0
  dinitp%sfcwater_depth(:) = 0.0
  dinitp%sfcwater_energy(:) = 0.0  ! THIS IS IN W/m2!!!  Convert at end.
  dinitp%sfcwater_mass(:) = 0.0
  dinitp%virtual_heat = 0.0
  dinitp%virtual_water = 0.0
  dinitp%virtual_depth = 0.0
  initp%extracted_water(:) = 0.0
  dinitp%avg_smoist_gc(:) = 0.0

  nsoil = csite%ntext_soil(nzg,ipa)
  call ed_grndvap(ksn,                                &
       nsoil,  &
       initp%soil_water       (nzg),  &
       initp%soil_energy      (nzg),  &
       initp%sfcwater_energy(max(1,ksn)), &
       rhos,  &
       initp%can_shv,  &
       initp%ground_shv,  &
       initp%surface_ssh)

  ! Calculate water available to vegetation (in meters)
  ! SLZ is specified in RAMSIN.  Each element of the array sets the value 
  ! of the bottom of a corresponding soil layer.
  ! Eg, SLZ = -2, -1, -0.5, -0.25.  There are four soil layers in this example;
  ! soil layer 1 goes from 2 meters below the 
  ! surface to 1 meter below the surface.
  ! ---------------------------------------------------------------------------

  initp%available_liquid_water(nzg) = dslz(nzg) * max(0.0,  &
       &  initp%soil_fracliq(nzg) * (sngl(initp%soil_water(nzg)) - soil(nsoil)%soilcp))

  do k = nzg - 1, lsl, -1
     nsoil = csite%ntext_soil(k,ipa)
     initp%available_liquid_water(k) = initp%available_liquid_water(k+1) +  &
          dslz(k) * max(0.0, (sngl(initp%soil_water(k)) - soil(nsoil)%soilcp) *  &
          initp%soil_fracliq(k))

  enddo

  ! Get derivatives of canopy variables
  call canopy_derivs_two_ar(initp, dinitp, csite, ipa,isi,ipy,hflxgc, wflxgc, dewgnd,  &
       wshed,qwshed, rhos, prss, pcpg, qpcpg, exner, geoht, atm_tmp, lsl)

  
  ! Calculate conductivities:
  !        soil


  ! We should not have any nsoil checks, instead bedrock should
  ! be given a zero conductivity parameter
  


  do k = lsl, nzg
     nsoil = csite%ntext_soil(k,ipa)
     if(nsoil <= 12)then
        wgpfrac = min(sngl(initp%soil_water(k)) / soil(nsoil)%slmsts,1.0)
        soilcond = soil(nsoil)%soilcond0 + wgpfrac * (soil(nsoil)%soilcond1  &
             + wgpfrac * soil(nsoil)%soilcond2)
     else
        soilcond=soil(nsoil)%soilcond0
     endif
     rfactor(k) = dslz(k) / soilcond
  enddo


  !        snow/surface water
  do k = 1, ksn
     if(initp%sfcwater_depth(k) > 0.0)then
        snden = initp%sfcwater_mass(k) / initp%sfcwater_depth(k)
        rfactor(k+nzg) = initp%sfcwater_depth(k)  &
             / (ss(1) * exp(ss(2) * initp%sfcwater_tempk(k))   &
             * (ss(3) + snden * (ss(4) + snden * (ss(5) + snden * ss(6)))))
     else
        rfactor(k+nzg) = 0.0
     endif
  enddo

  ! Calculate the sensible heat fluxes
  !        find soil and sfcwater internal sensible heat fluxes (hfluxgsc) [W/m2]
  !MLO - Uncommented here, hfluxgsc(1) is used to determine the soil heat later on...
  hfluxgsc(:) = 0.0
  do k = lsl+1, nzg
     hfluxgsc(k) = - (initp%soil_tempk(k) - initp%soil_tempk(k-1))   &
                   / ((rfactor(k)   + rfactor(k-1)) * .5)
     
     dinitp%avg_sensible_gg(k-1) = hfluxgsc(k)  !Diagnostic

  enddo

  if(ksn >= 1)then

     hfluxgsc(nzg+1) = - (initp%sfcwater_tempk(1) - initp%soil_tempk(nzg)) &
          / ((rfactor(nzg+1)   + rfactor(nzg)) * .5)

     do k = 2,ksn
        hfluxgsc(nzg+k) = - (initp%sfcwater_tempk(k) -   &
             initp%sfcwater_tempk(k-1)) / &
             ((rfactor(nzg+k) + rfactor(nzg+k-1)) * .5)
     enddo
     
     ! We will need the liquid fraction of the surface water 
     ! to partition the latent heat of sublimation
     ! and evaporation
     
     call qtk(initp%sfcwater_energy(ksn),tempk,fracliq)
     
  else

     ! We will need the liquid fraction of the soil water
     ! to partition the latent heat of sublimation
     ! and evaporation
     
     call qwtk8(initp%soil_energy(nzg),initp%soil_water(nzg)*1.d3, &
          soil(nsoil)%slcpd,tempk,fracliq)

  endif

  !      heat flux (hfluxgsc) at soil or sfcwater top from longwave, sensible
  !      [W/m^2]
  hfluxgsc(nzg+ksn+1) = hflxgc + wflxgc * (fracliq*alvl+(1.0-fracliq)*alvi) - csite%rlong_g(ipa) - csite%rlong_s(ipa)

  
  dinitp%avg_sensible_gg(nzg)=hfluxgsc(nzg+ksn+1) ! Diagnostic

  ! Update soil Q values [J/m3] from sensible heat, upward water vapor 
  ! (latent heat) and longwave fluxes. This excludes effects of dew/frost 
  ! formation, precipitation, shedding, and percolation

  do k = lsl,nzg
     dinitp%soil_energy(k) = dslzi(k) * (hfluxgsc(k)- hfluxgsc(k+1))
  enddo

  ! Update soil Q values [J/m3] from shortwave flux.
  dinitp%soil_energy(nzg) = dinitp%soil_energy(nzg) + dslzi(nzg) * csite%rshort_g(ipa)

  ! Update sfcwater Q values [J/kg] from sensible heat, upward water vapor 
  ! (latent heat), longwave, and shortwave fluxes.  This excludes effects 
  ! of dew/frost formation, precipitation, shedding and percolation.  

  do k = 1,ksn
     dinitp%sfcwater_energy(k) = hfluxgsc(k+nzg) - hfluxgsc(k+1+nzg) +   &
          csite%rshort_s(k,ipa)
  enddo

  ! Calculate the fluxes of water with their associated heat fluxes.
  ! Update top soil or snow moisture from evaporation only.

  ! New moisture, qw, and depth from dew/frost formation, precipitation,
  ! shedding, and percolation.  ksnnew is the layer that receives the new
  ! condensate that comes directly from the air above.  If there is no
  ! pre-existing snowcover, this is a temporary "snow" layer.  This weird
  ! factor is essentially the one used by default in RAMS 4.3.0.  I'd rather
  ! use a different factor, huh?

  qw_flux(nzg+ksn+1) = - dewgnd * (fracliq*alvl+(1.0-fracliq)*alvi) - qwshed
  w_flux(nzg+ksn+1) = - dewgnd - wshed
  d_flux(ksn+1) = w_flux(nzg+ksn+1) * 0.001
  
  !! account for SNOW DENSITY
  !! fcn derived from CLM3.0 documentation which is based on
  !! Anderson 1975 NWS Technical Doc # 19 
  !! which I have yet to find   <mcd>
  if(w_flux(nzg+1) < 0.0)then ! gaining water

     ! DMM: Note that qw_flux(nzg+ksn+1) > 0 represents a flux of frozen water.
     
!     call qtk(qw_flux(nzg+ksn+1)/w_flux(nzg+ksn+1),tempk,fracliq)
!     if(fracliq < 0.9) then !snow is falling, make light and fluffy
     if(qw_flux(nzg+ksn+1) > 0.0) then !snow is falling, make light and fluffy
        snowdens = 50.0 
        tempk = atm_tmp !! set temperature to atm
        if(tempk > 275.15)tempk=275.15
        if(tempk > 258.15)snowdens=50.0+1.5*(tempk-258.15)**1.5
        d_flux(ksn+1) = d_flux(ksn+1)*1000.0/snowdens
     endif
  endif
  !!  else  !! loosing water
  !!     call qtk(initp%snow_heat(ksn)/initp%sfcwater_mass(ksn),tempk,fracliq)
  !!     if(tempk .lt. 273.15) then !snow is evaporating, preserve density
  !!        dflux(ksn+1) = dflux(ksn+1)*snden
  !!     endif
  !  endif

  dinitp%avg_vapor_gc  = wflxgc   ! Diagnostic
  dinitp%avg_dew_cg    = dewgnd   ! Diagnostic
  dinitp%avg_qwshed_vg = qwshed   ! Diagnostic
  dinitp%avg_wshed_vg  = wshed    ! Diagnostic

  ! Transfer water downward through snow layers by percolation.
  ! Fracliq is the fraction of liquid in the snowcover or surface water.  wfree
  ! is the quantity of that liquid in kg/m2 which is free (not attached to
  ! snowcover) and therefore available to soak into the layer below).
  ! soilcap is the capacity of the top soil layer in kg/m2 to accept surface
  ! water.  wfree in the lowest snow layer is limited by this value.
  ! flxliq is the effective soaking flux in kg/m2/s over the timestep which
  ! is used to update the moisture in the top soil layer.
  if(ksn == 0)then
     if(w_flux(nzg+1) < 0.0)then
        ! gaining water, transfer everything to the virtual layer
        dinitp%virtual_heat = -qw_flux(nzg+1)
        dinitp%virtual_water = -w_flux(nzg+1)
        dinitp%virtual_depth = -d_flux(1)
        qw_flux(nzg+1) = 0.0
        w_flux(nzg+1) = wflxgc
     else
        w_flux(nzg+1) = w_flux(nzg+1) + wflxgc
     endif
  else
     dinitp%sfcwater_mass(ksn) = -w_flux(nzg+ksn+1) - wflxgc
     dinitp%sfcwater_energy(ksn) = dinitp%sfcwater_energy(ksn) -  &
          qw_flux(nzg+ksn+1)
     dinitp%sfcwater_depth(ksn) = -d_flux(ksn+1)
  endif

  dinitp%avg_smoist_gg(nzg) = w_flux(nzg+ksn+1)  ! Diagnostic

  ! Compute gravitational potential plus moisture potential,  
  ! psi + z (psiplusz) [m], liquid water content (soil_liq) [m], 
  ! and 99% the remaining water capacity (soilair99) [m].

  do k = lsl, nzg
     nsoil = csite%ntext_soil(k,ipa)
     psiplusz(k) = slzt(k) + soil(nsoil)%slpots  &
          * (soil(nsoil)%slmsts / real(initp%soil_water(k))) ** soil(nsoil)%slbs

     ! Soil liquid water must be converted to meters of liquid water per layer.
     ! This requires multiplication of volumetric water content, m3(water)/m3
     ! must be multiplied by depth to get a depth of water.

     soil_liq(k) = max(0.0, (sngl(initp%soil_water(k)) - soil(nsoil)%soilcp) *  &
                            initp%soil_fracliq(k))

     soilair99(k) = 0.99 * soil(nsoil)%slmsts - real(initp%soil_water(k))

  enddo


  ! Find amount of water transferred between soil layers (w_flux) [m]
  ! modulated by the liquid water fraction
  w_flux(nzg+1) = w_flux(nzg+1) * 1.0e-3 ! now in m/s

  !! Alternate surface infiltration  (mcd)
  !! based on surface conductivity not capacity
  !!  if(.false.) then
  if(infiltration_method .gt. 0) then
     print*,"running alt infiltation when we shouldn't be"
     if(initp%virtual_water .ne. 0.0) then  !!process "virtural water" pool
        nsoil = csite%ntext_soil(nzg,ipa)
        if(nsoil.ne.13) then

           call qtk(initp%virtual_heat/initp%virtual_water,tempk,fracliq)
           infilt = -dslzi(nzg)* 0.5 * slcons1(nzg,nsoil)  &
                * (real(initp%soil_water(nzg)) / soil(nsoil)%slmsts)**(2. * soil(nsoil)%slbs + 3.)  &
                * (psiplusz(nzg) - real(initp%virtual_water)/2000.0)  &  !!difference in potentials
                * .5 * (initp%soil_fracliq(nzg)+ fracliq)  !! mean liquid fraction

           !! adjust other rates accordingly
           w_flux(nzg+1) = w_flux(nzg+1) + infilt
           qw_flux(nzg+1)= qw_flux(nzg+1)+ infilt * (cliq1000 * tempk + alli1000)
           dinitp%virtual_water = dinitp%virtual_water - infilt*1000.
           dinitp%virtual_heat  = dinitp%virtual_heat  - infilt*(cliq1000 * tempk + alli1000)
        endif
     endif  !! end virtual water pool
     if(initp%nlev_sfcwater .ge. 1) then   !!process "snow" water pool 
        call qtk(initp%sfcwater_energy(1),tempk,fracliq)
        surface_water = initp%sfcwater_mass(1)*fracliq*0.001 !(m/m2)
        !        surface_heat  = surface_water* (cliq1000*tempk + alli1000)
        nsoil = csite%ntext_soil(nzg,ipa)
        if(nsoil.ne.13) then
           !! calculate infiltration rate (m/s?)
           infilt = -dslzi(nzg) * 0.5 * slcons1(nzg,nsoil)  &
                * (real(initp%soil_water(nzg)) / soil(nsoil)%slmsts)**(2. * soil(nsoil)%slbs + 3.)  &
                * (psiplusz(nzg) - surface_water/2.0)  &  !!difference in potentials
                * .5 * (initp%soil_fracliq(nzg)+ fracliq)  !! mean liquid fraction
           !! adjust other rates accordingly
           w_flux(nzg+1) = w_flux(nzg+1) + infilt
           qw_flux(nzg+1)= qw_flux(nzg+1)+ infilt * (cliq1000 * tempk + alli1000)
           dinitp%sfcwater_mass(1) = dinitp%sfcwater_mass(1) - infilt*1000.0
           dinitp%sfcwater_energy(1) = dinitp%sfcwater_energy(1) - infilt * (cliq1000 * tempk + alli1000) 
           dinitp%sfcwater_depth(1) = dinitp%sfcwater_depth(1) - infilt
        endif
     endif  !! end snow water pool

  endif  !! end alternate infiltration

  ! keep qw_flux in J/m2/s
  do k = lsl+1, nzg
     nsoil = csite%ntext_soil(k,ipa)
     if(nsoil /= 13 .and. csite%ntext_soil(k-1,ipa) /= 13)then

        wgpmid = 0.5 * real(initp%soil_water(k) + initp%soil_water(k-1))
        freezeCor = 0.5 * (initp%soil_fracliq(k)+ initp%soil_fracliq(k-1))
        if(freezeCor .lt. 1.0) freezeCor = 10.0**(-freezeCoef*(1.0-freezeCor))
        w_flux(k) = dslzti(k) * slcons1(k,nsoil)  &
             * (wgpmid / soil(nsoil)%slmsts)**(2. * soil(nsoil)%slbs + 3.)  &
             * (psiplusz(k-1) - psiplusz(k)) * freezeCor

        ! Limit water transfers to prevent over-saturation and over-depletion
        ! Compute q transfers between soil layers (qw_flux) [J/m2]
        ! Added some "lids" that exist on LEAF-3

	! SLOW POINT 1

        if (w_flux(k) > 0.) then
           if(soil_liq(k-1) <= 0.0 .or. soilair99(k) <= 0.0) then
              w_flux(k) = 0.0
           else
              w_flux(k) = min(w_flux(k),soil_liq(k-1),soilair99(k)*dslz(k)*.5)
           end if
        else
           if(soil_liq(k) <= 0.0 .or. soilair99(k-1) <= 0.0) then
              w_flux(k) = 0.0
           else
              w_flux(k) = -min(-w_flux(k),soil_liq(k),  &
                   soilair99(k-1)*dslz(k-1)*.5)
           end if
        end if
     endif
     

     qw_flux(k) = w_flux(k) * (cliq1000 * initp%soil_tempk(k) + alli1000)
     
     dinitp%avg_smoist_gg(k-1) = w_flux(k)*1000.0   ! Diagnostic

  enddo

  nsoil = csite%ntext_soil(lsl,ipa)
  if (nsoil /= 13 .and. isoilbc == 1) then
     !----- Free drainage -----------------------------------------------------------------!
     wgpmid      = real(initp%soil_water(lsl))
     freezeCor   = initp%soil_fracliq(lsl)
     if(freezeCor .lt. 1.0) freezeCor = 10.0**(-freezeCoef*(1.0-freezeCor))
     w_flux(lsl) = dslzti(lsl) * slcons1(lsl,nsoil)  &
                 * (wgpmid / soil(nsoil)%slmsts)**(2. * soil(nsoil)%slbs + 3.)  &
                 * freezeCor
     if (soil_liq(lsl) == 0.) w_flux(lsl) = 0.
     qw_flux(lsl) = w_flux(lsl) * (cliq1000 * initp%soil_tempk(lsl) + alli1000)
  else
     !----- Bedrock -----------------------------------------------------------------------!
     w_flux(lsl) = 0.
     qw_flux(lsl) = 0.
  end if

  ! Finally, update soil moisture (impose minimum value of soilcp) and q value.
  do k = lsl,nzg
     dinitp%soil_water(k) = dinitp%soil_water(k) - dble(dslzi(k) *   &
          ( w_flux(k+1) -  w_flux(k)))
     dinitp%soil_energy(k) =  dinitp%soil_energy(k)  - dslzi(k) *   &
          (qw_flux(k+1) - qw_flux(k))

  enddo

  ! Update soil moisture from transpiration/root uptake
  if (csite%lai(ipa) > lai_min) then
     do k1 = lsl, nzg    ! loop over extracted water
        do k2=k1,nzg
           if (csite%ntext_soil(k2,ipa) /= 13) then
              if (initp%available_liquid_water(k1) > 0.0) then
                 
                 wloss = 0.001 * initp%extracted_water(k1)                                 &
                       * soil_liq(k2) / initp%available_liquid_water(k1)
                 
                 dinitp%soil_water(k2) = dinitp%soil_water(k2) - dble(wloss)
                 
                 qwloss = wloss * (cliq1000 * initp%soil_tempk(k2) + alli1000)
                 
                 dinitp%soil_energy(k2) = dinitp%soil_energy(k2) - qwloss
                 
                 dinitp%avg_smoist_gc(k2)=dinitp%avg_smoist_gc(k2)-1000.0*wloss
                 
                 dinitp%ebudget_latent = dinitp%ebudget_latent + qwloss
              end if
           end if
        end do
     end do
  end if

  ! If we have a thin layer of snow the heat derivatives will require 
  ! special treatment.
  if(initp%nlev_sfcwater == 1 .and.   &
       initp%sfcwater_mass(1) < water_stab_thresh)then
     ! Top soil layer and surface water proceed in equilibrium
     dqwt = dinitp%soil_energy(nzg) * dslz(nzg)   &
         + dinitp%sfcwater_energy(1)
     dinitp%soil_energy(nzg) = dqwt * dslzi(nzg)
     dinitp%sfcwater_energy(1) = 0.0

  else
     ! Layer is computationally stable.

     ! Convert dinitp%sfcwater_energy to W/kg by the quotient rule:
     !    J/kg = J/m2 / kg/m2
     !  D(J/kg) = (  kg/m2 D(J/m2) - J/m2 D(kg/m2)   ) / (kg/m2)**2
     do k = 1, initp%nlev_sfcwater
        if(initp%sfcwater_mass(k) >= min_sfcwater_mass)then
           dinitp%sfcwater_energy(k) = (initp%sfcwater_mass(k) *  &
                dinitp%sfcwater_energy(k) - initp%sfcwater_energy(k) *   &
                dinitp%sfcwater_mass(k)) / initp%sfcwater_mass(k)**2 
        endif
     enddo
  endif

  return
end subroutine leaftw_derivs_ar

!==================================================================
subroutine canopy_derivs_two_ar(initp, dinitp, csite,ipa,isi,ipy, hflxgc, wflxgc,  &
     dewgndflx, wshed_tot, qwshed_tot, rhos, prss, pcpg, qpcpg, exner, &
     geoht, atm_tmp, lsl)
  
  use ed_state_vars,only: rk4patchtype,sitetype,patchtype
  use consts_coms, only : alvl, cp, cpi, day_sec, grav, alvi, umol_2_kgC, mmdry, mmdryi,t3ple
  use grid_coms, only : nzg
  use soil_coms, only : soil, dewmax
  use canopy_radiation_coms, only: lai_min
  use therm_lib, only : qwtk,rslif
  use misc_coms, only: dtlsm
  use ed_misc_coms, only: fast_diagnostics
  use allometry, only: dbh2ca
  implicit none


  type(sitetype),target :: csite
  type(patchtype),pointer :: cpatch
  type(rk4patchtype), target :: initp,dinitp 
  integer :: ipa,ico,ipy,isi

  integer, intent(in) :: lsl
  real, intent(in) :: rhos
  real, intent(in) :: atm_tmp
  real, intent(in) :: prss
  real, intent(in) :: pcpg
  real, intent(in) :: qpcpg
  real, intent(in) :: exner
  real, intent(in) :: geoht
  real, intent(out) :: hflxgc,wflxgc
  real, intent(out) :: wshed_tot,qwshed_tot
  real, intent(out) :: dewgndflx
  real, parameter :: leaf_h2o_thick = 0.11 !! mm
  real :: can_frac  !! total fractional canopy coverage
  logical,parameter :: debug = .true.

  real :: transp
  real :: cflxac
  real :: wflxac
  real :: hflxac

  real :: c2,c3

  real :: hflxvc
  real :: rasgnd
  real :: rbi
  real :: rd
  real :: sigmaw

  real :: wflxvc
  real :: wshed
  real :: qwshed
  real :: zoveg,zveg
  
  real :: wcapcan
  real :: wcapcani,hcapcani
  real :: cflxgc
  real :: laii
  real :: wflx
  
  real :: hflxvc_tot,transp_tot,cflxvc_tot,wflxvc_tot

  real :: rho_ustar
  real :: rdi
  real :: gpp_tot
  real :: storage_decay

  real :: leaf_flux
  real,external :: vertical_vel_flux

  real :: sat_shv
  real :: veg_temp,fracliq
  real :: max_leaf_water
  real :: maxfluxrate
  real :: qveg_water
  real :: intercepted,qintercepted
  real :: qwflxvc

  ! Canopies with LAI less than this number are assumed to be
  ! open, ie, some fraction of the rain-drops can reach
  ! the soil/litter layer unimpeded. 
  real,parameter :: lai_to_cover = 1.5

  ! Fluxes from atmosphere to canopy

  rho_ustar = rhos * initp%ustar

  ! Sensible Heat flux
  hflxac = rho_ustar * initp%tstar * exner

  ! Water flux 
  wflxac = rho_ustar * initp%rstar
  ! CO2 flux
  cflxac = rho_ustar * initp%cstar

  ! Surface roughness parameters.
  ! Eventually I should account for snow factors here.
  zoveg = csite%veg_rough(ipa)
  zveg = csite%veg_height(ipa)

  ! Capacities of the canopy air space
  wcapcan = rhos * max(zveg,3.5)
  wcapcani = 1.0 / wcapcan
  hcapcani = cpi * wcapcani

  ! The following value of ground-canopy resistance for the
  ! nonvegetated (bare soil or water) surface is from John Garratt.
  ! It is 5/ustar and replaces the one from old leaf.
  if(debug .and. abs(initp%ustar) < tiny(1.0)) print*,"USTAR = 0"
  rasgnd = 5. / initp%ustar


  cpatch => csite%patch(ipa)
  !! Calculate fraction of open canopy  
  can_frac = 1.0
  do ico = 1,cpatch%ncohorts
     if(cpatch%lai(ico) .gt. lai_min) then
        can_frac = can_frac*(1.0-min(1.0,cpatch%nplant(ico)*dbh2ca(cpatch%dbh(ico),cpatch%pft(ico))))
     endif
  enddo
  can_frac = 1.0 - can_frac

  if (csite%lai(ipa) > lai_min) then

     ! If vegetation is sufficiently abundant and not covered by snow, compute
     ! heat and moisture fluxes from vegetation to canopy, and flux resistance
     ! from soil or snow to canopy.
     c2 = max(0.,min(1., 0.509 * csite%lai(ipa)))
     rd = rasgnd * (1. - c2) + initp%rasveg * c2

     laii = 1.0/csite%lai(ipa)
     ! If the canopy does not cover all of the ground, then it should
     ! not intercept all of the water
     
     if(pcpg>0.0) then
        !! scale interception by canopy openess (MCD 01-12-09)
        intercepted  = pcpg * can_frac
        qintercepted = qpcpg * can_frac
        wshed_tot = pcpg-intercepted
        qwshed_tot = (wshed_tot/pcpg)*qpcpg
     else
        intercepted=0.0
        qintercepted = 0.0
        wshed_tot = 0.0
        qwshed_tot = 0.0
     endif

  else
     ! If the LAI is very small or if snow mostly covers the vegetation, bypass
     ! vegetation computations.  Set heat and moisture flux resistance rd 
     ! between the "canopy" and snow or soil surface to its bare soil value.  
     ! Set shed precipitation heat and moisture to unintercepted values.
     rd = rasgnd
     wshed_tot = pcpg
     qwshed_tot = qpcpg
     intercepted = 0.0
     qintercepted = 0.0

     ! Note: If the condition of low LAI for the entire patch was met, then
     ! it does not matter what the individual cohorts are normalized by,
     ! because they are effectively zero. So make sure the inverse patch LAI
     ! is a nominal non-zero/non-infinite number. This will only be used when
     ! parsing out intercepted leaf water into shed water; in which case
     ! the intercepted water is zero anyway. So this is just to prevent FPEs.
     
     laii = 1.0
     
  endif

  rdi = rhos / rd
  
  ! Compute sensible heat and moisture fluxes between top soil or snow surface
  ! and canopy air.  wflxgc [kg/m2/s] is the upward vapor flux from soil or 
  ! snow evaporation and dewgnd is the mass of dew that forms on the snow/soil
  ! surface this timestep; both are defined as always positive or zero.
  ! ---------------------------------------------------------------

  if(initp%nlev_sfcwater == 0)then
     hflxgc = cp * (initp%soil_tempk(nzg)   &
          - initp%can_temp) * rdi
  else

     hflxgc = cp * (initp%sfcwater_tempk(initp%nlev_sfcwater)   &
          - initp%can_temp) * rdi
  endif
  
  wflx = (initp%surface_ssh - initp%can_shv) * rdi
  
  ! Calculate the dew flux, and impose a dew cap
  ! ---------------------------------------------------------------

  dewgndflx = min(dewmax, max(0.0, -wflx))
  
  ! Final evap check, make sure that the projected integrated 
  ! evaporative mass flux does not exceed 55% of the available
  ! liquid water in the top soil layer. I don't like this 
  ! because it imposes layer thickness dependence, but this
  ! should make the model more stable.
  ! ----------------------------------------------------------------

  maxfluxrate = 0.55 * 1000. * initp%available_liquid_water(nzg)/dtlsm
  ! IF NO SURFACE WATER AND NOT DRY - EVAPORATE FROM SOIL PORES
  if(initp%nlev_sfcwater == 0 .and. initp%soil_water(nzg)  &
       > dble(soil(csite%ntext_soil(nzg,ipa))%soilcp)) then
     wflxgc = min(max(0.0, (initp%ground_shv - initp%can_shv) * rdi) &
          ,maxfluxrate)
  ! IF NO SURFACE WATER AND REALLY DRY - DONT EVAPORATE AT ALL
  else if ( initp%nlev_sfcwater == 0 .and. initp%soil_water(nzg)  &
       <= dble(soil(csite%ntext_soil(nzg,ipa))%soilcp)) then
     wflxgc = 0.0
  else
     wflxgc = max(0.,wflx)
  endif

  ! -----------------------------------------------------------------
  ! IF THERE IS ONLY 1mm OF WATER IN THE SOIL LAYER DONT
  ! EXTRACT.
  ! PROBLEM, THIS CALCULATION IS LAYER THICKNESS DEPENDANT.
  ! IF YOUR TOP LAYER IS ONLY 5 mm, THEN THIS WILL NULLIFY
  ! WATER FLUXES WHEN RELATIVE WATER CONTENT IS <20%, BUT
  ! IF YOU TOP LAYER IS 20mm THICK, THEN THIS THRESHOLD IS
  ! ABOUT THE SAME AS THE RESIDUAL WATER CONTENT.
  ! IN OTHER WORDS THIS CONDITION IS A KLUGE AND MAY NOT BE
  ! APPROPRIATE, UNLESS SOMEONE THINKS IT SHOULD STAY.
  ! RGK 11-2-08
  ! THE NEXT FEW LINES WILL BE REMOVED PENDING SENSITIVITY ANALYSIS

  if(initp%available_liquid_water(nzg) <= 0.01 .and. wflxgc > 0.0)then
     wflxgc = 0.0
  end if


  ! Loop over the cohorts in the patch. Calculate energy fluxes with
  ! surrounding canopy air space, integrate cohort energy, calculate
  ! precipitation throughfall and sum fluxes to the patch level.
  ! Initialize variables used to store sums over cohorts.
  ! ----------------------------------------------------------------

  hflxvc_tot = 0.0
  wflxvc_tot = 0.0
  cflxvc_tot = csite%cwd_rh(ipa)
  transp_tot = 0.0

  cflxgc = csite%rh(ipa) - csite%cwd_rh(ipa)
  gpp_tot = 0.0
  
  do ico = 1,cpatch%ncohorts
     
     cflxgc = cflxgc + cpatch%root_respiration(ico)
     
     ! calculate 'decay' term of storage (same for all)
     ! need to convert units from kgC/plant/day to umolC/m2/s.
     storage_decay = (cpatch%growth_respiration(ico) + cpatch%storage_respiration(ico)   &
          + cpatch%vleaf_respiration(ico)) & 
          * cpatch%nplant(ico) / (day_sec * umol_2_kgC)
     cflxvc_tot = cflxvc_tot + storage_decay
     
     
     ! See if this cohort has leaves, if not set the leaf energy derivatives to 
     ! zero, and pass all throughfall to the ground. Later, those small cohorts will
     ! have there leaf energy set to equilibrium with the canopy air space (temperature)
     ! ---------------------------------------------------------------------------------
     if(cpatch%lai(ico) > lai_min)then

        call qwtk (initp%veg_energy(ico),initp%veg_water(ico), &
        	   cpatch%hcapveg(ico),veg_temp,fracliq)

        !  Calculate leaf-level flux
        leaf_flux = cpatch%gpp(ico) - cpatch%leaf_respiration(ico)
        
        ! Update CO2 flux from vegetation to canopy air space.
        cflxvc_tot = cflxvc_tot - leaf_flux
        
        !  Calculate fraction of leaves covered with water
        if(initp%veg_water(ico) > 1.0e-12)then
           sigmaw = min(1.0,((initp%veg_water(ico) / (leaf_h2o_thick*cpatch%lai(ico)))**0.6667))
        else
           sigmaw = 0.0
        endif
                
        ! Here we need to be cautious though. The temperature may be off, but it can't be too off, like
        ! below 0K because this violates the basic laws of thermodynamics, and that does not make any sense.
        
        ! Use a temperature to calculate leaf surface saturation vapor pressure
        ! This conditioning prevents craziness with hot or cold leaves during 
        ! half-steps.  The vapor pressure sky-rockets when temperatures exceed
        ! the boiling point, and if this happens, it is possible to get extremely
        ! high evaporation rates, and then extremely high cooling rates, and then
        ! wild temperature fluctuations in small cohorts.
        ! So dont let the vapor pressure exceed that for nominal min and max temps
        ! But do not let the temperatures be capped below or above the canopy
        ! temperature, because if that happens, we will innapropriately chnage the sign of
        ! the fluxes.
        
        sat_shv=rslif(prss,veg_temp)

        c3 = cpatch%lai(ico) * rhos * (sat_shv - initp%can_shv)
        if(cpatch%rb(ico) < tiny(1.0)) then 
           cpatch%rb(ico) = 25.0*csite%lai(ipa)
           print*,"***WARNING*** uninitialized RB, setting to",cpatch%rb(ico)
        endif
        rbi = 1.0 / cpatch%rb(ico)
        
        if (c3 >= 0.) then  

           ! evapotranspiration
           ! Evaporation area factor being changed to 1.2 - RGK 11-2008
           wflxvc = c3 * sigmaw * 1.2 * rbi

           cpatch%Psi_open(ico)   = c3 / (cpatch%rb(ico) + cpatch%rsw_open(ico)  )
           cpatch%Psi_closed(ico) = c3 / (cpatch%rb(ico) + cpatch%rsw_closed(ico))
           if(initp%available_liquid_water(cpatch%krdepth(ico)) > 0.0 .and. veg_temp >= t3ple )then
              transp = cpatch%fsw(ico) * cpatch%Psi_open(ico) + (1.0 - cpatch%fsw(ico)) * cpatch%Psi_closed(ico)
           else
              transp = 0.0
           endif

        else   

           ! dew formation
           ! Dew area factor being changed to 1.2 - RGK 11-2008
           wflxvc = c3 * 1.2 * rbi
           transp = 0.0
           cpatch%Psi_open(ico) = 0.0
        endif

        dinitp%ebudget_latent = dinitp%ebudget_latent     + &
             wflxvc * (fracliq * alvl + (1.-fracliq) * alvi) + &
             transp * alvl

        ! We need to extract water from the soil equal to the transpiration
        initp%extracted_water(cpatch%krdepth(ico)) =   &
             initp%extracted_water(cpatch%krdepth(ico)) + transp

        ! Calculate vegetation-to-canopy fluxes 
        !   Energy (sensible heat)
        ! cpatchi%veg_temp = current vegetation temperature
        ! initp%can_temp = temperature of canopy air space
        ! rbi is the conductance
        ! rhos is the air density
        ! cp is the specific heat
        ! 2.2 accounts for stems and branches.
        ! --------------------------------------------------------

        hflxvc = 2.2 * cpatch%lai(ico) * cp * rhos * rbi  &
             * (veg_temp - initp%can_temp)
        
        ! If there is more leaf water (kg) than this threshold
        ! then no more water may be allowed to collect on the leaf
        max_leaf_water = leaf_h2o_thick*cpatch%lai(ico)

        !----------------------------------------------------------------
        !
        ! Calculate interception by leaves
        ! Added RGK 11-2008, comments welcomed
        !
        ! wflxvc accounts for evaporation and dew formation.  If
        ! the leaf has more water than the carrying capacity, then
        ! it must flux all precipitation and dew. The leaf may
        ! evaporation in every condition.
        ! ---------------------------------------------------------------

!! for a saturated leaf, why does the energy from dew transfer but not from rain?
!! seems like we should do both or not worry about either -- MCD
        
        ! Case 1: Leaf has no space for rain - evaporation dominates
        ! Assumptions: cooling from evaporation is removed from leaf and
        ! leaf water. Evaporation comes off of leaf water.  Rainfall
        ! and its internal energy immediately bypasses the leaf towards the ground.
        if(initp%veg_water(ico) >= max_leaf_water  .and. wflxvc >= 0. )then
           
           wshed = intercepted*cpatch%lai(ico)*laii
           qwshed = qintercepted*cpatch%lai(ico)*laii
           dinitp%veg_water(ico) = - wflxvc
           qveg_water = 0.

           ! Case 2: Leaf has no space for rain - dew dominates
           ! Assumptions: heating from the dew goes directly into the leaf
           ! and the leaf water, BUT, the mass of the dew goes into the
           ! shed water.  Rainfall and its internal energy bypass the leaf.
        else if(initp%veg_water(ico) >= max_leaf_water  .and. wflxvc < 0. )then
           
           wshed = intercepted*cpatch%lai(ico)*laii - wflxvc
           qwshed = qintercepted*cpatch%lai(ico)*laii ! The heat from dew is not shed
           dinitp%veg_water(ico) = 0.
           qveg_water = 0.

           ! Case 3: Leaf has space for rain - evaporation dominates
           ! Assumptions: Evaporation is removed from the leaf water, its cooling
           ! effects the leaf.  Rainfall and its internal energy accumulate
           ! on the leaf.
	   ! AND
	   ! Case 4: Leaf has space for rain - dew dominates
           ! Assumptions: Dew and its heating are applied to the leaf. Rainfall
           ! and its internal energy are applied to the leaf.
        else
           
           wshed = 0.0
           qwshed = 0.0
           dinitp%veg_water(ico) = -wflxvc + intercepted*cpatch%lai(ico)*laii
           qveg_water = qintercepted*cpatch%lai(ico)*laii
           
        endif

        !Latent heat of evap and dew
        qwflxvc = wflxvc * (fracliq * alvl + (1.-fracliq) * alvi)
        dinitp%veg_energy(ico) = &
             cpatch%rshort_v(ico)     &   ! Absorbed short wave radiation
             + cpatch%rlong_v(ico)    &   ! Net thermal radiation
             - hflxvc                 &   ! Sensible heat flux
             - qwflxvc                &   ! Evaporative phase cooling 
             - transp * alvl          &   ! Transpirative phase cooling
             + qveg_water                 ! Internal energy of intercepted water

!!$if(abs(dinitp%veg_energy(ico))>1000.0) then 
!!$   print*,"dVegE",dinitp%veg_energy(ico),&
!!$   cpatch%rshort_v(ico)&   ! Absorbed short wave radiation
!!$             , cpatch%rlong_v(ico)    &   ! Net thermal radiation
!!$             , -hflxvc                 &   ! Sensible heat flux
!!$             , -leflxvc                &   ! Evaporative phase cooling 
!!$             , -transp * alvl          &   ! Transpirative phase cooling
!!$             , heat_intercept_rate    !   ! Precipitation
!!$   print*,"hflxvc =",2.2 , cpatch%lai(ico) , cp , rhos , rbi  &
!!$        , veg_temp , initp%can_temp
!!$endif


        wflxvc_tot=wflxvc_tot+wflxvc
        hflxvc_tot=hflxvc_tot+hflxvc
        transp_tot=transp_tot+transp

        ! wshed:  Water passing through vegetated canopy to soil surface 
        ! (enters virtual layer first), [kg/m2/s]

        wshed_tot = wshed_tot + wshed
        qwshed_tot = qwshed_tot + qwshed


        if (.false.) then !ipy == 1 .and. ipa == 1 .and. ico == 1) then
           write (unit=62,fmt='(a)') '----------------------------------------------------------------------'
           write (unit=62,fmt='(a,1x,i12)')    ' - PFT:             ',cpatch%pft(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - NPLANT:          ',cpatch%nplant(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - BDEAD:           ',cpatch%bdead(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - BALIVE:          ',cpatch%balive(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - LAI:             ',cpatch%lai(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - DBH:             ',cpatch%dbh(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - HEIGHT:          ',cpatch%hite(ico)
           write (unit=62,fmt='(a)') ' '
           write (unit=62,fmt='(a,1x,es12.5)') ' - RHOS:            ',rhos
           write (unit=62,fmt='(a,1x,es12.5)') ' - PRSS:            ',prss
           write (unit=62,fmt='(a,1x,es12.5)') ' - ATM_TEMP:        ',atm_tmp
           write (unit=62,fmt='(a)') ' '
           write (unit=62,fmt='(a,1x,es12.5)') ' - CAN_TEMP:        ',initp%can_temp
           write (unit=62,fmt='(a,1x,es12.5)') ' - CAN_SHV:         ',initp%can_shv
           write (unit=62,fmt='(a)') ' '
           write (unit=62,fmt='(a,1x,es12.5)') ' - SIGMAW:          ',sigmaw
           write (unit=62,fmt='(a,1x,es12.5)') ' - VEG_WATER:       ',initp%veg_water(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - D(VEG_WATER)/DT: ',dinitp%veg_water(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - RBI(CONDUC):     ',rbi
           write (unit=62,fmt='(a,1x,es12.5)') ' - TRANSP:          ',transp
           write (unit=62,fmt='(a,1x,es12.5)') ' - WFLXVC:          ',wflxvc
           write (unit=62,fmt='(a,1x,es12.5)') ' - WSHED:           ',wshed
           write (unit=62,fmt='(a,1x,es12.5)') ' - INTERCEPTED:     ',intercepted*cpatch%lai(ico)*laii
           write (unit=62,fmt='(a)') ' '
           write (unit=62,fmt='(a,1x,es12.5)') ' - HCAPVEG:         ',cpatch%hcapveg(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - VEG_TEMP:        ',veg_temp
           write (unit=62,fmt='(a,1x,es12.5)') ' - VEG_ENERGY:      ',initp%veg_energy(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - D(VEG_ENERGY)/DT:',dinitp%veg_energy(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - RSHORT_V:        ',cpatch%rshort_v(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - RLONG_V:         ',cpatch%rlong_v(ico)
           write (unit=62,fmt='(a,1x,es12.5)') ' - HFLXVC:          ',hflxvc
           write (unit=62,fmt='(a,1x,es12.5)') ' - QWFLXVC:         ',qwflxvc
           write (unit=62,fmt='(a,1x,es12.5)') ' - QTRANSP:         ',transp*alvl
           write (unit=62,fmt='(a,1x,es12.5)') ' - QWSHED:          ',qwshed
           write (unit=62,fmt='(a,1x,es12.5)') ' - QINTERCEPT:      ',qintercepted*cpatch%lai(ico)*laii
           write (unit=62,fmt='(a)') '----------------------------------------------------------------------'
           write (unit=62,fmt='(a)') ' '
        end if

     else

        ! If there are no leaves, 
        dinitp%veg_energy(ico) = 0.0
        dinitp%veg_water(ico) = 0.0

        ! Allow the complete bypass of precipitation if there are no leaves
        ! Added RGK 11-2008, comments welcomed
        ! This will cause small deviations if the patch has LAI larger
        ! than lai_min, but only some of the cohorts do. If the ones that
        ! do not exceed the threshold of lai_min make up an appreciable
        ! fraction of the total lai, then we better account for it.
        wshed_tot = wshed_tot + intercepted*cpatch%lai(ico)*laii
        qwshed_tot = qwshed_tot + qintercepted*cpatch%lai(ico)*laii


     endif

  enddo  !  cohorts

  ! Update temperature and moisture of canopy.  hcapcan [J/m2/K] and
  ! wcapcan [kg_air/m2] are the heat and moisture capacities of   
  ! the canopy.
  ! --------------------------------------------------------------------------


  dinitp%can_temp = (hflxgc + hflxvc_tot + hflxac) * hcapcani

  
  dinitp%can_shv = (wflxgc - dewgndflx + wflxvc_tot + transp_tot +   &
       wflxac) * wcapcani


  ! Update co2 concentration in the canopy

  dinitp%can_co2 = ( (cflxgc + cflxvc_tot)*mmdry + cflxac) * wcapcani

  ! Integrate diagnostic variables - These are not activated
  ! unless fast file-type outputs are selected. This will speed up
  ! the integrator
  ! --------------------------------------------------------------------------

  if (fast_diagnostics) then

     dinitp%wbudget_loss2atm = - wflxac
     dinitp%ebudget_loss2atm = - hflxac
     dinitp%ebudget_latent = dinitp%ebudget_latent + (-dewgndflx + wflxgc) * alvi
     
     dinitp%co2budget_loss2atm = - cflxac * mmdryi
     dinitp%avg_gpp = gpp_tot
     dinitp%avg_carbon_ac = cflxac * mmdryi

     dinitp%avg_sensible_vc   = hflxvc_tot                             ! Sensible heat        vegetation -> canopy air
     dinitp%avg_sensible_2cas = hflxgc+hflxac+hflxvc_tot               ! Sensible heat        everywhere -> canopy air
     dinitp%avg_vapor_vc      = alvl*wflxvc_tot                        ! Latent heat          vegetation -> canopy air
     dinitp%avg_sensible_gc   = hflxgc                                 ! Sensible heat        ground     -> canopy air
     dinitp%avg_sensible_ac   = hflxac / exner                         ! Sensible heat        canopy air -> atmosphere
     dinitp%avg_vapor_ac      = alvl*wflxac                            ! Latent heat          canopy air -> atmosphere
     dinitp%avg_transp        = alvl*transp_tot                        ! Transpiration
     dinitp%avg_evap          = alvl*(wflxgc - dewgndflx + wflxvc_tot) ! Evaporation
     dinitp%avg_sensible_tot  = (hflxgc + hflxvc_tot)                  ! Sensible heat

     dinitp%avg_netrad = csite%rlong_g(ipa) + csite%rlong_s(ipa) + csite%rshort_g(ipa) + &
          sum(csite%rshort_s(1:initp%nlev_sfcwater,ipa),1)

     ! Auxillary variable
     !  dinitp%aux = dinitp%aux

  endif
  
  ! These variables below are virtual copies of the variables above, but are here for 
  ! for use in the coupled model. They form the set of canopy-atmospher fluxes that are
  ! used for turbulent closure. These variables are also zeroed and normalized
  ! every dtlsm timestep, the others are likely averaged over the analysis period
  
  dinitp%upwp = -(initp%ustar**2)
  dinitp%rpwp = -(initp%ustar*initp%rstar)
  dinitp%tpwp = -(initp%ustar*initp%tstar)
  if(debug .and. abs(atm_tmp) < tiny(1.0)) print*,"atm_tmp = 0"
  dinitp%wpwp = vertical_vel_flux(grav * geoht * cpi * exner / atm_tmp &
       ,initp%tstar,initp%ustar)

  return
end subroutine canopy_derivs_two_ar
!==================================================================
subroutine ed_stars(tha, rva, chia, thv, zpm, um, rough, ustar, rstar,   &
     tstar, cstar, can_shv, can_co2)

  use consts_coms, only: grav, vonk
  use canopy_air_coms, only: ustmin, ubmin
  implicit none

  real, intent(in) :: rough
  real, intent(out) :: ustar
  real, intent(out) :: rstar
  real, intent(out) :: tstar
  real, intent(out) :: cstar
  real, intent(in) :: can_shv
  real, intent(in) :: can_co2
  real :: b,csm,csh,d,a2,c1,ri,fm,fh,c2,cm,ch,c3,thv
  real :: rva,tha,zpm,um,chia,vels_pat
  ! Routine to compute Louis (1981) surface layer parameterization.

  b=5.
  csm=7.5
  csh=5.
  d=5.

  ! a2 is the drag coefficient in neutral conditions, here same for h/m
  ! ri is the bulk richardson numer, eq. 3.45 in Garratt

  vels_pat = max(um,ubmin)
  
  
  a2 = (vonk / log(zpm / rough)) ** 2
  c1 = a2 * vels_pat

  ri = grav * zpm * (tha - thv)  / (.5 * (tha + thv) * vels_pat * vels_pat)
  if (tha - thv .gt. 0.) then
     ! STABLE CASE
     fm = 1. / (1. + (2. * b * ri / sqrt(1. + d * ri)))
     fh = 1. / (1. + (3. * b * ri * sqrt(1. + d * ri)))
  else
     ! UNSTABLE CASE
     c2 = b * a2 * sqrt(zpm / rough * (abs(ri)))
     cm = csm * c2
     ch = csh * c2
     fm = (1. - 2. * b * ri / (1. + 2. * cm))
     fh = (1. - 3. * b * ri / (1. + 3. * ch))
  endif

  ustar = max(ustmin,sqrt(c1 * vels_pat * fm))
  c3 = c1 * fh / ustar
  rstar = c3 * (rva - can_shv)
  tstar = c3 * (tha - thv)
  cstar = c3 * (chia - can_co2)

  return
end subroutine ed_stars

!====================================================================
real function vertical_vel_flux(gzotheta,tstar,ustar)
  
  
!  include 'rconstants.h'

!  real :: vt2dd,vt2de,tstar,ustar,psin
!  real, parameter :: cc=4.7,wtol=1.0e-20
!  real :: sqrarg,numarg
  use consts_coms, only: vonk
  
  implicit none
  ! Arguments:
  real, intent(in)    :: ustar, tstar, gzotheta

  ! Local Variables:
  real :: zoverl
  real :: wtol,cx,psin
  
  data wtol/1e-20/
  
  !  if(ustar.gt.0.0)then
  !     vt2de = vt2dd * vonk * tstar / (ustar * ustar)
  !  else
  !     vt2de = 0.0
  !  endif
  !  if (vt2de .lt. 0.0)then
  !     vt2de = vt2de * sqrt(sqrt(1.0 - 15.0 * vt2de))
  !  else
  !     vt2de = vt2de / (1.0 + cc * vt2de)
  !  endif
  !  numarg = (1.0 + vt2de * (-5.39 + vt2de * 6.998 ))
  !  if(numarg .ne. 0.0)then
  !     sqrarg = (1.0 - 2.86 * vt2de) / numarg
  !  else
  !     sqrarg = 0.0
  !  endif
  !  if(sqrarg .gt.0.0)then
  !     psin = sqrt(sqrarg)
  !  else
  !     psin = 0.0
  !  endif

  zoverl = gzotheta * vonk * tstar / (ustar * ustar)
  
  if (zoverl < 0.)then
     cx = zoverl * sqrt(sqrt(1. - 15. * zoverl))
  else
     cx = zoverl / (1.0 + 4.7 * zoverl)
  endif
  
  psin = sqrt((1.-2.86 * cx) / (1. + cx * (-5.39 + cx * 6.998 )))
  vertical_vel_flux =  (0.27 * max(6.25 * (1. - cx) * psin,wtol)  &
       - 1.18 * cx * psin) * ustar * ustar
  
  return
end function vertical_vel_flux
