subroutine leaf_derivs_ar(initp, dinitp, csite, ipa,isi,ipy, rhos, prss, pcpg, qpcpg,   &
     atm_tmp, exner, geoht, vels, atm_shv, atm_co2, lsl)
  
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

  integer :: k

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
  use consts_coms, only : alvl, cliq1000, cpi, alvi, alli1000, t3ple,tsupercool
  use grid_coms, only: nzg,nzs
  use soil_coms, only : soil, slz, dslz, dslzi, water_stab_thresh,  &
       infiltration_method, dslzti, slcons1, slzt, min_sfcwater_mass,ss
  use misc_coms, only: current_time
  use canopy_radiation_coms, only: lai_min

  use ed_state_vars,only:sitetype,patchtype,rk4patchtype
  
  use therm_lib, only : qtk, qwtk, qwtk8

  implicit none

  type(sitetype),target :: csite
  type(patchtype),pointer :: cpatch
  integer :: ipa,ico,isi,ipy
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
  real, dimension(nzgmax,nzgmax) :: thick
  integer :: ksnnew,kk
  real :: qwshed,qwfree,wfree,depthgain,schar4c,dqw,dw,w,qw,wfreeb,depthloss
  
  !----------
  ! This is the free surface water transfer inverse time.  I made this up.  It 
  ! requires testing.
  !  real, parameter :: taui = 1.0/30.0 !!!!  Standard
  real, parameter :: taui = 1.0/300.0
  !----------
  real :: sndenmin,sndenmax,soilcap,totsnow,wgpmid,transp,wloss,fracw,zibar
  real :: dqwt,qwliq0,soilfrac,qwt

  real :: tgpsum,wgpsum,slcpdsum,tempktopm,fracliq,wateradd,tempk
  integer :: ksat,k1,k2
  real :: qwloss
  real :: surface_water !! pool of availible liquid water on the soil surface (kg/m2?)
  real :: infilt !! surface infiltration rate
  real :: snowdens !! snow density (kg/m2)
  real :: freezeCor !! correction to conductivity for partially frozen soil
  real :: sum_hflux

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
  w_flux = 0.0
  qw_flux = 0.0
  d_flux = 0.0

  ! Initialize derivatives to zero
  ! initp contains the current values of the state variables of the patch
  ! dinitp contains the derivative of the state variable in the patch
  do k=lsl,nzg
     dinitp%soil_energy(k) = 0.0
     dinitp%soil_water(k) = 0.0
  enddo

  ! Initialize derivative of snow layers
  do k=1,nzs
     dinitp%sfcwater_depth(k) = 0.0
     dinitp%sfcwater_energy(k) = 0.0  ! THIS IS IN W/m2!!!  Convert at end.
     dinitp%sfcwater_mass(k) = 0.0
  enddo
  dinitp%virtual_heat = 0.0
  dinitp%virtual_water = 0.0
  dinitp%virtual_depth = 0.0

  ! Calculate diagnostic mixing ratios
  ksn = initp%nlev_sfcwater
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
  ! SLZ is specified in RAMSIN.  Each element of the array sets the value of the bottom of a corresponding soil layer.
  ! Eg, SLZ = -2, -1, -0.5, -0.25.  There are four soil layers in this example; soil layer 1 goes from 2 meters below the 
  ! surface to 1 meter below the surface.
  nsoil = csite%ntext_soil(nzg,ipa)
  initp%available_liquid_water(nzg) = dslz(nzg) * max(0.0,  &
       initp%soil_fracliq(nzg) * (real(initp%soil_water(nzg)) - soil(nsoil)%soilcp))

  ! initialized to zero
  initp%extracted_water(nzg) = 0.0
  do k = nzg - 1, lsl, -1
     nsoil = csite%ntext_soil(k,ipa)
     initp%available_liquid_water(k) = initp%available_liquid_water(k+1) +  &
          dslz(k) * max(0.0, (real(initp%soil_water(k)) - soil(nsoil)%soilcp) *  &
          initp%soil_fracliq(k))
     initp%extracted_water(k) = 0.0
  enddo

  ! Get derivatives of canopy variables
  call canopy_derivs_two_ar(initp, dinitp, csite, ipa,isi,ipy,hflxgc, wflxgc, dewgnd,  &
       wshed,qwshed, rhos, prss, pcpg, qpcpg, exner, geoht, atm_tmp, lsl)

  
  ! Calculate conductivities:
  !        soil
  do k = lsl, nzg
     nsoil = csite%ntext_soil(k,ipa)
     if(nsoil <= 12)then
        wgpfrac = min(real(initp%soil_water(k)) / soil(nsoil)%slmsts,1.0)
        soilcond = soil(nsoil)%soilcond0 + wgpfrac * (soil(nsoil)%soilcond1  &
             + wgpfrac * soil(nsoil)%soilcond2)
     else
        soilcond=soil(nsoil)%soilcond0
     endif
     rfactor(k) = dslz(k) / soilcond
  enddo

  !        water
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
          soil(csite%ntext_soil(nzg,ipa))%slcpd,tempk,fracliq)

  endif

  !      heat flux (hfluxgsc) at soil or sfcwater top from longwave, sensible
  !      [W/m^2]
  hfluxgsc(nzg+ksn+1) = hflxgc + wflxgc * (fracliq*alvl+(1-fracliq)*alvi) - csite%rlong_g(ipa) - csite%rlong_s(ipa)
  
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

  qw_flux(nzg+ksn+1) = - dewgnd * alvi - qwshed
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
        if(tempk > 258.15)snowdens=50+1.5*(tempk-258.15)**1.5
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
          * (soil(nsoil)%slmsts / initp%soil_water(k)) ** soil(nsoil)%slbs

     ! Soil liquid water must be converted to meters of liquid water per layer.
     ! This requires multiplication of volumetric water content, m3(water)/m3
     ! must be multiplied by depth to get a depth of water.

     soil_liq(k) = max(0.0, (real(initp%soil_water(k)) - soil(nsoil)%soilcp) *  &
                            initp%soil_fracliq(k))

     soilair99(k) = 0.99 * soil(nsoil)%slmsts - initp%soil_water(k)

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
                * (initp%soil_water(nzg) / soil(nsoil)%slmsts)**(2. * soil(nsoil)%slbs + 3.)  &
                * (psiplusz(nzg) - initp%virtual_water/2000)  &  !!difference in potentials
                * .5 * (initp%soil_fracliq(nzg)+ fracliq)  !! mean liquid fraction

           !! adjust other rates accordingly
           w_flux(nzg+1) = w_flux(nzg+1) + infilt
           qw_flux(nzg+1)= qw_flux(nzg+1)+ infilt * cliq1000 * (tempk - tsupercool)
           dinitp%virtual_water = dinitp%virtual_water - infilt*1000
           dinitp%virtual_heat  = dinitp%virtual_heat  - infilt*cliq1000 * (tempk - tsupercool)
        endif
     endif  !! end virtual water pool
     if(initp%nlev_sfcwater .ge. 1) then   !!process "snow" water pool 
        call qtk(initp%sfcwater_energy(1),tempk,fracliq)
        surface_water = initp%sfcwater_mass(1)*fracliq*0.001 !(m/m2)
        !        surface_heat  = surface_water*(tempk - tsupercool)*cliq1000
        nsoil = csite%ntext_soil(nzg,ipa)
        if(nsoil.ne.13) then
           !! calculate infiltration rate (m/s?)
           infilt = -dslzi(nzg) * 0.5 * slcons1(nzg,nsoil)  &
                * (initp%soil_water(nzg) / soil(nsoil)%slmsts)**(2. * soil(nsoil)%slbs + 3.)  &
                * (psiplusz(nzg) - surface_water/2)  &  !!difference in potentials
                * .5 * (initp%soil_fracliq(nzg)+ fracliq)  !! mean liquid fraction
           !! adjust other rates accordingly
           w_flux(nzg+1) = w_flux(nzg+1) + infilt
           qw_flux(nzg+1)= qw_flux(nzg+1)+ infilt * cliq1000 * (tempk - tsupercool)
           dinitp%sfcwater_mass(1) = dinitp%sfcwater_mass(1) - infilt*1000
           dinitp%sfcwater_energy(1) = dinitp%sfcwater_energy(1) - infilt * cliq1000 * (tempk - tsupercool) 
           dinitp%sfcwater_depth(1) = dinitp%sfcwater_depth(1) - infilt
        endif
     endif  !! end snow water pool

  endif  !! end alternate infiltration

  ! keep qw_flux in J/m2/s
  do k = lsl+1, nzg
     nsoil = csite%ntext_soil(k,ipa)
     if(nsoil /= 13 .and. csite%ntext_soil(k-1,ipa) /= 13)then

        wgpmid = 0.5 * (initp%soil_water(k) + initp%soil_water(k-1))
        freezeCor = 0.5 * (initp%soil_fracliq(k)+ initp%soil_fracliq(k-1))
        if(freezeCor .lt. 1.0)freezeCor = 10**(-7*(1-freezeCor))
        w_flux(k) = dslzti(k) * slcons1(k,nsoil)  &
             * (wgpmid / soil(nsoil)%slmsts)**(2. * soil(nsoil)%slbs + 3.)  &
             * (psiplusz(k-1) - psiplusz(k)) * freezeCor

        

        ! Limit water transfers to prevent over-saturation and over-depletion
        ! Compute q transfers between soil layers (qw_flux) [J/m2]
        ! Added some "lids" that exist on LEAF-3
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

     qw_flux(k) = w_flux(k) * cliq1000 * (initp%soil_tempk(k) - tsupercool)
     
     dinitp%avg_smoist_gg(k-1) = w_flux(k)*1000   ! Diagnostic

  enddo

  ! Finally, update soil moisture (impose minimum value of soilcp) and q value.
  do k = lsl,nzg
     dinitp%soil_water(k) = dinitp%soil_water(k) - dslzi(k) *   &
          ( w_flux(k+1) -  w_flux(k))
     dinitp%soil_energy(k) =  dinitp%soil_energy(k)  - dslzi(k) *   &
          (qw_flux(k+1) - qw_flux(k))

     ! MLO- Use the loop to zero transpiration flux...
     dinitp%avg_smoist_gc(k) = 0.0
     
     ! Auxillary Diangnostic Variable for soil
     ! Set as soil width for a check
      dinitp%aux_s(k) = 0.0 ! dslz(k)
  enddo

  ! Update soil moisture from transpiration
  !  if(csite%lai(ipa) > lai_min)then
  do k1 = lsl, nzg    ! loop over extracted water
     
     do k2=k1,nzg
        if(csite%ntext_soil(k2,ipa) /= 13) then
           if(initp%available_liquid_water(k1) > 0.0)then
              
              wloss = 0.001 * initp%extracted_water(k1)   &
                   * soil_liq(k2) / initp%available_liquid_water(k1)
              
              dinitp%soil_water(k2) = dinitp%soil_water(k2) - wloss
              
              qwloss = wloss * (cliq1000 * (initp%soil_tempk(k2) - t3ple) + &
                   alli1000)
              dinitp%soil_energy(k2) = dinitp%soil_energy(k2) - qwloss
              
              dinitp%avg_smoist_gc(k2)=dinitp%avg_smoist_gc(k2)-1000*wloss
              
              dinitp%ebudget_latent = dinitp%ebudget_latent + qwloss
              
           elseif(initp%extracted_water(k1) > 0.0)then
              print*,initp%extracted_water(k1),  &
                   initp%available_liquid_water(k1)
              print*,k1,initp%available_liquid_water(lsl:nzg)
              print*,current_time%time,initp%soil_fracliq(lsl:nzg)
              print*,'Model is trying to extract water '
              print*,'(via transpiration) from soil '
              print*,'that is either frozen or otherwise has no water.'
              stop
           endif
        endif
     enddo
  enddo
     !  endif

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
  
  use ed_state_vars,only: rk4patchtype,sitetype,patchtype,edgrid_g
 
  use consts_coms, only : alvl, cp, cpi, day_sec, grav, alvi,   &
       alli, cliq, cice, t3ple, umol_2_kgC, mmdry, mmdryi
  use grid_coms, only : nzg
  use soil_coms, only : soil, dewmax
  use canopy_radiation_coms, only: lai_min

  use canopy_air_coms, only: hcapveg_ref,heathite_min

  use pft_coms, only: q, qsw, water_conductance,leaf_width,rho

  use therm_lib, only : rslif,qwtk
  use misc_coms, only: dtlsm
  use ed_misc_coms, only: fast_diagnostics

  implicit none

  type(sitetype),target :: csite
  type(patchtype),pointer :: cpatch
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
  real :: transp,qpcpg_lw,dumarg,rsat,cflxac
  real :: wflxac,hflxac
  type(rk4patchtype), target :: initp,dinitp 
  integer :: k,iter_can
  real :: fac,aux,c2,c3,factv,fracliqv,hcapveg,hflxvc
  real :: qwtot,rasgnd,rasveg,rb,rbi,rd,rsatveg,sigmaw,tvegk
  real :: wflxvc,wshed,wtemp,zoveg,zdisp,zveg
  real, intent(out) :: dewgndflx
  real :: hcapcan,wcapcan,cflxgc,LAII,wflx  
  real :: hflxvc_tot,transp_tot,cflxvc_tot,wflxvc_tot,root_res_fac
  real :: A_op,P_op,A_cl,P_cl,gr_resp,A_net,P_net,cflxvc
  real :: leaf_flux_max,plant_flux_max,tvegc
  integer :: iipft
  real :: leaf_flux,plant_flux,qpcpg_k
  real :: leaf_flux_pot,plant_flux_pot,a_net_max
  real :: nitrogen_supply,nstepi
  real :: lhpc
  integer :: first_cohort
  real :: potential_water,qwshed,laicum,deltamr,rho_ustar,rdi,gpp_tot
  real :: storage_decay,vertical_vel_flux
  real :: heat_intercept_rate,dew_evap_latent_heat_loss
  real :: w_demand,w_supply,wcapcani,hcapcani,broot
  real :: sat_shv,veg_temp,fracliq
  real :: max_leaf_water
  real :: maxfluxrate
  real :: qveg_water

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
  wcapcan = rhos * zveg
  wcapcani = 1.0 / wcapcan
  hcapcani = cpi * wcapcani

  ! The following value of ground-canopy resistance for the
  ! nonvegetated (bare soil or water) surface is from John Garratt.
  ! It is 5/ustar and replaces the one from old leaf.
  rasgnd = 5. / initp%ustar

  laii = 1.0/csite%lai(ipa)

  if (csite%lai(ipa) > lai_min) then

     ! If vegetation is sufficiently abundant and not covered by snow, compute
     ! heat and moisture fluxes from vegetation to canopy, and flux resistance
     ! from soil or snow to canopy.
     c2 = max(0.,min(1., 0.509 * csite%lai(ipa)))
     rd = rasgnd * (1. - c2) + initp%rasveg * c2
     qwshed_tot = 0.0
     wshed_tot = 0.0
  else
     ! If the LAI is very small or if snow mostly covers the vegetation, bypass
     ! vegetation computations.  Set heat and moisture flux resistance rd 
     ! between the "canopy" and snow or soil surface to its bare soil value.  
     ! Set shed precipitation heat and moisture to unintercepted values.
     rd = rasgnd
     wshed_tot = pcpg
     qwshed_tot = qpcpg
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
       > soil(csite%ntext_soil(nzg,ipa))%soilcp) then
     wflxgc = min(max(0.0, (initp%ground_shv - initp%can_shv) * rdi) &
          ,maxfluxrate)
  ! IF NO SURFACE WATER AND REALLY DRY - DONT EVAPORATE AT ALL
  else if ( initp%nlev_sfcwater == 0 .and. initp%soil_water(nzg)  &
       <= soil(csite%ntext_soil(nzg,ipa))%soilcp) then
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
  
  cpatch => csite%patch(ipa)
  
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


        ! Effective heat capacity of vegetation [J K-1] = [m3] * [J m-3 K-1]
        hcapveg = hcapveg_ref * max(cpatch%hite(1),heathite_min) * cpatch%lai(ico) * laii

        call qwtk (initp%veg_energy(ico),initp%veg_water(ico),hcapveg,veg_temp,fracliq)

        !  Calculate leaf-level flux
        leaf_flux = cpatch%gpp(ico) - cpatch%leaf_respiration(ico)

        ! Update CO2 flux from vegetation to canopy air space.
        cflxvc_tot = cflxvc_tot - leaf_flux

        !  Calculate fraction of leaves covered with water
        if(initp%veg_water(ico) > 1.0e-12)then
           sigmaw = min(1.0,((initp%veg_water(ico) / (0.22*cpatch%lai(ico)))**0.6667))
        else
           sigmaw = 0.0
        endif
        
        ! Removed the check on vegetation temperature - RGK 11-2008
        ! Reason: State variables will be strange during partial steps, best to look
        ! at veg temperatures after integrations are complete.
        ! if (veg_temp < 183.15) then


        sat_shv=rslif(prss,veg_temp)
        c3 = cpatch%lai(ico) * rhos * (sat_shv - initp%can_shv)
        
        
        rbi = 1.0 / cpatch%rb(ico)

        if (c3 >= 0.) then  
           ! evapotranspiration
           ! Evaporation area factor being changed to 1.2 - RGK 11-2008
           wflxvc = c3 * sigmaw * 1.2 * rbi
           cpatch%Psi_open(ico)   = c3 / (cpatch%rb(ico) + cpatch%rsw_open(ico)  )
           cpatch%Psi_closed(ico) = c3 / (cpatch%rb(ico) + cpatch%rsw_closed(ico))
           if(initp%available_liquid_water(cpatch%krdepth(ico)) > 0.0)then
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
        ! 2.2 acpatchounts for stems and branches.
        ! --------------------------------------------------------

        hflxvc = 2.2 * cpatch%lai(ico) * cp * rhos * rbi  &
             * (veg_temp - initp%can_temp)

        
        ! If there is more leaf water (kg) than this threshold
        ! then no more water may be allowed to collect on the leaf
        max_leaf_water = 0.22*cpatch%lai(ico)


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


        ! Case 1: Leaf has no space for rain - evaporation dominates
        ! Assumptions: cooling from evaporation is removed from leaf and
        ! leaf water. Evaporation comes off of leaf water.  Rainfall
        ! and its internal energy immediately bypasses the leaf towards the ground.
        if(initp%veg_water(ico) >= max_leaf_water  .and. wflxvc >= 0. )then
           
           wshed = pcpg*cpatch%lai(ico)*laii
           qwshed = qpcpg*cpatch%lai(ico)*laii
           dinitp%veg_water(ico) = - wflxvc
           qveg_water = 0.

           ! Case 2: Leaf has no space for rain - dew dominates
           ! Assumptions: heating from the dew goes directly into the leaf
           ! and the leaf water, BUT, the mass of the dew goes into the
           ! shed water.  Rainfall and its internal energy bypass the leaf.
        else if(initp%veg_water(ico) >= max_leaf_water  .and. wflxvc < 0. )then
           
           wshed = pcpg*cpatch%lai(ico)*laii - wflxvc
           qwshed = qpcpg*cpatch%lai(ico)*laii ! The heat from dew is not shed
           dinitp%veg_water(ico) = 0.
           qveg_water = 0.

           ! Case 3: Leaf has space for rain - evaporation dominates
           ! Assumptions: Evaporation is removed from the leaf water, its cooling
           ! effects the leaf.  Rainfall and its internal energy accumulate
           ! on the leaf.
        else if(initp%veg_water(ico) < max_leaf_water  .and. wflxvc >= 0. )then
           
           wshed = 0.0
           qwshed = 0.0
           dinitp%veg_water(ico) = -wflxvc + pcpg*cpatch%lai(ico)*laii
           qveg_water = qpcpg*cpatch%lai(ico)*laii

           ! Case 4: Leaf has space for rain - dew dominates
           ! Assumptions: Dew and its heating are applied to the leaf. Rainfall
           ! and its internal energy are applied to the leaf.
        else if (initp%veg_water(ico) < max_leaf_water  .and. wflxvc < 0. )then
           
           wshed = 0.0
           qwshed = 0.0
           dinitp%veg_water(ico) = -wflxvc + pcpg*cpatch%lai(ico)*laii
           qveg_water = qpcpg*cpatch%lai(ico)*laii
           
        endif

        dinitp%veg_energy(ico) = &
             cpatch%rshort_v(ico)     &   ! Absorbed short wave radiation
             + cpatch%rlong_v(ico)    &   ! Net thermal radiation
             - hflxvc                 &   ! Sensible heat flux
             - wflxvc * (fracliq * alvl + (1.-fracliq) * alvi) &
                                          ! Latent heat of evap and dew
             - transp * alvl          &   ! Transpirative phase cooling (assume liquid stomata?)
             + qveg_water                 ! Internal energy of intercepted water

        wflxvc_tot=wflxvc_tot+wflxvc
        hflxvc_tot=hflxvc_tot+hflxvc
        transp_tot=transp_tot+transp

        ! wshed:  Water passing through vegetated canopy to soil surface 
        ! (enters virtual layer first), [kg/m2/s]

        wshed_tot = wshed_tot + wshed
        qwshed_tot = qwshed_tot + qwshed

        if(initp%veg_energy(ico) .ne.initp%veg_energy(ico) )then
           call fatal_error('initp%veg_energy is NaN','canopy_derivs_two_ar','rk4_derivs.F90')
        end if

        if(dinitp%veg_energy(ico) .ne.dinitp%veg_energy(ico) ) then
           print*, 'dinitp%veg_energy is nan'
           print*, dinitp%veg_energy(ico), initp%veg_energy(ico), veg_temp, &
                   dinitp%veg_water(ico), cpatch%rshort_v(ico),             &
                   cpatch%rlong_v(ico), hflxvc, wflxvc, transp, alvl,       &
                   heat_intercept_rate
           call fatal_error('dinitp%veg_energy is NaN','canopy_derivs_two_ar','rk4_derivs.F90')
        endif

     else

        ! If there are no leaves, 
        dinitp%veg_energy(ico) = 0.0
        dinitp%veg_water(ico) = 0.0

        ! Allow the complete bypass of precipitation if there are no leaves
        ! Added RGK 11-2008, comments welcomed
        wshed_tot = wshed_tot + pcpg*cpatch%lai(ico)*laii
        qwshed_tot = qwshed_tot + qpcpg*cpatch%lai(ico)*laii


     endif

  enddo  !  cohorts

  ! Update temperature and moisture of canopy.  hcapcan [J/m2/K] and
  ! wcapcan [kg_air/m2] are the heat and moisture capacities of   
  ! the canopy.
  ! Modified by RGK 11-2008, comments welcomed (reverting it back to original,
  ! by removing the vapor heating cooling)
  ! --------------------------------------------------------------------------
  
  !  dinitp%can_temp = (hflxgc + hflxvc_tot + leflxvc_tot + hflxac) * hcapcani

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
     
     ! Auxillary variable
     
     !  dinitp%aux = dinitp%aux

  endif

  
  ! These variables below are virtual copies of the variables above, but are here for 
  ! consistency's sake. They form the set of canopy-atmospher fluxes that are
  ! used for turbulent closure. These variables are also zeroed and normalized
  ! every dtlsm timestep, the others are likely averaged over the analysis period
  
  dinitp%upwp = -(initp%ustar**2)
  dinitp%rpwp = -(initp%ustar*initp%rstar)
  dinitp%tpwp = -(initp%ustar*initp%tstar)
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
  real :: rva,tha,zpm,um,chistar,chia,vels_pat
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
     fm = 1. / (1. + (2 * b * ri / sqrt(1 + d * ri)))
     fh = 1. / (1. + (3 * b * ri * sqrt(1 + d * ri)))
  else
     ! UNSTABLE CASE
     c2 = b * a2 * sqrt(zpm / rough * (abs(ri)))
     cm = csm * c2
     ch = csh * c2
     fm = (1. - 2 * b * ri / (1. + 2 * cm))
     fh = (1. - 3 * b * ri / (1. + 3 * ch))
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
  real :: wtol,cosine1,sine1,vtscr,cx,psin
  
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
