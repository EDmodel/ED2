module rk4_driver_ar
  
contains

! Main driver of short-time scale dynamics of the land surface model.
!-------------------------------------------------------
  subroutine rk4_timestep_ar(cgrid,ifm,integration_buff)

    use ed_state_vars,only:integration_vars_ar,edtype,polygontype,sitetype,patchtype
    use grid_coms, only: nzg
    use max_dims, only : n_dbh
    use misc_coms, only: dtlsm
    use consts_coms, only: umol_2_kgC

    implicit none

    type(integration_vars_ar), target :: integration_buff

    type(edtype),target       :: cgrid
    type(polygontype),pointer :: cpoly
    type(sitetype),pointer    :: csite
    type(patchtype),pointer   :: cpatch
    integer :: ifm,ipy,isi,ipa
        
    integer, dimension(nzg) :: ed_ktrans
    real :: sum_lai_rbi
    real, external :: compute_netrad_ar
    real :: gpp
    real, dimension(n_dbh) :: gpp_dbh
    real :: plant_respiration

   ! print*,cgrid%polygon(1)%site(1)%patch(1)%lai(1)
    
    polygonloop: do ipy = 1,cgrid%npolygons
       
       cpoly => cgrid%polygon(ipy)
       
       siteloop: do isi = 1,cpoly%nsites
          
          csite => cpoly%site(isi)
          
          patchloop: do ipa = 1,csite%npatches
             
             cpatch => csite%patch(ipa)


             ! Get velocity for aerodynamic resistance
             if(csite%can_temp(ipa) < cpoly%met(isi)%atm_tmp)then
                cpoly%met(isi)%vels = cpoly%met(isi)%vels_stab
             else
                cpoly%met(isi)%vels = cpoly%met(isi)%vels_unstab
             endif

             ! Get photosynthesis, stomatal conductance, and transpiration
             call canopy_photosynthesis_ar(csite, ipa, cpoly%met(isi)%vels,   &
                  cpoly%met(isi)%rhos,   &
                  cpoly%met(isi)%prss,  &
                  ed_ktrans, csite%ntext_soil(:,ipa), csite%soil_water(:,ipa),   &
                  csite%soil_fracliq(:,ipa), cpoly%lsl(isi), sum_lai_rbi,  &
                  cpoly%leaf_aging_factor(:,isi),  &
                  cpoly%green_leaf_factor(:,isi))

             ! Compute root and heterotrophic respiration
             call soil_respiration_ar(csite,ipa)

             csite%wbudget_precipgain(ipa) = csite%wbudget_precipgain(ipa) +   &
                  cpoly%met(isi)%pcpg * dtlsm
             csite%ebudget_precipgain(ipa) = csite%ebudget_precipgain(ipa) +   &
                  cpoly%met(isi)%qpcpg * dtlsm
             csite%ebudget_netrad(ipa) = csite%ebudget_netrad(ipa) +  &
                  compute_netrad_ar(csite,ipa) * dtlsm

             call sum_plant_cfluxes_ar(csite,ipa,gpp,gpp_dbh,plant_respiration)

             csite%co2budget_gpp(ipa) = csite%co2budget_gpp(ipa) + gpp * dtlsm
             csite%co2budget_gpp_dbh(:,ipa) = csite%co2budget_gpp_dbh(:,ipa) + gpp_dbh(:) *dtlsm
             csite%co2budget_plresp(ipa) = csite%co2budget_plresp(ipa) +  &
                  plant_respiration * dtlsm
             csite%co2budget_rh(ipa) = csite%co2budget_rh(ipa) + csite%rh(ipa) * dtlsm
             cgrid%cbudget_nep(ipy) = cgrid%cbudget_nep(ipy) + cpoly%area(isi) * csite%area(ipa) *   &
                  dtlsm * (gpp - plant_respiration - csite%rh(ipa)) * umol_2_kgC


             ! Calculate transfers of sensible, latent
             ! -------------------------------------------------------
             call integrate_patch_ar(csite,   &
                  ipa,                     &
                  isi,                     &
                  ipy,                     &
                  ifm,                     &
                  integration_buff,        &
                  cpoly%met(isi)%rhos,     &
                  cpoly%met(isi)%vels,     &
                  cpoly%met(isi)%atm_tmp,  &
                  cpoly%met(isi)%atm_shv,  &
                  cpoly%met(isi)%atm_co2,  &
                  cpoly%met(isi)%geoht,    &
                  cpoly%met(isi)%exner,    &
                  cpoly%met(isi)%pcpg,     &
                  cpoly%met(isi)%qpcpg,    &
                  cpoly%met(isi)%dpcpg,    &
                  cpoly%met(isi)%prss,     &
                  cpoly%met(isi)%atm_shv,  &
                  cpoly%met(isi)%geoht,    &
                  cpoly%lsl(isi))

             ! Update the minimum monthly temperature, based on canopy temperature
!             if ( cpoly%min_monthly_temp(isi) < cpoly%site(isi)%can_tmp(ipa) ) &
!                  cpoly%min_monthly_temp(isi)=cpoly%site(isi)%can_tmp(ipa)


          end do patchloop
                    
       end do siteloop

    end do polygonloop

    return
  end subroutine rk4_timestep_ar

!==============================================================

  subroutine integrate_patch_ar(csite,ipa,isi,ipy,ifm, integration_buff, rhos,  &
       vels, atm_tmp, rv, atm_co2, zoff, exner, pcpg, qpcpg, dpcpg, prss,  &
       atm_shv, geoht, lsl)

    use ed_state_vars,only:sitetype,patchtype,integration_vars_ar,rk4patchtype

    use misc_coms, only: dtlsm
    use soil_coms, only: soil_rough
    use consts_coms, only: vonk, cp
    use rk4_coms, only: tbeg,tend,dtrk4,dtrk4i

    implicit none

    type(sitetype),target   :: csite
    integer :: ifm,ipy,isi,ipa
    integer, intent(in) :: lsl
    type(integration_vars_ar), target :: integration_buff

    real, intent(in) :: rhos
    real, intent(in) :: vels
    real, intent(in) :: atm_tmp
    real, intent(in) :: atm_shv
    real, intent(in) :: rv
    real, intent(in) :: atm_co2
    real, intent(in) :: zoff
    real, intent(in) :: exner
    real, intent(in) :: pcpg
    real, intent(in) :: qpcpg
    real, intent(in) :: dpcpg
    real, intent(in) :: prss
    real, intent(in) :: geoht

    type(rk4patchtype), pointer :: initp
    real :: hbeg
    real :: factv
    real :: aux
    real, parameter :: exar=2.5
    real :: zveg
    real :: zdisp
    real, parameter :: snowrough=0.001
    
    logical, save :: first_time=.true.
    
    
    if (first_time) then
       first_time = .false.
       tbeg   = 0.0
       tend   = dtlsm
       dtrk4  = tend - tbeg
       dtrk4i = 1./dtrk4
    end if

    !---------------------------------
    ! Set up the integration patch
    !---------------------------------

    initp => integration_buff%initp

    call copy_patch_init_ar(csite,ipa, initp, lsl)

    ! initial step size.  experience has shown that giving this too large a 
    ! value causes the integrator to fail (e.g., soil layers become
    ! supersaturated).
    hbeg = csite%htry(ipa)


    ! This is Bob Walko's recommended way of calculating the resistance.
    ! Note that temperature, not potential temperature, is input here.

    initp%rough = max(soil_rough, csite%veg_rough(ipa)) * (1.0 - csite%snowfac(ipa)) +  &
         snowrough
    zveg = csite%veg_height(ipa) * (1.0 - csite%snowfac(ipa))
    zdisp = 0.63 * zveg

    ! Zero the canopy-atmosphere flux values.  These values are updated
    ! every dtlsm, so they must be zeroed at each call.
    ! -----------------------------------------------------------------

    initp%upwp = 0.
    initp%tpwp = 0.
    initp%rpwp = 0.
    initp%wpwp = 0.


    call ed_stars(atm_tmp, atm_shv, atm_co2, initp%can_temp, geoht, vels,  &
         initp%rough, initp%ustar, initp%rstar, initp%tstar, initp%cstar,  &
         initp%can_shv, initp%can_co2)

    initp%tstar = initp%tstar * cp / exner  

!    if (ipa.eq.1) then
!       print*,"AR ",atm_tmp, atm_shv, initp%can_temp, geoht, vels,  &
!         initp%rough, initp%ustar, initp%rstar, initp%tstar,  &
!         initp%can_shv
!    endif

    factv = 1.0 / (vonk * initp%ustar)
    aux = exp(exar * (1. - (zdisp + initp%rough) / zveg))
    initp%rasveg = factv * zveg / (exar * (zveg - zdisp)) * (exp(exar) - aux)
    
    !aux = exp(0.925-1.575*csite%veg_rough/zdisp)
    !initp%rasveg = 1.081 * log(geoht/csite%veg_rough) * (12.182 - aux)  &
    !     / (vonk**2 * vels) 


    ! Go into the integrator
    call odeint_ar(hbeg, csite,ipa,isi,ipy,ifm,  &
         integration_buff, rhos, vels, atm_tmp, atm_shv, atm_co2, geoht,  &
         exner, pcpg, qpcpg, dpcpg, prss, lsl)

    ! Normalize canopy-atmosphere flux values.  These values are updated
    ! every dtlsm, so they must be normalized every time.
    ! -----------------------------------------------------------------
    initp%upwp = rhos*initp%upwp/dtlsm
    initp%tpwp = rhos*initp%tpwp/dtlsm
    initp%rpwp = rhos*initp%rpwp/dtlsm
    initp%wpwp = rhos*initp%wpwp/dtlsm
    
    
    !------------------------
    ! Move the state variables from the integrated patch to the model patch
    !------------------------
    call initp2modelp_ar(tend-tbeg, initp, csite, ipa,isi,ipy, rhos, lsl)

    return
  end subroutine integrate_patch_ar

  !====================================================================

  subroutine initp2modelp_ar(hdid, initp, csite, ipa,isi,ipy, rhos, lsl)

    use ed_state_vars,only:sitetype,patchtype,rk4patchtype,edgrid_g
    use consts_coms, only: day_sec,t3ple
    use ed_misc_coms,only: fast_diagnostics
    use soil_coms, only: soil, slz, min_sfcwater_mass
    use ed_misc_coms,only: fast_diagnostics
    use grid_coms, only: nzg, nzs
    use canopy_radiation_coms, only: lai_min, veg_temp_min
    use therm_lib, only: qwtk
    use ed_therm_lib, only: calc_hcapveg,ed_grndvap
    implicit none

    integer, intent(in) :: lsl
    real, intent(in) :: rhos

    type(sitetype),target :: csite
    type(patchtype),pointer :: cpatch
    type(rk4patchtype), target :: initp

    integer :: ipa,ico,ipy,isi
    integer :: k,ksn,nsoil, nlsw1
    real :: hdid
    real :: available_water
    real :: surface_temp, surface_fliq
    real, parameter :: tendays_sec=10.*day_sec
    real :: fracliq

    csite%can_temp(ipa) = initp%can_temp
    csite%can_shv(ipa) = initp%can_shv
    csite%can_co2(ipa) = initp%can_co2

    csite%ustar(ipa) = initp%ustar
    csite%tstar(ipa) = initp%tstar
    csite%rstar(ipa) = initp%rstar
    csite%cstar(ipa) = initp%cstar

    csite%upwp(ipa) = initp%upwp
    csite%wpwp(ipa) = initp%wpwp
    csite%tpwp(ipa) = initp%tpwp
    csite%rpwp(ipa) = initp%rpwp

    if(fast_diagnostics) then
       
       csite%wbudget_loss2atm(ipa) = initp%wbudget_loss2atm
       csite%ebudget_loss2atm(ipa) = initp%ebudget_loss2atm
       csite%co2budget_loss2atm(ipa) = initp%co2budget_loss2atm
       csite%ebudget_latent(ipa) = initp%ebudget_latent
       csite%nlev_sfcwater(ipa) = initp%nlev_sfcwater
       !    csite%avg_gpp(ipa) = initp%avg_gpp
       csite%avg_vapor_vc(ipa)       = initp%avg_vapor_vc  
       csite%avg_dew_cg(ipa)         = initp%avg_dew_cg    
       csite%avg_vapor_gc(ipa)       = initp%avg_vapor_gc  
       csite%avg_wshed_vg(ipa)       = initp%avg_wshed_vg  
       csite%avg_vapor_ac(ipa)       = initp%avg_vapor_ac
       csite%avg_transp(ipa)         = initp%avg_transp
       csite%avg_evap(ipa)           = initp%avg_evap
       csite%avg_netrad(ipa)         = initp%avg_netrad
       csite%aux(ipa)                = initp%aux
       csite%avg_sensible_vc(ipa)    = initp%avg_sensible_vc  
       csite%avg_sensible_2cas(ipa)  = initp%avg_sensible_2cas
       csite%avg_qwshed_vg(ipa)      = initp%avg_qwshed_vg    
       csite%avg_sensible_gc(ipa)    = initp%avg_sensible_gc  
       csite%avg_sensible_ac(ipa)    = initp%avg_sensible_ac  
       csite%avg_sensible_tot(ipa)   = initp%avg_sensible_tot
       csite%avg_carbon_ac(ipa)      = initp%avg_carbon_ac

       do k = lsl, nzg
          csite%avg_sensible_gg(k,ipa) = initp%avg_sensible_gg(k)
          csite%avg_smoist_gg(k,ipa)   = initp%avg_smoist_gg(k)
          csite%avg_smoist_gc(k,ipa)   = initp%avg_smoist_gc(k)
          csite%aux_s(k,ipa)  = initp%aux_s(k)
       enddo

    endif

       
    ! The following is not a pure diagnostic, it is used for phenology
    ! and mortality functions, preserve this variable and its dependencies
    ! in all contexts

    csite%avg_daily_temp(ipa) = csite%avg_daily_temp(ipa) + csite%can_temp(ipa)

    ! [KIM - 10-day average of plant available water - paw_avg10d
    ! MLO - Added after the return statement to avoid computations over water.
    !      Changed the name from theta to available_water
    
    cpatch => csite%patch(ipa)

    do ico = 1,cpatch%ncohorts

       available_water = 0.0
       do k = cpatch%krdepth(ico), nzg - 1
          available_water = available_water                               &
               + real((initp%soil_water(k)-dble(soil(csite%ntext_soil(k,ipa))%soilcp))        &
               * dble(slz(k+1)-slz(k))/dble(soil(csite%ntext_soil(k,ipa))%slmsts &
               - soil(csite%ntext_soil(k,ipa))%soilcp) )
       enddo
       available_water = available_water + real((initp%soil_water(nzg)  &
            -dble(soil(csite%ntext_soil(nzg,ipa))%soilcp)) &
            *dble(-1.0*slz(nzg))              &
            /dble(soil(csite%ntext_soil(nzg,ipa))%slmsts &
            -soil(csite%ntext_soil(nzg,ipa))%soilcp)) 
       available_water = available_water/(-1.0*slz(cpatch%krdepth(ico)))
       cpatch%paw_avg10d(ico) = cpatch%paw_avg10d(ico)*(1.0-hdid/tendays_sec)  &
            + available_water*hdid/tendays_sec

    enddo

    
    do k = lsl, nzg
       csite%soil_water(k,ipa) = initp%soil_water(k)
       csite%soil_energy(k,ipa) = initp%soil_energy(k)
       csite%soil_tempk(k,ipa) = initp%soil_tempk(k)
       csite%soil_fracliq(k,ipa) = initp%soil_fracliq(k)
    enddo
    

    do k = 1, csite%nlev_sfcwater(ipa)
       !-----------------------------------------------------------------------------------!
       !    Surface water energy is computed in J/m² inside the integrator. Converting it  !
       ! back to J/kg.                                                                     !
       !-----------------------------------------------------------------------------------!
       if (initp%sfcwater_mass(k) > min_sfcwater_mass) then
           csite%sfcwater_depth(k,ipa)   = initp%sfcwater_depth(k)
           csite%sfcwater_mass(k,ipa)    = initp%sfcwater_mass(k)
           csite%sfcwater_tempk(k,ipa)   = initp%sfcwater_tempk(k)
           csite%sfcwater_fracliq(k,ipa) = initp%sfcwater_fracliq(k)
           csite%sfcwater_energy(k,ipa)  = initp%sfcwater_energy(k)/initp%sfcwater_mass(k)
       elseif (k == 1) then
          csite%sfcwater_energy(k,ipa)  = 0.
          csite%sfcwater_mass(k,ipa)    = 0.
          csite%sfcwater_depth(k,ipa)   = 0.
          csite%sfcwater_fracliq(k,ipa) = csite%soil_fracliq(nzg,ipa)
          csite%sfcwater_tempk(k,ipa)   = csite%soil_tempk(nzg,ipa)
       else
          csite%sfcwater_energy(k,ipa)  = 0.
          csite%sfcwater_mass(k,ipa)    = 0.
          csite%sfcwater_depth(k,ipa)   = 0.
          csite%sfcwater_fracliq(k,ipa) = csite%sfcwater_fracliq(k-1,ipa)
          csite%sfcwater_tempk(k,ipa)   = csite%sfcwater_tempk(k-1,ipa)
       end if
    enddo
    

    do ico = 1,cpatch%ncohorts

       cpatch%veg_water(ico)  = initp%veg_water(ico)
       cpatch%veg_energy(ico) = initp%veg_energy(ico)

       ! For plants with minimal foliage, fix the vegetation
       ! temperature to the canopy air space
       if (cpatch%lai(ico) < lai_min) then

          cpatch%veg_temp(ico)   = csite%can_temp(ipa)
          cpatch%veg_water(ico)  = 0. 
          cpatch%veg_energy(ico) = cpatch%hcapveg(ico) * cpatch%veg_temp(ico)
       end if

       call qwtk(cpatch%veg_energy(ico),cpatch%veg_water(ico),cpatch%hcapveg(ico)          &
                ,cpatch%veg_temp(ico),fracliq)
        
       if (cpatch%veg_temp(ico) < veg_temp_min .or. cpatch%veg_temp(ico) > 360.0  ) then
          print*,"==========================================================="
          
          write(unit=*,fmt='(a,1x,es14.7)') 'veg temp      :',cpatch%veg_temp(ico)
          write(unit=*,fmt='(a,1x,es14.7)') 'veg_temp min  :',veg_temp_min
          write(unit=*,fmt='(a,1x,es14.7)') 'lai           :',cpatch%lai(ico)
          write(unit=*,fmt='(a,1x,es14.7)') 'lai_min       :', lai_min
          print*,"veg_water",cpatch%veg_water(ico),cpatch%veg_water(ico)/cpatch%lai(ico)
          print*,"veg_energy",cpatch%veg_energy(ico),cpatch%veg_energy(ico)/cpatch%lai(ico)
          print*,"hcapveg",cpatch%hcapveg(ico),cpatch%hcapveg(ico)/cpatch%lai(ico)

          print*,"Polygon:",ipy," Site:",isi," Patch:",ipa," Cohort:",ico," of",cpatch%ncohorts
          
          print*,"LAI",cpatch%lai(ico)
          print*,"Height",cpatch%hite(ico)
          print*,"DBH",cpatch%dbh(ico)
          print*,"Phenology Status",cpatch%phenology_status(ico)
          print*,"PFT",cpatch%pft(ico)
          print*,"Leaf Biomass",cpatch%bleaf(ico)
          print*,"Density",cpatch%nplant(ico)
          print*,"Patch LAI",csite%lai(ipa)
          print*,"Patch Disturbance Type",csite%dist_type(ipa)
          print*,"Patch Canopy Temperature",csite%can_temp(ipa)
          
          write(unit=*,fmt='(a,1x,i5)')     '================== FATAL ERROR =================='
          write(unit=*,fmt='(a,1x,i5)')     ' IPY:',ipy
          write(unit=*,fmt='(a,1x,i5)')     ' ISI:',isi
          write(unit=*,fmt='(a,1x,i5)')     ' IPA:',ipa
          write(unit=*,fmt='(a,1x,i5)')     ' ICO:',ico
!!$       write(unit=*,fmt='(a,1x,f14.5)')  ' Longitude:',edgrid_g(1)%lon(ipy)
!!$       write(unit=*,fmt='(a,1x,f14.5)')  ' Latitude: ',edgrid_g(1)%lat(ipy)
!!$       write(unit=*,fmt='(a)')           ' '
!!$       write(unit=*,fmt='(a,1x,es14.7)') ' PRSS:     ',edgrid_g(1)%polygon(ipy)%met(isi)%prss
!!$       write(unit=*,fmt='(a,1x,es14.7)') ' ATM_TMP:  ',edgrid_g(1)%polygon(ipy)%met(isi)%atm_tmp
!!$       write(unit=*,fmt='(a,1x,es14.7)') ' RHOS:     ',edgrid_g(1)%polygon(ipy)%met(isi)%rhos
!!$       write(unit=*,fmt='(a,1x,es14.7)') ' PCPG:     ',edgrid_g(1)%polygon(ipy)%met(isi)%pcpg
!!$       write(unit=*,fmt='(a,1x,es14.7)') ' SHV:g/kg  ',edgrid_g(1)%polygon(ipy)%met(isi)%atm_shv*1000
!!$       write(unit=*,fmt='(a,1x,es14.7)') ' VELS:     ',edgrid_g(1)%polygon(ipy)%met(isi)%vels
!!$       write(unit=*,fmt='(a)')           ' '
!!$       write(unit=*,fmt='(a,1x,es14.7)') ' rshort_v: ',cpatch%rshort_v(ico)
!!$       write(unit=*,fmt='(a,1x,es14.7)') ' rlong_v:  ',cpatch%rlong_v(ico)
!!$       write(unit=*,fmt='(a)')           ' '
!!$       write(unit=*,fmt='(a,1x,es14.7)') ' can_temp :',csite%can_temp(ipa)
!!$       write(unit=*,fmt='(a,1x,es14.7)') ' can_shv g/kg :',csite%can_shv(ipa)*1000.
!!$       write(unit=*,fmt='(a,1x,es14.7)') ' gnd_shv g/kg :',csite%ground_shv(ipa)*1000.
!!$       write(unit=*,fmt='(a)')           ' '
!!$       write(unit=*,fmt='(a,1x,es14.7)') 'Lai_coh   :',cpatch%lai(ico)
!!$       
!!$       write(unit=*,fmt='(a,1x,es14.7)') 'veg_water :',cpatch%veg_water(ico)
!!$       write(unit=*,fmt='(a,1x,es14.7)') 'rb        :',cpatch%rb(ico)
!!$       write(unit=*,fmt='(a)')           ' '
!!$
          call print_patch_ar(initp, csite,ipa, lsl)
          
          call fatal_error('extreme vegetation temperature','initp2modelp','rk4_driver.f90')

       end if
    end do


    ksn = csite%nlev_sfcwater(ipa)
    nsoil = csite%ntext_soil(nzg,ipa)
    nlsw1 = max(1, ksn)

    call ed_grndvap(ksn, nsoil, csite%soil_water(nzg,ipa),   &
         csite%soil_energy(nzg,ipa), csite%sfcwater_energy(nlsw1,ipa), &
         rhos, &
         csite%can_shv(ipa),csite%ground_shv(ipa),csite%surface_ssh(ipa), &
         surface_temp,surface_fliq)

    return
  end subroutine initp2modelp_ar

  !******************************************************************************

  subroutine canopy_atm_fluxes_ar(csite,cpoly,ipa,isi)

    use ed_state_vars,only:polygontype,sitetype,patchtype
    use consts_coms, only: cpi

    implicit none

    type(polygontype),target :: cpoly
    type(sitetype),target :: csite
    integer :: ipa,isi

    print*,'decide how to set vels in canopy_atm_fluxes.'
    stop

    ! Calculate turbulent fluxes between atmosphere and canopy
    !    pis = cpoly%pi0 * cpi
    !    thetacan = pss%can_temp / pis
    
    !    if(thetacan.lt.cpoly%theta)then
    !       cpoly%vels = cpoly%vels_stab
    !    else
    !       cpoly%vels = cpoly%vels_unstab
    !    endif

    return
  end subroutine canopy_atm_fluxes_ar
  
end module rk4_driver_ar

!================================================================

real function compute_water_storage_ar(csite, lsl, rhos,ipa)

  use ed_state_vars,only:sitetype,patchtype
  use grid_coms, only: nzg
  use soil_coms, only: dslz

  implicit none
  type(sitetype), target :: csite
  type(patchtype),pointer :: cpatch
  integer :: k
  integer :: ipa,ico
  integer, intent(in) :: lsl
  real, intent(in) :: rhos

  compute_water_storage_ar = 0.0
  cpatch => csite%patch(ipa)

  do k = lsl, nzg
     compute_water_storage_ar = compute_water_storage_ar + &
          real(csite%soil_water(k,ipa)) * dslz(k) * 1000.0
  enddo

  do k = 1, csite%nlev_sfcwater(ipa)
     compute_water_storage_ar = compute_water_storage_ar +  &
          csite%sfcwater_mass(k,ipa)
  enddo
  
  compute_water_storage_ar = compute_water_storage_ar +  &
       csite%can_shv(ipa) * csite%veg_height(ipa) * rhos

  do ico = 1,cpatch%ncohorts
     compute_water_storage_ar = compute_water_storage_ar +  &
          cpatch%veg_water(ico)
     
  enddo

  return
end function compute_water_storage_ar
!================================================================

real function compute_netrad_ar(csite,ipa)

  use ed_state_vars,only:sitetype,patchtype

  implicit none

  type(sitetype), target :: csite
  type(patchtype),pointer :: cpatch
  integer :: ipa,ico
  integer :: k

  cpatch => csite%patch(ipa)

  compute_netrad_ar = 0.0
  compute_netrad_ar = csite%rshort_g(ipa) + csite%rlong_g(ipa) + csite%rlong_s(ipa)

  do k = 1, csite%nlev_sfcwater(ipa)
     compute_netrad_ar = compute_netrad_ar + csite%rshort_s(k,ipa)
  enddo

  do ico = 1,cpatch%ncohorts
     compute_netrad_ar = compute_netrad_ar + cpatch%rshort_v(ico) + cpatch%rlong_v(ico)
  enddo

  return
end function compute_netrad_ar

!=====================================================================

real function compute_energy_storage_ar(csite, lsl, rhos, ipa)

  use ed_state_vars,only:sitetype,patchtype
  use grid_coms, only: nzg
  use soil_coms, only: dslz
  use consts_coms, only: cp, cliq, cice, alli, t3ple
  use canopy_radiation_coms, only: lai_min

  implicit none
  
  type(sitetype),target :: csite
  type(patchtype),pointer :: cpatch
  integer :: ipa,ico
  real, intent(in) :: rhos
  integer, intent(in) :: lsl
  integer :: k
  real :: soil_storage, sfcwater_storage, cas_storage, veg_storage

  cpatch => csite%patch(ipa)

  soil_storage = 0.0
  do k = lsl, nzg
     soil_storage = soil_storage + csite%soil_energy(k,ipa) * dslz(k)
  enddo

  sfcwater_storage = 0.0
  do k = 1, csite%nlev_sfcwater(ipa)
     sfcwater_storage = sfcwater_storage + csite%sfcwater_energy(k,ipa) *   &
          csite%sfcwater_mass(k,ipa)
  enddo

  cas_storage = cp * rhos * csite%veg_height(ipa) * csite%can_temp(ipa)

  veg_storage = 0.0
  do ico = 1,cpatch%ncohorts

     !!!!!! ASSUMES THAT THE TALLEST COHORT IS IN BIN 1 !!!!!!!!!

     if(csite%lai(ipa) > lai_min)then
        ! MLO: I think we can simply add the vegetation energy instead of recalculating
        !      it from the temperature. This will avoid problems when veg_temp is t3ple
        !      and water and ice coexist.
        veg_storage = veg_storage + cpatch%veg_energy(ico)
     end if

  end do

  compute_energy_storage_ar = soil_storage + sfcwater_storage + cas_storage + &
       veg_storage

  return
end function compute_energy_storage_ar
!=====================================================================

subroutine sum_plant_cfluxes_ar(csite,ipa, gpp, gpp_dbh,plresp)

  use ed_state_vars,only:sitetype,patchtype
  use consts_coms, only: day_sec, umol_2_kgC
  use canopy_radiation_coms, only: lai_min
  use max_dims, only: n_dbh

  implicit none

  type(sitetype),target :: csite
  type(patchtype),pointer :: cpatch
  integer :: ipa,ico,idbh
  real, intent(out) :: gpp
  real, dimension(n_dbh), intent(out) :: gpp_dbh
  real, intent(out) :: plresp
  real, parameter             :: ddbh=1./real(n_dbh-1)
  real :: lrresp ! leaf and root respiration
  real :: sresp ! storage, growth, vleaf respiration
  logical :: forest
  
  !GPP by DBH is computed for forested areas only.
  forest = csite%dist_type(ipa) /= 1

  gpp = 0.0
  gpp_dbh= 0.0 
  lrresp = 0.0
  sresp = 0.0
  cpatch => csite%patch(ipa)

  do ico = 1,cpatch%ncohorts

     if(cpatch%lai(ico) > lai_min)then
        gpp = gpp + cpatch%gpp(ico)
        if (forest) then 
           idbh=max(1,min(n_dbh,ceiling(cpatch%dbh(ico)*ddbh)))
           gpp_dbh(idbh) = gpp_dbh(idbh) + cpatch%gpp(ico)
        end if

        lrresp = lrresp + cpatch%leaf_respiration(ico)

     endif
     lrresp = lrresp + cpatch%root_respiration(ico)
     sresp = sresp + (cpatch%growth_respiration(ico) + cpatch%storage_respiration(ico)   &
          + cpatch%vleaf_respiration(ico)) * cpatch%nplant(ico) / (day_sec * umol_2_kgC)
  enddo
  plresp = lrresp + sresp

  return
end subroutine sum_plant_cfluxes_ar


!=====================================================================

real function compute_co2_storage_ar(csite, rhos, ipa)

  use ed_state_vars,only: sitetype
  use consts_coms, only : mmdryi
  implicit none
  type(sitetype),target :: csite
  real, intent(in) :: rhos
  integer :: ipa
  compute_co2_storage_ar = csite%can_co2(ipa) * mmdryi * rhos * csite%veg_height(ipa)
  return
end function compute_co2_storage_ar


