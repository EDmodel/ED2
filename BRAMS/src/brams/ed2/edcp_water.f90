!==========================================================================================!
!==========================================================================================!
!     This subroutine will compute the fluxes between water and the air.  ED solves only   !
! the land regions, so this subroutine is a simplified version of LEAF-3 without consider- !
! ing land, only water.  It's called lake model but it is used for oceans and rivers too   !
! (using prescribed surface temperature.)                                                  !
!------------------------------------------------------------------------------------------!
subroutine simple_lake_model(time,dtlongest)
   use node_mod               , only : ja                & ! intent(in)
                                     , jz                & ! intent(in)
                                     , ia                & ! intent(in)
                                     , iz                & ! intent(in)
                                     , mynum             ! ! intent(in)
   use consts_coms            , only : stefan            & ! intent(in)
                                     , cpi               & ! intent(in)
                                     , vonk              & ! intent(in)
                                     , cp                & ! intent(in)
                                     , grav              & ! intent(in)
                                     , rdry              & ! intent(in)
                                     , t00               & ! intent(in)
                                     , p00               & ! intent(in)
                                     , p00i              & ! intent(in)
                                     , rocp              & ! intent(in)
                                     , cpor              & ! intent(in)
                                     , alvl              & ! intent(in)
                                     , mmdryi            & ! intent(in)
                                     , mmdry             ! ! intent(in)
   use canopy_air_coms        , only : ustmin            & ! intent(in)
                                     , ubmin             ! ! intent(in)
   use mem_edcp               , only : wgrid_g           ! ! structure
   use io_params              , only : ssttime1          & ! intent(in)
                                     , ssttime2          & ! intent(in)
                                     , iupdsst           ! ! intent(in)
   use mem_leaf               , only : leaf_g            & ! structure
                                     , dtleaf            ! ! intent(in)
   use mem_basic              , only : co2_on            & ! intent(in)
                                     , co2con            & ! intent(in)
                                     , basic_g           ! ! structure
   use mem_radiate            , only : radiate_g         ! ! structure
   use mem_cuparm             , only : cuparm_g          ! ! structure
   use mem_micro              , only : micro_g           ! ! structure
   use mem_grid               , only : zt                & ! intent(in)
                                     , grid_g            & ! structure
                                     , dzt               & ! intent(in)
                                     , zm                & ! intent(in)
                                     , if_adap           & ! intent(in)
                                     , jdim              & ! intent(in)
                                     , ngrid             & ! intent(in)
                                     , nzs               ! ! intent(in)
   use met_driver_coms        , only : atm_tmp_min       & ! intent(in)
                                     , atm_tmp_max       & ! intent(in)
                                     , atm_shv_min       & ! intent(in)
                                     , atm_shv_max       & ! intent(in)
                                     , prss_min          & ! intent(in)
                                     , prss_max          ! ! intent(in)
   use leaf_coms              , only : min_waterrough    & ! intent(in)
                                     , can_depth_min     ! ! intent(in)
   use rk4_coms               , only : rk4min_can_shv    & ! intent(in)
                                     , rk4max_can_shv    & ! intent(in)
                                     , rk4min_can_temp   & ! intent(in)
                                     , rk4max_can_temp   & ! intent(in)
                                     , rk4min_can_co2    & ! intent(in)
                                     , rk4max_can_co2    ! ! intent(in)
   use therm_lib              , only : rhovsil           & ! function
                                     , reducedpress      & ! function
                                     , ptqz2enthalpy     & ! function
                                     , hpqz2temp         & ! function
                                     , idealdenssh       ! ! function
   use canopy_struct_dynamics , only : ed_stars          & ! subroutine
                                     , vertical_vel_flux ! ! function
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real        , intent(in) :: dtlongest
   real(kind=8), intent(in) :: time
   !----- Local variables -----------------------------------------------------------------!
   integer                  :: k1w,k2w,k3w
   integer                  :: k2u,k3u,k2u_1,k3u_1
   integer                  :: k2v,k3v,k2v_1,k3v_1
   integer                  :: i,j,k,n
   integer                  :: niter_leaf
   real                     :: topma_t,wtw,wtu1,wtu2,wtv1,wtv2
   real                     :: dtll, dtll_factor,dtllowcc,dtlloccc
   real                     :: up_mean,vp_mean,exner_mean, rv_mean,rtp_mean
   real                     :: theta_mean,co2p_mean
   real                     :: cosz, geoht, vels
   real                     :: cosine, sine
   real                     :: atm_temp,atm_enthalpy,atm_theta,atm_shv,atm_co2,atm_prss
   real                     :: can_temp,can_enthalpy,can_theta,can_shv,can_co2,can_prss
   real                     :: can_exner,atm_rhos,can_rhos
   real                     :: water_temp,water_ssh,water_rough
   real                     :: ustar,tstar,estar,qstar,cstar
   real(kind=8)             :: angle
   real                     :: water_lc,water_ac
   real                     :: sensible_lc,enthalpy_ac
   real                     :: carbon_emission,carbon_uptake,carbon_ac
   real                     :: avg_water_lc,avg_water_ac
   real                     :: avg_sensible_lc,avg_enthalpy_ac
   real                     :: avg_carbon_emission,avg_carbon_uptake,avg_carbon_ac
   real                     :: prev_can_shv,prev_can_temp,prev_can_enthalpy,prev_can_co2
   real                     :: prev_can_rhos
   real                     :: hcapcan, wcapcan, ccapcan, rdi
   real                     :: sflux_u,sflux_v,sflux_w,sflux_t,sflux_r,sflux_c
   real                     :: gzotheta, timefac_sst
   real                     :: fm
   !----- Local constants -----------------------------------------------------------------!
   real, parameter          :: d0=0.
   real, parameter          :: z0fac_water = .016 / grav
   real, parameter          :: emiss_w = 0.97         ! emissivity of water
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Define lake time-split timesteps here.  This ensures that the lake model will    !
   ! never use a timestep longer than about 30 seconds, but the actual number depends on   !
   ! the user's choice and the actual time step.                                           !
   !---------------------------------------------------------------------------------------!
   niter_leaf  = max(1,nint(dtlongest/dtleaf+.4))
   dtll_factor = 1. / float(niter_leaf)
   dtll        = dtlongest * dtll_factor
   
   
   
   !---------------------------------------------------------------------------------------!
   !     If we have time-varying water surface temperature, we perform a linear interpol-  !
   ! ation over time now. Otherwise, we just use the initially prescribed value (which     !
   ! can be spatially heterogeneous.                                                       !
   !---------------------------------------------------------------------------------------!
   if (iupdsst == 0) then
      timefac_sst = 0. !----- This will force it to use the initial value only. -----------!
   else
      timefac_sst = sngl((time - ssttime1(ngrid)) / (ssttime2(ngrid) - ssttime1(ngrid)))
   end if



   !---------------------------------------------------------------------------------------!
   !    Big domain loops start here.                                                       !
   !---------------------------------------------------------------------------------------!
   jloop: do j=ja,jz
      iloop: do i=ia,iz

         !---------------------------------------------------------------------------------!
         !   We now retrieve some variables from the atmosphere.  We may need to do some   !
         ! interpolation in the vertical, for example, when we use adaptive coordinate.    !
         ! Therefore, we always check this.                                                !
         !---------------------------------------------------------------------------------!
         select case (if_adap)
         case (0)
            !------------------------------------------------------------------------------!
            !      Terrain-following coordinate (sigma_z).  Here we simply use the lowest  !
            ! predicted level for thermo and wind variables, and a simple average between  !
            ! levels one and two for Exner function and density.  I am not sure this is    !
            ! consistent with what is done for density in adaptive coordinate...           !
            !------------------------------------------------------------------------------!
            !----- W/T points, only the i,j points are needed...  -------------------------!
            k2w   = nint(grid_g(ngrid)%flpw(i,j))
            k1w   = k2w - 1
            k3w   = k2w + 1
            theta_mean = basic_g(ngrid)%theta(2,i,j)
            rv_mean    = basic_g(ngrid)%rv(2,i,j)
            rtp_mean   = basic_g(ngrid)%rtp(2,i,j)

            if (co2_on) then
               co2p_mean = basic_g(ngrid)%co2p(2,i,j)
            else 
               co2p_mean = co2con(1)
            end if

            up_mean    = (basic_g(ngrid)%up(2,i,j) + basic_g(ngrid)%up(2,i-1,j))     * 0.5
            vp_mean    = (basic_g(ngrid)%vp(2,i,j) + basic_g(ngrid)% vp(2,i,j-jdim)) * 0.5
            exner_mean = ( basic_g(ngrid)%pp(1,i,j) + basic_g(ngrid)%pp(2,i,j)             &
                         + basic_g(ngrid)%pi0(1,i,j)+basic_g(ngrid)%pi0(2,i,j))      * 0.5
         case (1)
            !------------------------------------------------------------------------------!
            !     Adaptive height coordinate (shaved-eta).  Here we need to do a weighted  !
            ! average using the lowest predicted points and the points avove them.         !
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !    Finding the lowest boundary levels, depending on which kind of variable   !
            ! we are averaging (staggered grid).                                           !
            !------------------------------------------------------------------------------!
            !----- U points, need to average between i-1 and i... -------------------------!
            k2u   = nint(grid_g(ngrid)%flpu(i,j))
            k3u   = k2u + 1
            k2u_1 = nint(grid_g(ngrid)%flpu(i-1,j))
            k3u_1 = k2u_1 + 1
            !----- V points, need to average between j-jdim and j... ----------------------!
            k2v   = nint(grid_g(ngrid)%flpv(i,j))
            k3v   = k2v+1
            k2v_1 = nint(grid_g(ngrid)%flpv(i,j-jdim))
            k3v_1 = k2v_1 + 1
            !----- W/T points, only the i,j points are needed...  -------------------------!
            k2w   = nint(grid_g(ngrid)%flpw(i,j))
            k1w   = k2w - 1
            k3w   = k2w + 1
            !------------------------------------------------------------------------------!
            topma_t = .25 * ( grid_g(ngrid)%topma(i,j) + grid_g(ngrid)%topma(i-1,j)            &
                            + grid_g(ngrid)%topma(i,j-jdim) + grid_g(ngrid)%topma(i-1,j-jdim))

            !------------------------------------------------------------------------------!
            !     Computing the weights for lowest predicted points, relative to points    !
            ! above them.                                                                  !
            !------------------------------------------------------------------------------!
            wtw  = (zm(k2w) - topma_t) * dzt(k2w)
            wtu1 = grid_g(ngrid)%aru(k2u_1,i-1,j)    / grid_g(ngrid)%aru(k3u_1,i-1,j)
            wtu2 = grid_g(ngrid)%aru(k2u,i,j)        / grid_g(ngrid)%aru(k3u,i,j)
            wtv1 = grid_g(ngrid)%arv(k2v_1,i,j-jdim) / grid_g(ngrid)%arv(k3v_1,i,j-jdim)
            wtv2 = grid_g(ngrid)%arv(k2v,i,j)        / grid_g(ngrid)%arv(k3v,i,j)

            theta_mean   =       wtw   * basic_g(ngrid)%theta(k2w,i,j)                     &
                         + (1. - wtw)  * basic_g(ngrid)%theta(k3w,i,j)

            rv_mean      =       wtw   * basic_g(ngrid)%rv(k2w,i,j)                        &
                         + (1. - wtw)  * basic_g(ngrid)%rv(k3w,i,j)

            rtp_mean     =       wtw   * basic_g(ngrid)%rtp(k2w,i,j)                       &
                         + (1. - wtw)  * basic_g(ngrid)%rtp(k3w,i,j)

            if (co2_on) then
               co2p_mean =       wtw   * basic_g(ngrid)%co2p(k2w,i,j)                       &
                         + (1. - wtw)  * basic_g(ngrid)%co2p(k3w,i,j)

            else 
               co2p_mean = co2con(1)
            end if

            up_mean      = (        wtu1  * basic_g(ngrid)%up(k2u_1,i-1,j)                 &
                           +  (1. - wtu1) * basic_g(ngrid)%up(k3u_1,i-1,j)                 &
                           +        wtu2  * basic_g(ngrid)%up(k2u,i,j)                     &
                           +  (1. - wtu2) * basic_g(ngrid)%up(k3u,i,j)       ) * .5
            vp_mean      = (        wtv1  * basic_g(ngrid)%vp(k2v_1,i,j-jdim)              &
                           +  (1. - wtv1) * basic_g(ngrid)%vp(k3v_1,i,j-jdim)              &
                           +        wtv2  * basic_g(ngrid)%vp(k2v,i,j)                     &
                           +  (1. - wtv2) * basic_g(ngrid)%vp(k3v,i,j)       ) * .5
            
            !------------------------------------------------------------------------------!
            !     Exner function and density, need to consider the layer beneath the       !
            ! ground to make the profile.                                                  !
            !------------------------------------------------------------------------------!
            if (wtw >= .5) then
               exner_mean  = (wtw - .5)                                                    &
                           * (basic_g(ngrid)%pp(k1w,i,j) + basic_g(ngrid)%pi0(k1w,i,j))    &
                           + (1.5 - wtw)                                                   &
                           * (basic_g(ngrid)%pp(k2w,i,j) + basic_g(ngrid)%pi0(k2w,i,j))
            else
               exner_mean  = (wtw + .5)                                                    &
                           * (basic_g(ngrid)%pp(k2w,i,j) + basic_g(ngrid)%pi0(k2w,i,j))    &
                           + (.5 - wtw)                                                    &
                           * (basic_g(ngrid)%pp(k3w,i,j) + basic_g(ngrid)%pi0(k3w,i,j))
            end if
         end select

         !----- Defining some variables that won't change during the sub time steps. ------!
         cosz         = radiate_g(ngrid)%cosz(i,j)
         geoht        = (zt(k2w)-zm(k1w))*grid_g(ngrid)%rtgt(i,j)
         atm_prss     = (exner_mean * cpi) ** cpor * p00
         atm_temp     = cpi * theta_mean * exner_mean
         atm_theta    = theta_mean
         atm_co2      = co2p_mean
         !---------------------------------------------------------------------------------!
         !   Most of ED expects specific humidity, not mixing ratio.  Since we will use    !
         ! ed_stars, which is set up for the former, not the latter, we locally solve      !
         ! everything for specific humidity, converting in the end.                        !
         !---------------------------------------------------------------------------------!
         atm_shv  = rv_mean / (1. + rtp_mean)
         !---------------------------------------------------------------------------------!

         !----- Finding the atmospheric enthalpy and density. -----------------------------!
         !if (atm_temp < atm_tmp_min .or. atm_temp > atm_tmp_max  .or.                      &
         !    atm_shv  < atm_shv_min .or. atm_shv  > atm_shv_max  .or.                      &
         !    atm_prss < prss_min    .or. atm_prss > prss_max   ) then
         !    write (unit=*,fmt='(a)') '======== Weird atm properties... ========'
         !    write (unit=*,fmt='(a,i5)')   'Node             = ',mynum
         !    write (unit=*,fmt='(a,i5)')   'X                = ',i
         !    write (unit=*,fmt='(a,i5)')   'Y                = ',j
         !    write (unit=*,fmt='(a,f8.2)') 'LONG      [degE] = ',grid_g(ngrid)%glon(i,j)
         !    write (unit=*,fmt='(a,f8.2)') 'LAT       [degN] = ',grid_g(ngrid)%glat(i,j)
         !    write (unit=*,fmt='(a,f7.2)') 'ATM_PRSS  [ hPa] = ',atm_prss    * 0.01
         !    write (unit=*,fmt='(a,f7.2)') 'ATM_TEMP  [degC] = ',atm_temp - t00
         !    write (unit=*,fmt='(a,f7.2)') 'ATM_SHV   [g/kg] = ',atm_shv  * 1.e3
         !    call fatal_error('Non-sense atm met values!!!'                                &
         !                    ,'simple_lake_model','edcp_water.f90')
         !end if

         atm_enthalpy = ptqz2enthalpy(atm_prss,atm_temp,atm_shv,geoht)
         atm_rhos     = idealdenssh(atm_prss,atm_temp,atm_shv)


         !---------------------------------------------------------------------------------!
         !    Finding the mean wind speed and mean wind direction, so we can split the     !
         ! flux of momentum according to this.                                             !
         !---------------------------------------------------------------------------------!
         vels     = max(ubmin,sqrt(up_mean*up_mean + vp_mean*vp_mean))
         angle    = datan2(dble(vp_mean),dble(up_mean))
         cosine   = sngl(dcos(angle))
         sine     = sngl(dsin(angle))
 

         !---------------------------------------------------------------------------------!
         !     Assigning initial values for "canopy" properties.                           !
         !---------------------------------------------------------------------------------!
         !------ First we copy those variables that do not change with pressure. ----------!
         can_co2      = leaf_g(ngrid)%can_co2(i,j,1)
         can_theta    = leaf_g(ngrid)%can_theta(i,j,1)
         !------ Converting mixing ratio to specific humidity. ----------------------------!
         can_shv      = leaf_g(ngrid)%can_rvap(i,j,1) / (1. +leaf_g(ngrid)%can_rvap(i,j,1))

         !------ Finding the derived properties. ------------------------------------------!
         can_prss     = reducedpress(atm_prss,atm_theta,atm_shv,geoht                      &
                                             ,can_theta,can_shv,can_depth_min)
         can_exner    = cp * (can_prss * p00i) ** rocp
         can_temp     = cpi * can_theta * can_exner
         can_enthalpy = ptqz2enthalpy(can_prss,can_temp,can_shv,can_depth_min)
         can_rhos     = idealdenssh(can_prss,can_temp,can_shv)
         !---------------------------------------------------------------------------------!

         !------ Resetting the surface fluxes. --------------------------------------------!
         sflux_u    = 0.0
         sflux_v    = 0.0
         sflux_w    = 0.0
         sflux_r    = 0.0
         sflux_c    = 0.0
         sflux_t    = 0.0

         !------ Resetting the flux components. -------------------------------------------!
         avg_water_lc        = 0.0
         avg_water_ac        = 0.0
         avg_sensible_lc     = 0.0
         avg_enthalpy_ac     = 0.0
         avg_carbon_emission = 0.0
         avg_carbon_uptake   = 0.0

         !---------------------------------------------------------------------------------!
         !     Assigning an initial value for ustar.  At the first time, ustar is probably !
         ! zero, so to avoid singularities, we impose the minimum ustmin, which is also    !
         ! used in ed_stars subroutine.                                                    !
         !---------------------------------------------------------------------------------!
         ustar      = max(ustmin,wgrid_g(ngrid)%ustar(i,j))

         !---------------------------------------------------------------------------------!
         !     Loop for the intermediate time steps.                                       !
         !---------------------------------------------------------------------------------!
         timeloop: do n = 1,niter_leaf

            !------------------------------------------------------------------------------!
            !     Finding standard values for vapour and heat capacity of our layer.  This !
            ! assumes a constant thickness canopy.                                         !
            !------------------------------------------------------------------------------!
            wcapcan = can_depth_min * can_rhos
            hcapcan = wcapcan * cp
            ccapcan = wcapcan * mmdryi
            
            !----- Finding effective water surface roughness. -----------------------------!
            water_rough = max(z0fac_water * ustar ** 2,min_waterrough)
            
            !------------------------------------------------------------------------------!
            !    Finding the current water temperature and saturation spec. humidity.      !
            !------------------------------------------------------------------------------!
            water_temp = leaf_g(ngrid)%seatp(i,j)                                          &
                       + timefac_sst*(leaf_g(ngrid)%seatf(i,j) - leaf_g(ngrid)%seatp(i,j))
            water_ssh  = rhovsil(water_temp) / can_rhos
            !---------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            !     This is the ED-2 stars subroutine (Euler method, single precision).      !
            !------------------------------------------------------------------------------!
            call ed_stars(atm_theta,atm_enthalpy,atm_shv,atm_co2                           &
                         ,can_theta,can_enthalpy,can_shv,can_co2                           &
                         ,geoht,d0,vels,water_rough                                        &
                         ,ustar,tstar,estar,qstar,cstar,fm)

            
            gzotheta = grav * geoht / atm_theta
            sflux_w  = sflux_w + vertical_vel_flux(gzotheta,tstar,ustar)

            !------------------------------------------------------------------------------!
            !      Finding the resistance, and updating the canopy temperature and         !
            ! specific humidity.                                                           !
            !------------------------------------------------------------------------------!
            rdi = .2 * ustar

            !------ Copying the values to scratch buffers. --------------------------------!
            prev_can_enthalpy = can_enthalpy
            prev_can_temp     = can_temp
            prev_can_shv      = can_shv
            prev_can_co2      = can_co2

            !------ Computing the lake->canopy and free atmosphere->canopy fluxes. --------!
            water_lc        = can_rhos * ( water_ssh - can_shv) * rdi
            water_ac        = can_rhos * ustar * qstar
            sensible_lc     = can_rhos * cp * (water_temp -  can_temp) * rdi 
            enthalpy_ac     = can_rhos * ustar * estar
            carbon_emission = 0. ! For the time being, no idea of what to put in here...
            carbon_uptake   = 0. ! For the time being, no idea of what to put in here...
            carbon_ac       = can_rhos * ustar * cstar * mmdryi

            !----- Finding the scales for fluxes. -----------------------------------------!
            dtllowcc = dtll / (can_depth_min * can_rhos)
            dtlloccc = mmdry * dtllowcc

            !------------------------------------------------------------------------------!
            !    Updating the prognostic variables, enthalpy, specific humidity, and CO2   !
            ! mixing ratio.                                                                !
            !------------------------------------------------------------------------------!
            can_enthalpy = prev_can_enthalpy                                               &
                         + dtllowcc * ( sensible_lc + enthalpy_ac + alvl * avg_water_lc)
            can_shv      = prev_can_shv  + dtllowcc * (water_lc    + water_ac)
            can_co2      = prev_can_co2                                                    &
                         + dtlloccc * (carbon_emission + carbon_uptake  + carbon_ac)


            !------------------------------------------------------------------------------!
            !     Updating the thermodynamic properties.  Notice that we assume pressure   !
            ! to remain constant over the full time step, which is a simplification...     !
            !------------------------------------------------------------------------------!
            can_temp     = hpqz2temp(can_enthalpy,can_prss,can_shv,can_depth_min)
            can_theta    = cp * can_temp / can_exner
            can_rhos     = idealdenssh(can_prss,can_temp,can_shv)




            !----- Sanity check -----------------------------------------------------------!
            if(can_shv  < sngl(rk4min_can_shv)  .or. can_shv  > sngl(rk4max_can_shv)  .or. &
               can_temp < sngl(rk4min_can_temp) .or. can_temp > sngl(rk4max_can_temp) .or. &
               can_co2  < sngl(rk4min_can_co2)  .or. can_co2  > sngl(rk4max_can_co2) ) then

               write(unit=*,fmt='(a)') '====== PROBLEMS IN THE LAKE MODEL!!! ====== '
               write(unit=*,fmt='(3(a,1x,i5,1x))') 'i=',i,'j=',j,'n=',n
               write(unit=*,fmt='(2(a,1x,f8.2,1x))') 'Lon: ',grid_g(ngrid)%glon(i,j) &
                                                    ,'Lat: ',grid_g(ngrid)%glat(i,j)
               write(unit=*,fmt='(a,1x,es12.5)') 'EXNER (PIO)        : ',exner_mean
               write(unit=*,fmt='(a,1x,es12.5)') 'ATM_PRSS           : ',atm_prss
               write(unit=*,fmt='(a,1x,es12.5)') 'ATM_DENS           : ',atm_rhos
               write(unit=*,fmt='(a,1x,es12.5)') 'ATM_TEMP           : ',atm_temp
               write(unit=*,fmt='(a,1x,es12.5)') 'ATM_QVAP           : ',atm_shv
               write(unit=*,fmt='(a,1x,es12.5)') 'ATM_CO2            : ',atm_co2
               write(unit=*,fmt='(a,1x,es12.5)') 'WATER_TEMP         : ',water_temp
               write(unit=*,fmt='(a,1x,es12.5)') 'WATER_SSH          : ',water_ssh
               write(unit=*,fmt='(a,1x,es12.5)') 'CANOPY_PRSS        : ',can_prss
               write(unit=*,fmt='(a,1x,es12.5)') 'CANOPY_DENS        : ',can_rhos
               write(unit=*,fmt='(a,1x,es12.5)') 'CANOPY_ENTHALPY    : ',can_enthalpy
               write(unit=*,fmt='(a,1x,es12.5)') 'PREV_CAN_ENTHALPY  : ',prev_can_enthalpy
               write(unit=*,fmt='(a,1x,es12.5)') 'CANOPY_TEMP        : ',can_temp
               write(unit=*,fmt='(a,1x,es12.5)') 'PREV_CAN_TEMP      : ',prev_can_temp
               write(unit=*,fmt='(a,1x,es12.5)') 'CANOPY_QVAP        : ',can_shv
               write(unit=*,fmt='(a,1x,es12.5)') 'PREV_CAN_QVAP      : ',prev_can_shv
               write(unit=*,fmt='(a,1x,es12.5)') 'CANOPY_CO2         : ',can_co2
               write(unit=*,fmt='(a,1x,es12.5)') 'PREV_CAN_CO2       : ',prev_can_co2
               write(unit=*,fmt='(a,1x,es12.5)') 'RDI                : ',rdi
               write(unit=*,fmt='(a,1x,es12.5)') 'DTLLOWCC           : ',dtllowcc
               write(unit=*,fmt='(a,1x,es12.5)') 'USTAR              : ',ustar
               write(unit=*,fmt='(a,1x,es12.5)') 'TSTAR              : ',tstar
               write(unit=*,fmt='(a,1x,es12.5)') 'QSTAR              : ',qstar
               write(unit=*,fmt='(a,1x,es12.5)') 'CSTAR              : ',cstar
               write(unit=*,fmt='(a,1x,es12.5)') 'ESTAR              : ',estar
               write(unit=*,fmt='(a,1x,es12.5)') 'WATER_LAKE2CAN     : ',water_lc
               write(unit=*,fmt='(a,1x,es12.5)') 'WATER_ATM2CAN      : ',water_ac
               write(unit=*,fmt='(a,1x,es12.5)') 'HEAT_LAKE2CAN      : ',sensible_lc
               write(unit=*,fmt='(a,1x,es12.5)') 'HEAT_ATM2CAN       : ',enthalpy_ac
               write(unit=*,fmt='(a,1x,es12.5)') 'CARBON_EMISSION    : ',carbon_emission
               write(unit=*,fmt='(a,1x,es12.5)') 'CARBON_UPTAKE      : ',carbon_uptake
               write(unit=*,fmt='(a,1x,es12.5)') 'CARBON_ATM2CAN     : ',carbon_ac
               call fatal_error('Lake model failed','simple_lake_model','edcp_water.f90')
            end if

            !------------------------------------------------------------------------------!
            !     Integrate the full fluxes of water, heat, and momentum.  They will be    !
            ! normalised outside the loop.                                                 !
            !------------------------------------------------------------------------------!
            sflux_u = sflux_u - ustar*ustar*cosine
            sflux_v = sflux_v - ustar*ustar*sine
            sflux_t = sflux_t - ustar*tstar
            sflux_r = sflux_r                                                              &
                    - ustar * qstar / ((1.-atm_shv)*(1. - 0.5*(can_shv+prev_can_shv)))
            sflux_c = sflux_c - ustar * cstar

            !------------------------------------------------------------------------------!
            !     Also integrate some of the energy/water/carbon budget components.        !
            !------------------------------------------------------------------------------!
            avg_water_lc        = avg_water_lc        + water_lc
            avg_water_ac        = avg_water_ac        + water_ac
            avg_sensible_lc     = avg_sensible_lc     + sensible_lc
            avg_enthalpy_ac     = avg_enthalpy_ac     + enthalpy_ac
            avg_carbon_emission = avg_carbon_emission + carbon_emission
            avg_carbon_uptake   = avg_carbon_uptake   + carbon_uptake
         end do timeloop

         !---------------------------------------------------------------------------------!
         !     Transfer model scalars back to global arrays.                               !
         !---------------------------------------------------------------------------------!
         !----- Stars, converting qstar to rstar ------------------------------------------!
         wgrid_g(ngrid)%ustar(i,j) = ustar
         wgrid_g(ngrid)%rstar(i,j) = qstar                                                 &
                                   / ((1.-atm_shv)*(1. - 0.5*(can_shv+prev_can_shv)))
         wgrid_g(ngrid)%tstar(i,j) = tstar
         wgrid_g(ngrid)%cstar(i,j) = cstar



         !----- Finding the fluxes towards the free atmosphere. ---------------------------!
         wgrid_g(ngrid)%sflux_u(i,j) = atm_rhos * sflux_u * dtll_factor
         wgrid_g(ngrid)%sflux_v(i,j) = atm_rhos * sflux_v * dtll_factor
         wgrid_g(ngrid)%sflux_w(i,j) = atm_rhos * sflux_w * dtll_factor
         wgrid_g(ngrid)%sflux_t(i,j) = atm_rhos * sflux_t * dtll_factor
         wgrid_g(ngrid)%sflux_r(i,j) = atm_rhos * sflux_r * dtll_factor
         wgrid_g(ngrid)%sflux_c(i,j) = atm_rhos * sflux_c * dtll_factor
         
         !----- Finding the sources and sinks components. ---------------------------------!
         leaf_g(ngrid)%sensible(i,j,1) = avg_sensible_lc     * dtll_factor
         leaf_g(ngrid)%evap(i,j,1)     = avg_water_lc * alvl * dtll_factor
         leaf_g(ngrid)%transp(i,j,1)   = 0.
         leaf_g(ngrid)%gpp(i,j,1)      = avg_carbon_uptake   * dtll_factor
         leaf_g(ngrid)%plresp(i,j,1)   = 0.
         leaf_g(ngrid)%resphet(i,j,1)  = avg_carbon_emission * dtll_factor

         !----- Finding some radiation properties. ----------------------------------------!
         wgrid_g(ngrid)%albedt(i,j)  = min(max(-.0139 + .0467 * tan(acos(cosz)),.03),.999)
         wgrid_g(ngrid)%rlongup(i,j) = emiss_w * stefan * water_temp**4

         !----- Finding some canopy air properties. ---------------------------------------!
         leaf_g(ngrid)%can_theta(i,j,1)    = can_theta
         leaf_g(ngrid)%can_rvap(i,j,1)     = can_shv / (1. - can_shv)
         leaf_g(ngrid)%can_co2(i,j,1)      = can_co2
         leaf_g(ngrid)%can_prss(i,j,1)     = can_prss


         !---------------------------------------------------------------------------------! 
         !     Unless one includes a PFT for water lilies, there is no vegetation on       !
         ! water...                                                                        !
         !---------------------------------------------------------------------------------! 
         leaf_g(ngrid)%veg_energy(i,j,1) = 0.
         leaf_g(ngrid)%veg_water (i,j,1) = 0.

         !---------------------------------------------------------------------------------!
         !    Our simple lake model allows no ice, making sure the values are properly     !
         ! zeroed.                                                                         !
         !---------------------------------------------------------------------------------!
         leaf_g(ngrid)%sfcwater_nlev     (i,j,1) = 0.
         do k=1, nzs
            leaf_g(ngrid)%sfcwater_energy (1,i,j,1) = 0.
            leaf_g(ngrid)%sfcwater_mass   (1,i,j,1) = 0.
            leaf_g(ngrid)%sfcwater_depth  (1,i,j,1) = 0.
         end do
      end do iloop
   end do jloop

   return
end subroutine simple_lake_model
!==========================================================================================!
!==========================================================================================!
