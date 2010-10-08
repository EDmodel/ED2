!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   LEAF-3 wrapper, this will call the main driver for each grid.                          !
!------------------------------------------------------------------------------------------!
subroutine sfclyr(mzp,mxp,myp,ia,iz,ja,jz,ibcon)
   use mem_all
   use teb_spm_start  , only : TEB_SPM    ! ! intent(in)
   USE mem_teb        , only : teb_g      & ! data type
                             , teb_vars   ! ! type
   USE mem_teb_common , only : tebc_g     & ! data type
                             , teb_common ! ! type
   implicit none
   
   !----- Arguments -----------------------------------------------------------------------!
   integer                 , intent(in) :: mzp,mxp,myp,ia,iz,ja,jz,ibcon
   !----- Local variables -----------------------------------------------------------------!
   integer                              :: ng
   real, dimension(mxp,myp)             :: l_thils2,l_ths2,l_rvs2,l_rtps2,l_pis2,l_dens2
   real, dimension(mxp,myp)             :: l_ups2,l_vps2,l_zts2,l_co2s2
   !---------------------------------------------------------------------------------------!
   
   if (nstbot == 0) return

   ng=ngrid


   !----- Calling LEAF-3 main driver ------------------------------------------------------!
   call leaf3(mzp,mxp,myp,nzg,nzs,npatch,ia,iz,ja,jz,leaf_g (ng),basic_g (ng),turb_g (ng)  &
             ,radiate_g(ng),grid_g (ng),cuparm_g(ng),micro_g(ng),l_thils2,l_ths2,l_rvs2    &
             ,l_rtps2,l_co2s2,l_pis2,l_dens2,l_ups2,l_vps2,l_zts2,teb_g(ng),tebc_g(ng))

   !---- Calling TOPMODEL if the user wants so --------------------------------------------!
   if (isfcl == 2) then
      call hydro(mxp,myp,nzg,nzs,npatch            , leaf_g(ng)%soil_water                 &
                , leaf_g(ng)%soil_energy           , leaf_g(ng)%soil_text                  &
                , leaf_g(ng)%sfcwater_mass         , leaf_g(ng)%sfcwater_energy            &
                , leaf_g(ng)%patch_area            , leaf_g(ng)%patch_wetind               )
   end if
   
   !----- Apply lateral boundary conditions to leaf3 arrays -------------------------------!
   call leaf_bcond(mxp,myp,nzg,nzs,npatch,jdim     , leaf_g(ng)%soil_water                 &
                  , leaf_g(ng)%sfcwater_mass       , leaf_g(ng)%soil_energy                &
                  , leaf_g(ng)%sfcwater_energy     , leaf_g(ng)%soil_text                  &
                  , leaf_g(ng)%sfcwater_depth      , leaf_g(ng)%ustar                      &
                  , leaf_g(ng)%tstar               , leaf_g(ng)%rstar                      &
                  , leaf_g(ng)%cstar               , leaf_g(ng)%zeta                       &
                  , leaf_g(ng)%ribulk              , leaf_g(ng)%veg_albedo                 &
                  , leaf_g(ng)%veg_fracarea        , leaf_g(ng)%veg_lai                    &
                  , leaf_g(ng)%veg_tai             , leaf_g(ng)%veg_rough                  &
                  , leaf_g(ng)%veg_height          , leaf_g(ng)%patch_area                 &
                  , leaf_g(ng)%patch_rough         , leaf_g(ng)%patch_wetind               &
                  , leaf_g(ng)%leaf_class          , leaf_g(ng)%soil_rough                 &
                  , leaf_g(ng)%sfcwater_nlev       , leaf_g(ng)%stom_resist                &
                  , leaf_g(ng)%ground_rsat         , leaf_g(ng)%ground_rvap                &
                  , leaf_g(ng)%ground_temp         , leaf_g(ng)%ground_fliq                &
                  , leaf_g(ng)%veg_water           , leaf_g(ng)%veg_hcap                   &
                  , leaf_g(ng)%veg_energy          , leaf_g(ng)%can_prss                   &
                  , leaf_g(ng)%can_theiv           , leaf_g(ng)%can_theta                  &
                  , leaf_g(ng)%can_rvap            , leaf_g(ng)%can_co2                    &
                  , leaf_g(ng)%sensible            , leaf_g(ng)%evap                       &
                  , leaf_g(ng)%transp              , leaf_g(ng)%gpp                        &
                  , leaf_g(ng)%plresp              , leaf_g(ng)%resphet                    &
                  , leaf_g(ng)%veg_ndvip           , leaf_g(ng)%veg_ndvic                  &
                  , leaf_g(ng)%veg_ndvif           )
   return
end subroutine sfclyr
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This is the LEAF-3 main driver.                                                       !
!------------------------------------------------------------------------------------------!
subroutine leaf3(m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz,leaf,basic,turb,radiate,grid,cuparm,micro &
                ,thils2,ths2,rvs2,rtps2,co2s2,pis2,dens2,ups2,vps2,zts2,teb,tebc)

   use mem_all
   use leaf_coms
   use rconstants
   use node_mod       , only : mynum         ! ! intent(in)
   use catt_start     , only : CATT          ! ! intent(in)
   use teb_spm_start  , only : TEB_SPM       ! ! intent(in)
   use mem_teb        , only : teb_vars      ! ! type
   use mem_teb_common , only : teb_common    ! ! type
   use therm_lib      , only : rslif         & ! function
                             , tslif         & ! function
                             , reducedpress  & ! function
                             , thetaeiv      & ! function
                             , thetaeiv2thil & ! function
                             , idealdenssh   & ! function
                             , qwtk          ! ! subroutine
   use mem_scratch    , only : scratch
   implicit none
   
   !----- Arguments -----------------------------------------------------------------------!
   integer               , intent(in)    :: m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz
   type(leaf_vars)                       :: leaf
   type(basic_vars)                      :: basic
   type(turb_vars)                       :: turb
   type(radiate_vars)                    :: radiate
   type(grid_vars)                       :: grid
   type(cuparm_vars)                     :: cuparm
   type(micro_vars)                      :: micro
   real, dimension(m2,m3), intent(inout) :: thils2,ths2,rvs2,rtps2,co2s2,pis2,dens2
   real, dimension(m2,m3), intent(inout) :: ups2,vps2,zts2
   !----- Local variables -----------------------------------------------------------------!
   type(teb_vars)                      :: teb
   type(teb_common)                    :: tebc
   integer                             :: i,j,k,ksn,ip,iter_leaf,nveg
   real                                :: dvelu,dvelv,velnew,sflux_uv,cosine1,sine1
   real                                :: psup1,psup2,depe,alt2,deze,dpdz,exn1st,airt
   real                                :: zh_town,zle_town,zsfu_town,zsfv_town
   real                                :: sfcw_energy_int
   !----- Saved variables -----------------------------------------------------------------!
   real                  , save        :: emis_town, alb_town, ts_town,g_urban
   logical               , save        :: firsttime = .true. 
   !---------------------------------------------------------------------------------------!

   !----- Assigning some TEB variables to zero. This will remain so if TEB is not used. ---!
   if (firsttime) then 
      G_URBAN    = 0.
      EMIS_TOWN  = 0.
      ALB_TOWN   = 0.
      TS_TOWN    = 0.
      firsttime  = .false.
   end if

   !----- Time interpolation factor for updating SST. -------------------------------------!
   if (iupdsst == 0) then
      timefac_sst = 0.
   else
      timefac_sst = sngl((time - ssttime1(ngrid)) / (ssttime2(ngrid) - ssttime1(ngrid)))
   endif

   !---------------------------------------------------------------------------------------!
   !    Define LEAF-3 time-split timesteps here.  This ensures that LEAF-3 will never use  !
   ! a timestep longer than about 30 seconds, but the actual number depends on the user's  !
   ! choice and the actual time step.                                                      !
   !---------------------------------------------------------------------------------------!
   niter_leaf  = max(1,nint(dtlt/dtleaf + .4))
   dtll_factor = 1. / float(niter_leaf)
   dtll        = dtlt * dtll_factor

   !----- Check whether we have CO2, and copy to an scratch array. ------------------------!
   if (co2_on) then
      call atob(m1*m2*m3,basic%co2p,scratch%vt3do)
   else
      call ae0(m1*m2*m3,scratch%vt3do,co2con(1))
   end if

   !----- Copy surface atmospheric variables into 2d arrays for input to LEAF. ------------!
   select case (if_adap)
   case (0)
      call sfc_fields(m1,m2,m3,ia,iz,ja,jz,jdim           , basic%thp     , basic%theta    &
                          , basic%rv      , basic%rtp     , scratch%vt3do , basic%up       &
                          , basic%vp      , basic%dn0     , basic%pp      , basic%pi0      &
                          , grid%rtgt     , zt,zm,thils2,ths2,rvs2,rtps2,co2s2,ups2,vps2   &
                          ,pis2,dens2,zts2)
   case (1)
      call sfc_fields_adap(m1,m2,m3,ia,iz,ja,jz,jdim      , grid%flpu     , grid%flpv      &
                          , grid%flpw     , grid%topma    , grid%aru      , grid%arv       &
                          , basic%thp     , basic%theta   , basic%rv      , basic%rtp      &
                          , scratch%vt3do , basic%up      , basic%vp      , basic%dn0      &
                          , basic%pp      , basic%pi0     , zt,zm,dzt,thils2,ths2,rvs2     &
                          ,rtps2,co2s2,ups2,vps2,pis2,dens2,zts2)
   end select

   !----- Big domain loop -----------------------------------------------------------------!
   jloop1: do j = ja,jz
      iloop1: do i = ia,iz

         !----- Copy surface variables to single-column values ----------------------------!
         atm_up       = ups2(i,j)
         atm_vp       = vps2(i,j)
         atm_thil     = thils2(i,j)
         atm_theta    = ths2(i,j)
         atm_rvap     = rvs2(i,j)
         atm_rtot     = rtps2(i,j)
         atm_shv      = atm_rvap / (1. + atm_rvap)
         geoht        = zts2(i,j)
         atm_exner    = pis2(i,j)
         atm_co2      = co2s2(i,j)
         atm_prss     = p00 * (cpi * atm_exner) ** cpor
         vels         = sqrt(atm_up ** 2 + atm_vp ** 2)
         atm_temp     = cpi * atm_theta * atm_exner
         atm_theiv    = thetaeiv(atm_thil,atm_prss,atm_temp,atm_rvap,atm_rtot,-8)

         !----- Update water internal energy from time-dependent SST. ---------------------!
         leaf%soil_energy(mzg,i,j,1) =  cliq * (leaf%seatp(i,j)                            &
                                     + (leaf%seatf(i,j) - leaf%seatp(i,j)) * timefac_sst   &
                                     - tsupercool)

         !----- Fill surface precipitation arrays for input to LEAF-3 ---------------------!
         call sfc_pcp(nnqparm(ngrid),i,j,cuparm,micro)

         !---------------------------------------------------------------------------------!
         !    Zero out albedo, upward surface longwave, and momentum, heat, and moisture   !
         ! flux arrays before summing over patches.                                        !
         !---------------------------------------------------------------------------------!
         if (ilwrtyp > 0 .or. iswrtyp > 0) then
            radiate%albedt(i,j)  = 0.
            radiate%rlongup(i,j) = 0.
         end if

         !----- Resetting surface fluxes --------------------------------------------------!
         turb%sflux_u(i,j) = 0.
         turb%sflux_v(i,j) = 0.
         turb%sflux_w(i,j) = 0.
         turb%sflux_t(i,j) = 0.
         turb%sflux_r(i,j) = 0.
         turb%sflux_c(i,j) = 0.

         !----- For no soil model (patch 2) fill "canopy" temperature and moisture. -------!
         if (isfcl == 0) then
            leaf%can_theta (i,j,2) = atm_theta - dthcon
            leaf%can_rvap  (i,j,2) = atm_rvap  - drtcon
            leaf%can_co2   (i,j,2) = atm_co2

            can_shv                = leaf%can_rvap(i,j,2) / (1. + leaf%can_rvap(i,j,2))

            leaf%can_prss  (i,j,2) = reducedpress(atm_prss,atm_theta,atm_shv,geoht         &
                                                 ,leaf%can_theta(i,j,ip),can_shv           &
                                                 ,can_depth_min)
            can_exner              = cp  * (p00i * leaf%can_prss(i,j,2)) ** rocp
            can_temp               = cpi * leaf%can_theta(i,j,2) * can_exner

            leaf%can_theiv (i,j,2) = thetaeiv(leaf%can_theta(i,j,2),leaf%can_prss(i,j,2)   &
                                             ,can_temp,leaf%can_rvap(i,j,2)                &
                                             ,leaf%can_rvap(i,j,2),-8)
            can_lntheiv            = log(leaf%can_theiv (i,j,2))
            leaf%patch_area(i,j,1) = min(1.0,max(tiny_parea,1.-pctlcon))
            leaf%patch_area(i,j,2) = 1.0 - leaf%patch_area(i,j,1)
         end if

         !----- Begin patch loop. ---------------------------------------------------------!
         ploop1: do ip = 1,np

            !----- Resetting all fluxes. --------------------------------------------------!
            leaf%sensible(i,j,ip) = 0.
            leaf%evap    (i,j,ip) = 0.
            leaf%transp  (i,j,ip) = 0.
            leaf%gpp     (i,j,ip) = 0.
            leaf%plresp  (i,j,ip) = 0.
            leaf%resphet (i,j,ip) = 0.

            !----- Update time-dependent vegetation LAI and fractional coverage -----------!
            if (ip >= 2 .and. leaf%patch_area(i,j,ip) >= tiny_parea .and. isfcl >= 1) then
               call vegndvi( ngrid                      , leaf%patch_area  (i,j,ip)        &
                           , leaf%leaf_class(i,j,ip)    , leaf%veg_fracarea(i,j,ip)        &
                           , leaf%veg_lai   (i,j,ip)    , leaf%veg_tai     (i,j,ip)        &
                           , leaf%veg_rough (i,j,ip)    , leaf%veg_height  (i,j,ip)        &
                           , leaf%veg_albedo(i,j,ip)    , leaf%veg_ndvip   (i,j,ip)        &
                           , leaf%veg_ndvic (i,j,ip)    , leaf%veg_ndvif   (i,j,ip)        )
            end if

            !------------------------------------------------------------------------------!
            !    Converting the surface water internal energy to an extensive. We will in- !
            ! tegrate the extensive variable rather than the intensive, since the latter   !
            ! has a strong dependence on sfcwater_mass and this increases errors. The      !
            ! extensive variable is more accurate because the heat balance is done using   !
            ! extensive fluxes (W/m²), not (W/kg). After each iteration the intensive      !
            ! quantity is updated so routines like sfcrad will work fine.                  !
            !------------------------------------------------------------------------------!
            ksn = nint(leaf%sfcwater_nlev(i,j,ip))
            sfcw_energy_ext(:) = 0.
            do k=1,ksn
               sfcw_energy_ext(k) = leaf%sfcwater_energy(k,i,j,ip)                         &
                                  * leaf%sfcwater_mass(k,i,j,ip)
            end do

            !----- Set initial value of specific humidity ---------------------------------!
            can_shv     = leaf%can_rvap(i,j,ip) / (leaf%can_rvap(i,j,ip) + 1.)

            !------------------------------------------------------------------------------!
            !    Update canopy air pressure here.  Canopy air pressure is assumed to       !
            ! remain constant during one LEAF full timestep, which means that heat equals  !
            ! to change in enthalpy.  Between time steps, atmospheric pressure changes so  !
            ! enthalpy is no longer a conserved variable.  Potential temperature and       !
            ! equivalent potential temperature, on the other hand, are conserved if no     !
            ! heat happens, which is the case between successive calls.  Therefore, we     !
            ! track theta and theta_eiv.                                                   !
            !------------------------------------------------------------------------------!
            leaf%can_prss(i,j,ip) = reducedpress(atm_prss,atm_theta,atm_shv,geoht          &
                                                ,leaf%can_theta(i,j,ip),can_shv,can_depth)

            !----- Define the canopy depth. -----------------------------------------------!
            nveg      = nint(leaf%leaf_class(i,j,ip))
            if (nveg == 0 .or. ip == 1) then
               can_depth = can_depth_min
            else
               can_depth = max(can_depth_min,veg_ht(nveg))
            end if

            !----- Begin LEAF-3 small timestep here. --------------------------------------!
            tloop: do iter_leaf = 1,niter_leaf


               !---------------------------------------------------------------------------!
               !    Calculate radiative fluxes between atmosphere, vegetation, and ground/ !
               ! /snow based on already-computed downward shortwave and longwave fluxes    !
               ! from the atmosphere.  Fill tempk array with soil and snow temperature (C) !
               ! and fracliq array with liquid fraction of water content in soil and snow. !
               ! Other snowcover properties are also computed here.                        !
               !---------------------------------------------------------------------------!

               if (iswrtyp > 0 .or. ilwrtyp > 0) then

                  if (ip == 1 .or. leaf%patch_area(i,j,ip) >= tiny_parea) then

                     ! If TEB is on, copy the values to single variables ------------------!
                     if (teb_spm == 1) then
                        g_urban   = leaf%g_urban(i,j,ip)
                        emis_town = tebc%emis_town(i,j)
                        alb_town  = tebc%alb_town(i,j)
                        ts_town   = tebc%ts_town(i,j)
                     end if
                     
                     
                     !---------------------------------------------------------------------!
                     !     Find the radiation terms, and update canopy air temperature,    !
                     ! density, and enthalpy, and vegetation temperature and liquid water. !
                     ! fraction.                                                           !
                     !---------------------------------------------------------------------!
                     call sfcrad(mzg,mzs,ip              , leaf%soil_energy     (:,i,j,ip) &
                       , leaf%soil_water      (:,i,j,ip) , leaf%soil_text       (:,i,j,ip) &
                       , leaf%sfcwater_energy (:,i,j,ip) , leaf%sfcwater_mass   (:,i,j,ip) &
                       , leaf%sfcwater_depth  (:,i,j,ip) , leaf%patch_area      (  i,j,ip) &
                       , leaf%leaf_class      (  i,j,ip) , leaf%veg_height      (  i,j,ip) &
                       , leaf%veg_fracarea    (  i,j,ip) , leaf%veg_albedo      (  i,j,ip) &
                       , leaf%sfcwater_nlev   (  i,j,ip) , leaf%veg_energy      (  i,j,ip) &
                       , leaf%veg_water       (  i,j,ip) , leaf%veg_hcap        (  i,j,ip) &
                       , leaf%can_prss        (  i,j,ip) , leaf%can_theiv       (  i,j,ip) &
                       , leaf%can_theta       (  i,j,ip) , leaf%can_rvap        (  i,j,ip) &
                       , radiate%rshort       (  i,j   ) , radiate%rlong        (  i,j   ) &
                       , radiate%albedt       (  i,j   ) , radiate%rlongup      (  i,j   ) &
                       , radiate%cosz         (  i,j   ) , g_urban                         &
                       , emis_town                       , alb_town                        &
                       , ts_town                         )
                  else
                     !----- Update canopy air temperature, density, and theta. ------------!
                     can_lntheiv            = log(leaf%can_theiv(i,j,ip))
                     can_exner              = cp  * (p00i * leaf%can_prss(i,j,ip)) ** rocp
                     can_temp               = cpi * leaf%can_theta(i,j,ip) * can_exner
                     can_shv                = leaf%can_rvap(i,j,ip)                        &
                                            / (1. + leaf%can_rvap(i,j,ip))
                     can_rhos               = idealdenssh (leaf%can_prss(i,j,ip)           &
                                                          ,can_temp,can_shv)
                     !----- Find vegetation temperature and surface liquid water fraction. !
                     call qwtk(leaf%veg_energy(i,j,ip),leaf%veg_water(i,j,ip)              &
                              ,leaf%veg_hcap(i,j,ip),veg_temp,veg_fliq)
                  end if
               else
                  !----- Update canopy air temperature, density, and theta. ---------------!
                  can_lntheiv            = log(leaf%can_theiv(i,j,ip))
                  can_exner              = cp  * (p00i * leaf%can_prss(i,j,ip)) ** rocp
                  can_temp               = cpi * leaf%can_theta(i,j,ip) * can_exner
                  can_shv                = leaf%can_rvap(i,j,ip)                           &
                                         / (1. + leaf%can_rvap(i,j,ip))                
                  can_rhos               = idealdenssh (leaf%can_prss(i,j,ip)              &
                                                       ,can_temp,can_shv)              

                  !----- Find vegetation temperature and surface liquid water fraction. ---!
                  call qwtk(leaf%veg_energy(i,j,ip),leaf%veg_water(i,j,ip)                 &
                           ,leaf%veg_hcap(i,j,ip),veg_temp,veg_fliq)
               end if

               !---------------------------------------------------------------------------!
               !    For water surface (patch 1), compute surface saturation mixing ratio   !
               ! and roughness length based on previous ustar.  For soil patches, compute  !
               ! roughness length based on vegetation and snow.                            !
               !---------------------------------------------------------------------------!
               if (ip == 1) then

                  leaf%ground_temp(i,j,ip) = tempk(mzg)
                  leaf%ground_rsat(i,j,ip) = rslif(leaf%can_prss(i,j,ip),tempk(mzg))
                  leaf%ground_rvap(i,j,ip) = leaf%ground_rsat(i,j,ip)
                  if (tempk(mzg) >= t3ple) then
                     leaf%ground_fliq(i,j,ip) = 1.
                  else
                     leaf%ground_fliq(i,j,ip) = 0.
                  end if
                  leaf%patch_rough(i,j,ip) = max( z0fac_water * leaf%ustar(i,j,ip) ** 2    &
                                                , min_waterrough)

               elseif (isfcl >= 1 .and. leaf%patch_area(i,j,ip) >= tiny_parea) then
                  leaf%patch_rough(i,j,ip) = max( grid%topzo(i,j)                          &
                                                , leaf%soil_rough(i,j,ip)                  &
                                                , leaf%veg_rough(i,j,ip))                  &
                                           * (1. - snowfac) + snowrough * snowfac
               end if

               !---------------------------------------------------------------------------!
               !    Imposing minimum wind speed.                                           !
               !---------------------------------------------------------------------------!
               vels_pat = max(ubmin,vels)


               !----- Compute the characteristic scales. ----------------------------------!
               call leaf_stars(atm_theta,atm_theiv,atm_shv,atm_rvap,atm_co2                &
                              ,leaf%can_theta(i,j,ip),leaf%can_theiv(i,j,ip),can_shv       &
                              ,leaf%can_rvap(i,j,ip),leaf%can_co2(i,j,ip)                  &
                              ,geoht,vels_pat,dtll                                         &
                              ,leaf%patch_rough(i,j,ip),leaf%ustar(i,j,ip)                 &
                              ,leaf%tstar(i,j,ip),estar,qstar,leaf%rstar(i,j,ip)           &
                              ,leaf%cstar(i,j,ip),leaf%zeta(i,j,ip),leaf%ribulk(i,j,ip)    &
                              ,leaf%R_aer(i,j,ip))


               if (teb_spm==1) g_urban = leaf%g_urban(i,j,ip)

               call sfclmcv( leaf%ustar        (i,j,ip) , leaf%tstar              (i,j,ip) &
                           , leaf%rstar        (i,j,ip) , leaf%cstar              (i,j,ip) &
                           , leaf%zeta         (i,j,ip) , vels                             &
                           , vels_pat                   , atm_up                           &
                           , atm_vp                     , leaf%patch_area         (i,j,ip) &
                           , turb%sflux_u      (i,j   ) , turb%sflux_v            (i,j   ) &
                           , turb%sflux_w      (i,j   ) , turb%sflux_t            (i,j   ) &
                           , turb%sflux_r      (i,j   ) , turb%sflux_c            (i,j   ) &
                           , g_urban                    )

               select case (ip)
               case (1) !----- Water. -----------------------------------------------------!
                  !------------------------------------------------------------------------!
                  !     For water patches, update temperature and moisture of "canopy"     !
                  ! from divergence of fluxes with water surface and atmosphere.           !
                  ! rdi = rho * ustar/5 is the viscous sublayer conductivity from Garratt  !
                  ! (1992), but multiplied by density.                                     !
                  !------------------------------------------------------------------------!
                  rho_ustar = leaf%ustar(i,j,1) * can_rhos
                  rdi       = .2 * rho_ustar
                  dtllowcc  = dtll / (can_depth * can_rhos)
                  dtllohcc  = dtll / (can_depth * can_rhos * cp * can_temp)
                  dtlloccc  = mmdry * dtllowcc

                  !----- Compute the fluxes from water body to canopy. --------------------!
                  hflxgc  = rdi * cp * (leaf%ground_temp(i,j,ip) - can_temp)
                  wflxgc  = rdi * (leaf%ground_rsat(i,j,ip) - leaf%can_rvap(i,j,ip))
                  qwflxgc = wflxgc * alvl
                  cflxgc  = 0. !----- No water carbon emission model available...

                  !----- Compute the fluxes from atmosphere to canopy air space. ----------!
                  eflxac  = rho_ustar * estar              * cp * can_temp
                  wflxac  = rho_ustar * leaf%rstar(i,j,ip)
                  cflxac  = rho_ustar * leaf%cstar(i,j,ip) * mmdryi

                  !----- Integrate the state variables. -----------------------------------!
                  can_lntheiv            = can_lntheiv                                     &
                                         + dtllohcc * (hflxgc + qwflxgc + eflxac)

                  leaf%can_rvap(i,j,ip)  = leaf%can_rvap(i,j,ip)                           &
                                         + dtllowcc * (wflxgc + wflxac)

                  leaf%can_co2(i,j,ip)   = leaf%can_co2(i,j,ip)                            &
                                         + dtlloccc * (cflxgc + cflxac)

                  !----- Integrate the fluxes. --------------------------------------------!
                  leaf%sensible(i,j,ip) = leaf%sensible(i,j,ip) + hflxgc * dtll_factor
                  leaf%evap    (i,j,ip) = leaf%evap    (i,j,ip) + wflxgc * dtll_factor
                  leaf%plresp  (i,j,ip) = leaf%plresp  (i,j,ip) + cflxgc * dtll_factor

                  !------------------------------------------------------------------------!
                  !    Update canopy air properties.  This is done inside the loop to      !
                  ! ensure that we leave here with the right canopy air potential temper-  !
                  ! ature.                                                                 !
                  !------------------------------------------------------------------------!
                  leaf%can_theiv(i,j,ip) = exp(can_lntheiv)
                  leaf%can_theta(i,j,ip) = thetaeiv2thil(leaf%can_theiv(i,j,ip)            &
                                                        ,leaf%can_prss (i,j,ip)            &
                                                        ,leaf%can_rvap (i,j,ip))
                  can_shv                = leaf%can_rvap(i,j,ip)                           &
                                         / (leaf%can_rvap(i,j,ip) + 1.)
                  can_temp               = cpi * leaf%can_theta(i,j,ip) * can_exner
                  can_rhos               = idealdenssh(leaf%can_prss(i,j,ip)               &
                                                      ,can_temp,can_shv)

                  !----- KML - drydep -----------------------------------------------------!
                  if (catt == 1) leaf%R_aer(i,j,ip) = .2 * leaf%ustar(i,j,ip)


               case default !---- Land patch. ---------------------------------------------!

                  if (teb_spm == 1) then

                     !----- TEB: Urban canopy parameterization starts here ----------------!
                     if (nint(leaf%g_urban(i,j,ip))/=0.) then

                        !edmilson
                        psup1  = ((basic%pi0(1,i,j)+basic%pp(1,i,j))*cpi)**cpor*p00
                        psup2  = ((basic%pi0(2,i,j)+basic%pp(2,i,j))*cpi)**cpor*p00
                        depe   = psup2-psup1
                        alt2   = zt(1)*grid%rtgt(i,j)
                        deze   = geoht-alt2
                        dpdz   = depe/deze
                        exn1st = (psup2/p00)**rocp
                        
                        airt= basic%theta(2,i,j)*exn1st 

                        ! TEB - defining pointers
                        g_urban   = leaf%g_urban(i,j,ip)

                        call leaf3_teb_interface(istp,dtlt,dtll                            &
                           , radiate%cosz(i,j)          , geoht                            &
                           , radiate%rlong(i,j)         , radiate%rshort(i,j)              &
                           , psup2                      , airt                             &
                           , atm_up                     , atm_vp                           &
                           , basic%rv(2,i,j)            , pcpgl/dtlt                       &
                           , teb%fuso(i,j)              , teb%t_canyon(i,j)                &
                           , teb%r_canyon(i,j)          , teb%ts_roof(i,j)                 &
                           , teb%ts_road(i,j)           , teb%ts_wall(i,j)                 &
                           , teb%ti_road(i,j)           , teb%ti_bld(i,j)                  &
                           , teb%ws_roof(i,j)           , teb%ws_road(i,j)                 &
                           , teb%t_roof(2:4,i,j)        , teb%t_road(2:4,i,j)              &
                           , teb%t_wall(2:4,i,j)        , zh_town                          &
                           , zle_town                   , tebc%emis_town(i,j)              &
                           , zsfu_town                  , zsfv_town                        &
                           , tebc%ts_town(i,j)          , tebc%alb_town(i,j)               &
                           , nint(g_urban)              , teb%h_traffic(i,j)               &
                           , teb%h_industry(i,j)        , teb%le_traffic(i,j)              &
                           , teb%le_industry(i,j)       , teb%t2m_town(i,j)                &
                           , teb%r2m_town(i,j)          , time                             &
                           , itimea                     , dpdz                             &
                           , can_rhos                   )
                        
                        turb%sflux_u(i,j) = turb%sflux_u(i,j)                              &
                                          + leaf%patch_area(i,j,ip) * zsfu_town
                        turb%sflux_v(i,j) = turb%sflux_v(i,j)                              &
                                          + leaf%patch_area(i,j,ip) * zsfv_town
                        turb%sflux_t(i,j) = turb%sflux_t(i,j)                              &
                                          + leaf%patch_area(i,j,ip) * zh_town              &
                                          / (cp * can_rhos)
                        turb%sflux_r(i,j) = turb%sflux_r(i,j)                              &
                                          + leaf%patch_area(i,j,ip) * zle_town             &
                                          / (alvl * can_rhos)
                     end if
                  end if

                  if (isfcl >= 1) then
                     !---------------------------------------------------------------------!
                     !     For soil model patches, update temperature and moisture of      !
                     ! soil, vegetation, and canopy                                        !
                     !---------------------------------------------------------------------!
                     if (leaf%patch_area(i,j,ip) >= tiny_parea) then
                        call leaftw(mzg,mzs,np                                             &
                          , leaf%soil_water     (:,i,j,ip), leaf%soil_energy    (:,i,j,ip) &
                          , leaf%soil_text      (:,i,j,ip), leaf%sfcwater_mass  (:,i,j,ip) &
                          , leaf%sfcwater_depth (:,i,j,ip), leaf%ustar            (i,j,ip) &
                          , leaf%tstar            (i,j,ip), leaf%rstar            (i,j,ip) &
                          , leaf%cstar            (i,j,ip), leaf%veg_albedo       (i,j,ip) &
                          , leaf%veg_fracarea     (i,j,ip), leaf%veg_lai          (i,j,ip) &
                          , leaf%veg_tai          (i,j,ip), leaf%veg_rough        (i,j,ip) &
                          , leaf%veg_height       (i,j,ip), leaf%patch_area       (i,j,ip) &
                          , leaf%patch_rough      (i,j,ip), leaf%patch_wetind     (i,j,ip) &
                          , leaf%leaf_class       (i,j,ip), leaf%soil_rough       (i,j,ip) &
                          , leaf%sfcwater_nlev    (i,j,ip), leaf%stom_resist      (i,j,ip) &
                          , leaf%ground_rsat      (i,j,ip), leaf%ground_rvap      (i,j,ip) &
                          , leaf%ground_temp      (i,j,ip), leaf%ground_fliq      (i,j,ip) &
                          , leaf%veg_water        (i,j,ip), leaf%veg_hcap         (i,j,ip) &
                          , leaf%veg_energy       (i,j,ip), leaf%can_prss         (i,j,ip) &
                          , leaf%can_theiv        (i,j,ip), leaf%can_theta        (i,j,ip) &
                          , leaf%can_rvap         (i,j,ip), leaf%can_co2          (i,j,ip) &
                          , leaf%sensible         (i,j,ip), leaf%evap             (i,j,ip) &
                          , leaf%transp           (i,j,ip), leaf%gpp              (i,j,ip) &
                          , leaf%plresp           (i,j,ip), leaf%resphet          (i,j,ip) &
                          , leaf%veg_ndvip        (i,j,ip), leaf%veg_ndvic        (i,j,ip) &
                          , leaf%veg_ndvif        (i,j,ip), radiate%rshort        (i,j)    &
                          , radiate%cosz          (i,j)   , ip,i,j)
                     end if
                  end if

                  !------------------------------------------------------------------------!
                  !    Converting the surface water internal energy back to an intensive   !
                  ! variable.                                                              !
                  !------------------------------------------------------------------------!
                  ksn = nint(leaf%sfcwater_nlev(i,j,ip))
                  leaf%sfcwater_energy(:,i,j,ip) = 0.
                  do k=1,ksn
                     leaf%sfcwater_energy(k,i,j,ip) = sfcw_energy_ext(k)                   &
                                                    / leaf%sfcwater_mass(k,i,j,ip)
                  end do
                  !------------------------------------------------------------------------!
               end select
            end do tloop
         end do ploop1
      end do iloop1
   end do jloop1

   !---------------------------------------------------------------------------------------!
   !     Normalize accumulated fluxes and albedo seen by atmosphere over model timestep    !
   ! dtlt.                                                                                 !
   !---------------------------------------------------------------------------------------!
   do j = ja,jz
      do i = ia,iz
         turb%sflux_u(i,j) = turb%sflux_u(i,j) * dtll_factor * dens2(i,j)
         turb%sflux_v(i,j) = turb%sflux_v(i,j) * dtll_factor * dens2(i,j)
         turb%sflux_w(i,j) = turb%sflux_w(i,j) * dtll_factor * dens2(i,j)
         turb%sflux_t(i,j) = turb%sflux_t(i,j) * dtll_factor * dens2(i,j)
         turb%sflux_r(i,j) = turb%sflux_r(i,j) * dtll_factor * dens2(i,j)
         turb%sflux_c(i,j) = turb%sflux_c(i,j) * dtll_factor * dens2(i,j)
      end do
   end do

   if (ilwrtyp > 0 .or. iswrtyp > 0) then
      do j = ja,jz
         do i = ia,iz
            radiate%albedt (i,j) = radiate%albedt (i,j) * dtll_factor
            radiate%rlongup(i,j) = radiate%rlongup(i,j) * dtll_factor
         end do
      end do
   end if


   return
end subroutine leaf3
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine leaftw(mzg,mzs,np,soil_water, soil_energy,soil_text,sfcwater_mass               &
                 ,sfcwater_depth,ustar,tstar,rstar,cstar,veg_albedo,veg_fracarea,veg_lai   &
                 ,veg_tai,veg_rough,veg_height,patch_area,patch_rough,patch_wetind         &
                 ,leaf_class,soil_rough,sfcwater_nlev,stom_resist,ground_rsat,ground_rvap  &
                 ,ground_temp,ground_fliq,veg_water,veg_hcap,veg_energy,can_prss,can_theiv &
                 ,can_theta,can_rvap,can_co2,sensible,evap,transp,gpp,plresp,resphet       &
                 ,veg_ndvip,veg_ndvic,veg_ndvif,rshort,cosz,ip,i,j)

   use leaf_coms
   use mem_leaf
   use rconstants
   use mem_scratch
   use therm_lib   , only : qwtk        & ! subroutine
                          , qtk         & ! subroutine
                          , hpqz2temp   & ! function
                          , idealdenssh ! ! function
   use catt_start  , only : CATT        ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                , intent(in)    :: mzg
   integer                , intent(in)    :: mzs
   integer                , intent(in)    :: np
   real   , dimension(mzg), intent(inout) :: soil_water
   real   , dimension(mzg), intent(inout) :: soil_energy
   real   , dimension(mzg), intent(inout) :: soil_text
   real   , dimension(mzs), intent(inout) :: sfcwater_mass
   real   , dimension(mzs), intent(inout) :: sfcwater_depth
   real                   , intent(in)    :: ustar
   real                   , intent(in)    :: tstar
   real                   , intent(in)    :: rstar
   real                   , intent(in)    :: cstar
   real                   , intent(in)    :: veg_albedo
   real                   , intent(in)    :: veg_fracarea
   real                   , intent(in)    :: veg_lai
   real                   , intent(in)    :: veg_tai
   real                   , intent(in)    :: veg_rough
   real                   , intent(in)    :: veg_height
   real                   , intent(in)    :: patch_area
   real                   , intent(in)    :: patch_rough
   real                   , intent(in)    :: patch_wetind
   real                   , intent(in)    :: leaf_class
   real                   , intent(in)    :: soil_rough
   real                   , intent(inout) :: sfcwater_nlev
   real                   , intent(inout) :: stom_resist
   real                   , intent(inout) :: ground_rsat
   real                   , intent(inout) :: ground_rvap
   real                   , intent(inout) :: ground_temp
   real                   , intent(inout) :: ground_fliq
   real                   , intent(inout) :: veg_water
   real                   , intent(inout) :: veg_hcap
   real                   , intent(inout) :: veg_energy
   real                   , intent(inout) :: can_prss
   real                   , intent(inout) :: can_theiv
   real                   , intent(inout) :: can_theta
   real                   , intent(inout) :: can_rvap
   real                   , intent(inout) :: can_co2
   real                   , intent(inout) :: sensible
   real                   , intent(inout) :: evap
   real                   , intent(inout) :: transp
   real                   , intent(inout) :: gpp
   real                   , intent(inout) :: plresp
   real                   , intent(inout) :: resphet
   real                   , intent(in)    :: veg_ndvip
   real                   , intent(in)    :: veg_ndvic
   real                   , intent(in)    :: veg_ndvif
   real                   , intent(in)    :: rshort
   real                   , intent(in)    :: cosz
   !----- Local variables. ----------------------------------------------------------------!
   integer                                :: ip
   integer                                :: k
   integer                                :: nveg
   integer                                :: ksn
   integer                                :: nsoil
   integer                                :: ksnnew
   integer                                :: newlayers
   integer                                :: nlayers
   integer                                :: kold
   integer                                :: ktrans
   integer                                :: nsl
   integer                                :: kzs
   integer                                :: i
   integer                                :: j
   real                                   :: sfcw_energy_int
   real                                   :: stretch
   real                                   :: thik
   real                                   :: pf
   real                                   :: snden
   real                                   :: vegfracc
   real                                   :: qwfree
   real                                   :: wfree
   real                                   :: depthgain
   real                                   :: totsnow
   real                                   :: qw
   real                                   :: w
   real                                   :: qwt
   real                                   :: wt
   real                                   :: soilhcap
   real                                   :: fac
   real                                   :: wfreeb
   real                                   :: depthloss
   real                                   :: soilcap
   real                                   :: sndenmax
   real                                   :: sndenmin
   real                                   :: snowmin
   real                                   :: hold
   real                                   :: wtnew
   real                                   :: wtold
   real                                   :: wdiff
   real                                   :: watermid
   real                                   :: available_water
   real   , dimension(mzg)                :: available_layer
   real                                   :: extracted_water
   real                                   :: wloss
   real                                   :: qwloss
   real                                   :: soilcond
   real                                   :: waterfrac
   !----- Locally saved variables. --------------------------------------------------------!
   logical                , save          :: ncall=.true.
   real, dimension(20)    , save          :: thicknet
   real, dimension(20,20) , save          :: thick
   !---------------------------------------------------------------------------------------!


   !----- Initialise some levels. ---------------------------------------------------------!
   do k = 1,mzg
      dslz   (k) = slz(k+1) - slz(k)
      dslzi  (k) = 1. / dslz(k)
      dslzidt(k) = dslzi(k) * dtll
      slzt   (k) = .5 * (slz(k) + slz(k+1))
   end do
   do k = 2,mzg
      dslzt   (k) = slzt(k) - slzt(k-1)
      dslzti  (k) = 1. / dslzt(k)
      dslztidt(k) = dslzti(k) * dtll
   end do

   !----- Initialize snow thickness scaling array. ----------------------------------------!
   if (ncall) then
      ncall = .false.
      stretch = 2.0
      do kzs = 1,mzs
         thik = 1.0
         thicknet(kzs) = 0.0
         do k = 1,(kzs+1)/2
            thick(k,kzs) = thik
            thick(kzs+1-k,kzs) = thik
            thicknet(kzs) = thicknet(kzs) + 2. * thik
            thik = thik * stretch
         end do
         if ((kzs+1)/2 /= kzs/2) thicknet(kzs) = thicknet(kzs) - thik/stretch
         do k = 1,kzs
            thick(k,kzs) = thick(k,kzs) / thicknet(kzs)
         end do
      end do
   end if

   nveg = nint(leaf_class)
   ksn = nint(sfcwater_nlev)


   !---------------------------------------------------------------------------------------!
   !    Remove water from all layers up to the root depth for transpiration.  Bottom       !
   ! k-index in root zone is kroot(nveg).  Limit soil moisture to values above soilwp.     !
   ! The available water profile is given in kg/m2.                                        !
   !---------------------------------------------------------------------------------------!
   nsl    = nint(soil_text(mzg))
   ktrans = kroot(nveg)
   available_layer(:)   = 0.
   available_water      = 0.
   do k = mzg,ktrans,-1
      nsoil              = nint(soil_text(k))
      available_layer(k) = dslz(k) * wdns                                                  &
                         * max(0.0, fracliq(k) * (soil_water(k)-soilwp(nsoil)))
      available_water    = available_water + available_layer(k)
   end do


   !---------------------------------------------------------------------------------------!
   !      Evaluate any exchanges of heat and moisture to or from vegetation, apply moist-  !
   ! ure and heat changes to vegetation, and evaluate the resistance parameter rd between  !
   ! canopy air and the top soil or snow surface.                                          !
   !---------------------------------------------------------------------------------------!
   call leaf_canopy(mzg,mzs,ksn,soil_energy,soil_water,soil_text,sfcwater_mass,ustar       &
                   ,tstar,rstar,cstar,soil_rough,veg_rough,veg_height,veg_lai,veg_tai      &
                   ,veg_water,veg_hcap,veg_energy,leaf_class,veg_fracarea,stom_resist      &
                   ,can_prss,can_theiv,can_theta,can_rvap,can_co2,sensible,evap,transp,gpp &
                   ,plresp,resphet,ground_rsat,ground_rvap,ground_temp,ground_fliq         &
                   ,available_water,rshort,i,j,ip)

   !----- Compute soil and effective snow heat resistance times layer depth (rfactor). ----!
   do k = 1,mzg
      nsoil     = nint(soil_text(k))
      waterfrac = soil_water(k) / slmsts(nsoil)
      soilcond  = soilcond0(nsoil)                                                         &
                + waterfrac * (soilcond1(nsoil) + waterfrac * soilcond2(nsoil))
      rfactor(k) = dslz(k) / soilcond
   end do

   do k = 1,ksn
      snden = sfcwater_mass(k) / sfcwater_depth(k)
      rfactor(k+mzg) = sfcwater_depth(k) / (1.093e-3 * exp(.028 * tempk(k+mzg))            &
                    * (.030 + snden * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))
   end do

   !----- Find soil and snow internal sensible heat fluxes [W/m2]. ------------------------!
   hfluxgsc(1) = 0.
   do k = 2,mzg+ksn
      hfluxgsc(k) = - (tempk(k) - tempk(k-1)) / ((rfactor(k) + rfactor(k-1)) * .5)
   end do

   !---------------------------------------------------------------------------------------!
   !     Heat flux at soil or snow top from longwave, sensible, and upward latent heat     !
   ! fluxes [W/m^2].                                                                       !
   !---------------------------------------------------------------------------------------!
   hfluxgsc(ksn+1+mzg) = hflxgc + qwflxgc - rlonga_gs - rlongv_gs + rlonggs_v + rlonggs_a

   !---------------------------------------------------------------------------------------!
   !      Update soil internal energy values [J/m3] and snow/surface water internal energy !
   ! values [J/m2] from sensible heat, upward water vapor (latent heat), longwave, and     !
   ! shortwave fluxes.  This excludes effects of dew/frost formation, precipitation, shed- !
   ! ding, and percolation.  Update top soil or snow moisture from evaporation only.       !
   !---------------------------------------------------------------------------------------!
   do k = 1,mzg
      soil_energy(k) = soil_energy(k) + dslzidt(k) * (hfluxgsc(k) - hfluxgsc(k+1))
   end do
   soil_energy(mzg)  = soil_energy(mzg) + dslzidt(mzg) * rshort_g

   do k = 1,ksn
      sfcw_energy_ext(k) = sfcw_energy_ext(k)                                              &
                         + dtll * (hfluxgsc(k+mzg) - hfluxgsc(k+1+mzg) + rshort_s(k)) 
   end do

   if (ksn == 0) then
      soil_water(mzg) = soil_water(mzg) - wdnsi * wflxgc * dslzidt(mzg)
   else
      sfcwater_mass(ksn) = max(0.,sfcwater_mass(ksn) - wflxgc * dtll)
   endif

   !---------------------------------------------------------------------------------------!
   !      New moisture, qw, and depth from dew/frost formation, precipitation, shedding,   !
   ! and percolation.  ksnnew is the layer that receives the new condensate that comes     !
   ! directly from the air above.  If there is no pre-existing snowcover, this is a        !
   ! temporary snow/surface water layer.                                                   !
   !---------------------------------------------------------------------------------------!
   if (pcpgl + wshed_tot + dewgnd_tot > min_sfcwater_mass) then
      ksnnew = max(ksn,1)
      vegfracc = 1. - veg_fracarea
      qwfree    = qpcpgl * vegfracc + (qdewgnd_tot + qwshed_tot ) * veg_fracarea
      wfree     =  pcpgl * vegfracc + (dewgnd_tot  + wshed_tot  ) * veg_fracarea
      depthgain = dpcpgl * vegfracc + (dewgnd_tot  + wshed_tot  ) * veg_fracarea * wdnsi
   else
      ksnnew    = ksn
      qwfree    = 0.
      wfree     = 0.
      depthgain = 0.
   end if

   if (ksnnew > 0) then

      !------------------------------------------------------------------------------------!
      !     Transfer water downward through snow layers by percolation. Fracliq is the     !
      ! fraction of liquid in the snowcover or surface water. Wfree is the quantity of     !
      ! that liquid in kg/m2 which is free (not attached to snowcover) and therefore       !
      ! available to soak into the layer below). Soilcap is the capacity of the top soil   !
      ! layer in kg/m2 to accept surface water.  Wfree in the lowest snow layer is limited !
      ! by this value.                                                                     !
      !------------------------------------------------------------------------------------!
      totsnow = 0.
      nsoil = nint(soil_text(mzg))

      do k = ksnnew,1,-1
         qw = sfcw_energy_ext(k) + qwfree
         w  = sfcwater_mass(k)   + wfree

         !---------------------------------------------------------------------------------!
         !     If (only) snow layer is too thin for computational stability, bring it to   !
         ! thermal equilibrium with the top soil layer by exchanging heat between them.    !
         !---------------------------------------------------------------------------------!
         if (ksnnew == 1 .and. sfcwater_mass(k) < 3.) then
            qwt      = qw + soil_energy(mzg) * dslz(mzg)
            wt       = w  + soil_water(mzg) * wdns * dslz(mzg)
            soilhcap = slcpd(nsoil) * dslz(mzg)
            call qwtk(qwt,wt,soilhcap,tempk(k+mzg),fracliq(k+mzg))
            qw = w * ( fracliq(k+mzg)*cliq*(tempk(k+mzg)-tsupercool)                       &
                     + (1.-fracliq(k+mzg))*cice*tempk(k+mzg)  )
            tempk(mzg)       = tempk(k+mzg)
            fracliq(mzg)     = fracliq(k+mzg)
            soil_energy(mzg) = (qwt - qw) * dslzi(mzg)
         else
            call qtk(qw/w,tempk(k+mzg),fracliq(k+mzg))
         end if

         !---------------------------------------------------------------------------------!
         !    Shed liquid in excess of a 1:9 liquid-to-ice ratio.  Limit this shed amount  !
         ! (wfreeb) in lowest snow layer to amount top soil layer can hold.                !
         !---------------------------------------------------------------------------------!
         wfreeb    = max (0., w * (fracliq(k+mzg) - .1) / 0.9)
         depthloss = wfreeb * wdnsi
         if (k == 1) then
            soilcap = wdns * max (0.,-slz(mzg) * (slmsts(nsoil) - soil_water(mzg)))
            wfreeb  = min(wfreeb, soilcap)
            
            qwfree  = wfreeb * cliq * (tempk(k+mzg) - tsupercool)
            soil_water(mzg) = soil_water(mzg) + wdnsi * wfreeb * dslzi(mzg)
            soil_energy(mzg) = soil_energy(mzg) + qwfree * dslzi(mzg)
         else
            qwfree = wfreeb * cliq * (tempk(k+mzg) - tsupercool)
         end if

         sfcwater_mass(k)  = w - wfreeb
         sfcwater_depth(k) = sfcwater_depth(k) + depthgain - depthloss

         if (sfcwater_mass(k) < min_sfcwater_mass) then
            sfcw_energy_ext(k) = 0.
            sfcwater_mass(k)   = 0.
            sfcwater_depth(k)  = 0.
         else
            totsnow = totsnow + sfcwater_mass(k)
            sfcw_energy_ext(k) = qw - qwfree
         end if

         !----- Temporary simple evolution of snow layer depth and density. ---------------!
         sfcwater_depth(k) = sfcwater_depth(k) * (1. - dtll / 1.e5)
         snden             = sfcwater_mass(k)  / max(min_sfcwater_depth,sfcwater_depth(k))
         sndenmax          =  wdns
         sndenmin          = max( 30., 200. * (wfree + wfreeb)                             &
                                     / max(min_sfcwater_mass,sfcwater_mass(k)))
         snden             = min (sndenmax, max (sndenmin, snden))
         sfcwater_depth(k) = sfcwater_mass(k) / snden
         wfree             = wfreeb
         depthgain         = depthloss
      end do

      !----- Re-distribute snow layers to maintain prescribed distribution of mass. -------!
      if (totsnow < min_sfcwater_mass) then
         sfcwater_nlev      = 0.
         sfcwater_mass(1)   = 0.
         sfcw_energy_ext(1) = 0.
         sfcwater_depth(1)  = 0.
      else
         nlayers   = ksnnew
         snowmin   = 3.0
         newlayers = 1
         do k = 2,mzs
            if (snowmin * thicknet(k) <= totsnow .and.                                     &
                sfcw_energy_ext(k) < qliqt3*sfcwater_mass(k)) then
               newlayers = newlayers + 1
            end if
         end do
         newlayers     = min (newlayers, mzs, nlayers+1)
         sfcwater_nlev = float(newlayers)

         kold  = 1
         wtnew = 1.
         wtold = 1.
         do k = 1,newlayers
            vctr14(k) = totsnow * thick(k,newlayers)
            vctr16(k) = 0.
            vctr18(k) = 0.
            inner_loop: do
               wdiff = wtnew * vctr14(k) - wtold * sfcwater_mass(kold)
               if (wdiff > 0.) then
                  vctr16(k) = vctr16(k) + wtold * sfcw_energy_ext(kold)
                  vctr18(k) = vctr18(k) + wtold * sfcwater_depth(kold)
                  wtnew     = wtnew - wtold * sfcwater_mass(kold) / vctr14(k)
                  kold      = kold + 1
                  wtold     = 1.
                  if (kold <= ksn) cycle inner_loop
               else
                  vctr16(k) = vctr16(k) + wtnew * vctr14(k) * sfcw_energy_ext(kold)        &
                            / max(min_sfcwater_mass,sfcwater_mass(kold))
                  vctr18(k) = vctr18(k) + wtnew * vctr14(k) * sfcwater_depth(kold)         &
                            / max(min_sfcwater_mass,sfcwater_mass(kold))
                  wtold     = wtold - wtnew * vctr14(k) / sfcwater_mass(kold)
                  wtnew     = 1.
               end if
               exit inner_loop
            end do inner_loop
         end do

         do k = 1,newlayers
            sfcwater_mass(k)   = vctr14(k)
            sfcw_energy_ext(k) = vctr16(k) 
            sfcwater_depth(k)  = vctr18(k)
         end do
      end if
   end if

   !---------------------------------------------------------------------------------------!
   !     Compute gravitational potential plus moisture potential z + psi [m], liquid water !
   ! content [m], and half the remaining water capacity [m].                               !
   !---------------------------------------------------------------------------------------!
   do k = 1,mzg
      nsoil = nint(soil_text(k))
      psiplusz(k) = slzt(k) + slpots(nsoil) * (slmsts(nsoil)/soil_water(k))**slbs(nsoil)
      soil_liq(k) = dslz(k) * max(0.,(soil_water(k) - soilcp(nsoil))*fracliq(k))
      half_soilair(k) = (slmsts(nsoil) - soil_water(k)) * dslz(k) * .5
   end do

   !---------------------------------------------------------------------------------------!
   !     Find amount of water transferred between soil layers (wflux) [m] modulated by the !
   ! liquid water fraction.                                                                !
   !---------------------------------------------------------------------------------------!
   wflux(1) = 0.
   wflux(mzg+1) = 0.
   qwflux(1) = 0.
   qwflux(mzg+1) = 0.

   do k = 2,mzg
      nsoil    = nint(soil_text(k))
      watermid = 0.5 * (soil_water(k) + soil_water(k-1))
      wflux(k) = dslztidt(k) * slcons1(k,nsoil)                                            &
               * (watermid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)                     &
               * (psiplusz(k-1) - psiplusz(k)) * .5 * (fracliq(k) + fracliq(k-1))

      !------------------------------------------------------------------------------------!
      !     Limit water transfers to prevent over-saturation and over-depletion. Compute q !
      ! transfers between soil layers (qwflux) [J/m2].                                     !
      !------------------------------------------------------------------------------------!
      if (wflux(k) > 0.) then
         wflux(k) = min(wflux(k),soil_liq(k-1),half_soilair(k))
      else
         wflux(k) = - min(-wflux(k),soil_liq(k),half_soilair(k-1))
      endif
      qwflux(k) = wflux(k) * cliqvlme * (tempk(k) - tsupercool)

   end do

   !----- Update soil moisture (impose minimum value of soilcp) and q value. --------------!
   do k = 1,mzg
      nsoil = nint(soil_text(k))
      soil_water(k)  = max(soilcp(nsoil),soil_water(k) - dslzi(k) * (wflux(k+1)-wflux(k)))
      soil_energy(k) = soil_energy(k) - dslzi(k) * (qwflux(k+1)-qwflux(k))
   enddo

   !---------------------------------------------------------------------------------------!
   !     Update soil moisture by extracting the water lost due to transpiration.           !
   !---------------------------------------------------------------------------------------!
   extracted_water = transp_tot * dtll
   if (extracted_water > 0. .and. available_water > 0.) then
      do k = ktrans,mzg
         wloss  = wdnsi * dslzi(k) * extracted_water * available_layer(k) / available_water
         qwloss = wloss * cliqvlme * (tempk(k) - tsupercool)
         
         soil_water(k)  = soil_water(k)  - wloss
         soil_energy(k) = soil_energy(k) - qwloss
      end do
   end if

   !---------------------------------------------------------------------------------------!
   ! Compute ground vap mxrat for availability on next timestep; put into ground_rsat.     !
   !---------------------------------------------------------------------------------------!
   newlayers = nint(sfcwater_nlev)
   if (newlayers == 0) then
      sfcw_energy_int = 0.
   else
      if (sfcwater_mass(newlayers) > min_sfcwater_mass) then
         sfcw_energy_int = sfcw_energy_ext(newlayers) / sfcwater_mass(newlayers)
      else 
         sfcw_energy_int = 0.
      end if
   end if
   call leaf_grndvap(soil_energy(mzg),soil_water(mzg),soil_text(mzg),sfcw_energy_int       &
                    ,sfcwater_nlev,can_rvap,can_prss,ground_rsat,ground_rvap,ground_temp   &
                    ,ground_fliq)

   return
end subroutine leaftw
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the variables that depend on heat, water, and (eventually)  !
! carbon exchanges with the canopy air space and vegetation.                               !
!------------------------------------------------------------------------------------------!
subroutine leaf_canopy(mzg,mzs,ksn,soil_energy,soil_water,soil_text,sfcwater_mass          &
                      ,ustar,tstar,rstar,cstar,soil_rough,veg_rough,veg_height,veg_lai     &
                      ,veg_tai,veg_water,veg_hcap,veg_energy,leaf_class,veg_fracarea       &
                      ,stom_resist,can_prss,can_theiv,can_theta,can_rvap,can_co2,sensible  &
                      ,evap,transp,gpp,plresp,resphet,ground_rsat,ground_rvap,ground_temp  &
                      ,ground_fliq,available_water,rshort,i,j,ip)
   use leaf_coms
   use rconstants
   use therm_lib , only : eslif          & ! function
                        , rslif          & ! function
                        , thetaeiv2thil  & ! function
                        , idealdenssh    & ! function
                        , tslif          & ! function
                        , qwtk           ! ! subroutine
   use catt_start, only : CATT
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer             , intent(in)    :: mzg
   integer             , intent(in)    :: mzs
   integer             , intent(in)    :: ksn
   real, dimension(mzg), intent(in)    :: soil_energy
   real, dimension(mzg), intent(in)    :: soil_water
   real, dimension(mzg), intent(in)    :: soil_text
   real, dimension(mzs), intent(in)    :: sfcwater_mass
   real                , intent(in)    :: ustar
   real                , intent(in)    :: tstar
   real                , intent(in)    :: rstar
   real                , intent(in)    :: cstar
   real                , intent(in)    :: soil_rough
   real                , intent(in)    :: veg_rough
   real                , intent(in)    :: veg_height
   real                , intent(in)    :: veg_lai
   real                , intent(in)    :: veg_tai
   real                , intent(in)    :: leaf_class
   real                , intent(in)    :: veg_fracarea
   real                , intent(in)    :: ground_rsat
   real                , intent(in)    :: ground_rvap
   real                , intent(in)    :: ground_temp
   real                , intent(in)    :: ground_fliq
   real                , intent(in)    :: available_water
   real                , intent(in)    :: rshort
   real                , intent(inout) :: veg_water
   real                , intent(inout) :: veg_hcap
   real                , intent(inout) :: veg_energy
   real                , intent(inout) :: stom_resist
   real                , intent(inout) :: can_prss
   real                , intent(inout) :: can_theiv
   real                , intent(inout) :: can_theta
   real                , intent(inout) :: can_rvap
   real                , intent(inout) :: can_co2
   real                , intent(inout) :: sensible
   real                , intent(inout) :: evap
   real                , intent(inout) :: transp
   real                , intent(inout) :: gpp
   real                , intent(inout) :: plresp
   real                , intent(inout) :: resphet
   !----- Local arguments. ----------------------------------------------------------------!
   integer                             :: k
   integer                             :: kk
   integer                             :: nveg
   integer                             :: nsoil
   integer                             :: i
   integer                             :: j
   integer                             :: ip
   real                                :: aux
   real                                :: c2
   real                                :: c3
   real                                :: c4
   real                                :: dsm
   real                                :: es
   real                                :: fac
   real                                :: factv
   real                                :: fthi
   real                                :: ftlo
   real                                :: frad
   real                                :: fswp
   real                                :: fvpd
   real                                :: qwtot
   real                                :: rasgnd
   real                                :: rasveg
   real                                :: rleaf
   real                                :: rsatveg
   real                                :: sigmaw 
   real                                :: slai
   real                                :: stai
   real                                :: slpotv
   real                                :: swp
   real                                :: vpd
   real                                :: wtemp
   real                                :: wtroot
   real                                :: x
   real                                :: zognd
   real                                :: zoveg
   real                                :: zdisp
   real                                :: zveg
   real                                :: wflx
   real                                :: dewgndflx
   real                                :: qdewgndflx
   real                                :: transp_test
   real                                :: transp_wilt
   real                                :: transp_wilti
   real                                :: rc
   real                                :: rc_inf
   real                                :: wshed_loc
   real                                :: qwshed_loc
   real                                :: old_veg_water
   real                                :: old_veg_energy
   real                                :: old_can_rvap
   real                                :: old_can_theiv
   real                                :: old_can_shv
   real                                :: old_can_theta
   real                                :: old_can_temp
   real                                :: old_can_rhos
   !----- Local constants. ----------------------------------------------------------------!
   character(len=9)      , parameter   :: fmti='(a,1x,i6)'
   character(len=13)     , parameter   :: fmtf='(a,1x,es12.5)'
   character(len=3)      , parameter   :: fmtc='(a)'
   character(len=9)      , parameter   :: fmtl='(a,1x,l1)'
   !---------------------------------------------------------------------------------------!



   !----- Find the time step auxiliary variables. -----------------------------------------!
   dtllowcc = dtll / (can_depth * can_rhos)
   dtllohcc = dtll / (can_depth * can_rhos * cp * can_temp)
   dtlloccc = dtllowcc * mmdry
   !---------------------------------------------------------------------------------------!


   !----- Save the previous values for debugging. -----------------------------------------!
   old_can_theiv    = can_theiv
   old_can_rvap     = can_rvap
   old_can_shv      = can_shv
   old_can_theta    = can_theta
   old_can_temp     = can_temp
   old_can_rhos     = can_rhos



   !---------------------------------------------------------------------------------------!
   !     Compute ground-canopy resistance rd.  Assume zognd not affected by snow.  Assume  !
   ! (zoveg,zdisp) decrease linearly with snow depth, attaining the values (zognd,0) when  !
   ! vegetation is fully buried in snow.                                                   !
   !---------------------------------------------------------------------------------------!

   zognd = soil_rough
   zoveg = veg_rough * (1.-snowfac) + zognd * snowfac
   zdisp = veg_height * (1.-snowfac)
   zveg = zdisp / 0.63

   !---------------------------------------------------------------------------------------!
   ! The following value of ground-canopy resistance for the nonvegetated (bare soil or    !
   ! water) surface is from John Garratt.  It is 5/ustar and replaces the one from old     !
   ! LEAF.                                                                                 !
   !---------------------------------------------------------------------------------------!
   rasgnd = 5. / ustar

   if (veg_tai >= .1 .and. snowfac < .9) then

      !------------------------------------------------------------------------------------!
      !    If vegetation is sufficiently abundant and not covered by snow, compute  heat   !
      ! and moisture fluxes from vegetation to canopy, and flux resistance from soil or    !
      ! snow to canopy.                                                                    !
      !------------------------------------------------------------------------------------!
      factv       = log(geoht / zoveg) / (vonk * vonk * vels)
      aux         = exp(exar * (1. - (zdisp + zoveg) / zveg))
      rasveg      = factv * zveg / (exar * (zveg - zdisp)) * (exp(exar) - aux)
      c2          = max(0.,min(1., 1.1 * veg_tai / covr))
      rd          = rasgnd * (1. - c2) + rasveg * c2
      wshed_tot   = 0.
      qwshed_tot  = 0.
      transp_tot  = 0.
   else
      !------------------------------------------------------------------------------------!
      !     If the TAI is very small or if snow mostly covers the vegetation, bypass       !
      ! vegetation computations.  Set heat and moisture flux resistance rd between the     !
      ! "canopy" and snow or soil surface to its bare soil value.  Set shed precipitation  !
      ! heat and moisture to unintercepted values.                                         !
      !------------------------------------------------------------------------------------!
      wshed_tot  = pcpgl
      qwshed_tot = qpcpgl
      transp_tot = 0.
      rd         = rasgnd
   end if

   !---------------------------------------------------------------------------------------!
   !     Compute sensible heat and moisture fluxes between top soil or snow surface and    !
   ! canopy air.  wflxgc [kg/m2/s] is the upward vapor flux from soil or snow evaporation  !
   ! and dewgnd is the mass of dew that forms on the snow/soil surface this timestep; both !
   ! are defined as always positive or zero.                                               !
   !---------------------------------------------------------------------------------------!
   rdi = can_rhos / rd

   hflxgc     = cp * rdi * (ground_temp - can_temp)
   wflx       =      rdi * (ground_rsat - can_rvap)
   
   if (wflx >= 0.) then
      dewgndflx = 0.
   else
      dewgndflx = min(-wflx,(can_rvap - toodry) / dtllowcc)
   end if
   dewgnd_tot = dewgndflx * dtll

   if (ksn == 0) then
      wflxgc = max(0.,rdi * (ground_rvap - can_rvap))
   else
      wflxgc = max(0.,min(wflx,sfcwater_mass(ksn)/dtll))
   end if

   qdewgnd_tot = (alvi - alli * ground_fliq) * dewgnd_tot
   qdewgndflx  = (alvi - alli * ground_fliq) * dewgndflx
   qwflxgc     = (alvi - alli * ground_fliq) * wflxgc

   !---------------------------------------------------------------------------------------!
   !     If vegetation is sufficiently abundant and not covered by snow, compute heat and  !
   ! moisture fluxes from vegetation to canopy, and flux resistance from soil or snow to   !
   ! canopy.                                                                               !
   !---------------------------------------------------------------------------------------!
   if (veg_tai >= .1 .and. snowfac < .9) then


      !------------------------------------------------------------------------------------!
      !     Compute veg-canopy resistance rb.  Note that rb and rc are both defined        !
      ! WITHOUT the LAI factor; the LAI factor is included later in computing  fluxes that !
      ! involve rb and/or rc.                                                              !
      !------------------------------------------------------------------------------------!
      c4 = .01 * sqrt(ustar * c1)
      rb  = (1. + .5 * veg_tai) / c4

      !----- Soil water potential factor for stomatal control. ----------------------------!
      swp = -200.0
      nveg = nint(leaf_class)

      do k = kroot(nveg),mzg
         nsoil  = nint(soil_text(k))
         slpotv = slpots(nsoil) * (slmsts(nsoil) / soil_water(k)) ** slbs(nsoil)
         if (slpotv > swp) swp = slpotv
      end do
      swp = swp * grav * wdns

      !------------------------------------------------------------------------------------!
      !      Begin canopy air computations.                                                !
      !------------------------------------------------------------------------------------!

      !----- Calculate the saturation vapor pressure at TVEG. -----------------------------!
      es      = eslif(veg_temp)
      rsatveg = rslif(can_prss,veg_temp)
      
      !----- Compute mixing ratio at leaf surface using previous rc -----------------------!
      rc    = stom_resist
      rleaf = (rb * rsatveg + rc * can_rvap) / (rb + rc)
      vpd   = max((es - rleaf * can_prss / (ep + rleaf)),0.)

      !----- Evaluate 5 environmental factors and new rc ----------------------------------!
      ftlo = 1. + exp(-stlo * (veg_temp - btlo))
      fthi = 1. + exp(-sthi * (veg_temp - bthi))
      frad = 1. + exp(-srad * (rshort   - brad))
      fswp = 1. + exp(-sswp * (swp      - bswp))
      fvpd = 1. + exp(-svpd * (vpd      - bvpd))

      !----- 15-minute response time for stomatal conductance (must have dtll <= 900.). ---!
      rc_inf = ftlo * fthi * frad * fvpd * fswp * rcmin(nveg)
      rc     = 1. / (1. / rc + .0011 * dtll * (1. / rc_inf - 1. / rc))

      !------------------------------------------------------------------------------------!
      !    Transpiration can only happen if there is water.                                !
      !------------------------------------------------------------------------------------!
      if (available_water > 0.) then
         !---------------------------------------------------------------------------------!
         !    Limit maximum transpiration to be <= transp_max and less than the maximum    !
         ! water available by increasing rc if necessary.                                  !
         !---------------------------------------------------------------------------------!
         transp_test  = can_rhos * veg_lai * (rsatveg - can_rvap) / (rb + rc)
         transp_wilt  = min(available_water / dtll, transp_max)
         transp_wilti = 1./ transp_wilt
         if (transp_test > transp_wilt) then
            rc          = (rb + rc) * transp_test * transp_wilti - rb
         end if
      else
         transp_test = 0.
      end if
      stom_resist = rc

      !----- Flux of heat and moisture from vegetation to canopy. -------------------------!
      stai = veg_tai * (1. - snowfac)
      slai = veg_lai * (1. - snowfac)

      !----- Sensible heat flux from leaf to canopy. --------------------------------------!
      hflxvc = 2. * stai * cp * can_rhos * (veg_temp - can_temp) / rb

      c3     = can_rhos * (rsatveg - can_rvap)

      if (c3 >= 0.) then
         !----- Flow will go towards the canopy, allow evaporation and transpiration. -----!
         sigmaw     = min(1.,(veg_water / (.2 * stai)) ** .66667)
         wflxvc     = min(c3 * 2. * stai * sigmaw / rb,veg_water/dtll)
         if (transp_test > 0.) then
            transp_loc = max(0.,c3 * slai * (1.-sigmaw) / (rb + rc))
         else
            transp_loc = 0.
         end if
      else
         !----- Flow is towards the leaf, dew/frost formation and no transpiration. -------!
         wflxvc     = max(c3 * 2. * stai / rb,min(0.,(toodry-can_rvap)/ dtllowcc))
         transp_loc = 0.
      end if
      transp_tot = transp_tot + transp_loc

      !------------------------------------------------------------------------------------!
      !     Update vegetation moisture from precipitation, vapor flux with canopy,         !
      ! and shedding of excess moisture.                                                   !
      !------------------------------------------------------------------------------------!
      old_veg_water = veg_water
      wtemp         = veg_water - wflxvc * dtll + pcpgl * vf
      !----- If this result will lead to negative solution, then reduce evaporation. ------!
      if (wtemp < 0.0005 * stai) then
         wshed_loc = 0.
         wflxvc    = (veg_water + pcpgl * vf) / dtll
         veg_water = 0.
      else
         wshed_loc = max(0.,wtemp - 0.2 * stai)
         wshed_tot = wshed_tot + wshed_loc
         veg_water = wtemp - wshed_loc
      end if
      !------------------------------------------------------------------------------------!

      !------ Find the associated latent heat flux from vegetation to canopy. -------------!
      qwflxvc     = wflxvc     * (alvi - alli * veg_fliq)
      qtransp_loc = transp_loc *  alvl !----- Liquid phase only in transpiration. ---------!
      !------------------------------------------------------------------------------------!


      !----- Exchange heat between vegetation and precipitation in implicit scheme --------!
      qwshed_loc = wshed_loc * (     veg_fliq  * cliq * (veg_temp - tsupercool)            &
                               + (1.-veg_fliq) * cice *  veg_temp )
      qwshed_tot = qwshed_tot + qwshed_loc

      !------------------------------------------------------------------------------------!
      !      Update vegetation internal energy from radiative fluxes and sensible and      !
      ! latent heat transfer with canopy air.                                              !
      !------------------------------------------------------------------------------------!
      old_veg_energy = veg_energy
      veg_energy     = veg_energy                                                          &
                     + dtll * ( rshort_v + rlonga_v + rlonggs_v - rlongv_gs - rlongv_a     &
                              - hflxvc   - qwflxvc   - qtransp_loc)                        &
                     + qpcpgl * vf - qwshed_loc

      !------------------------------------------------------------------------------------!
      !     Find the atmosphere -> canopy fluxes.                                          !
      !------------------------------------------------------------------------------------!
      rho_ustar  = can_rhos  * ustar
      eflxac     = rho_ustar * estar * cp * can_temp ! Enthalpy exchange
      wflxac     = rho_ustar * rstar                 ! Water vapour exchange
      cflxac     = rho_ustar * cstar * mmdryi        ! CO2 exchange
      !------------------------------------------------------------------------------------!




      !------------------------------------------------------------------------------------!
      !      Update enthalpy, CO2, and canopy mixing ratio.                                !
      !------------------------------------------------------------------------------------!
      can_lntheiv      = can_lntheiv                                                       &
                       + dtllohcc * ( hflxgc + hflxvc + eflxac                             &
                                    + qwflxgc - qdewgndflx + qwflxvc + qtransp_loc)
      can_rvap         = can_rvap                                                          &
                       + dtllowcc * (wflxgc - dewgndflx + wflxvc + transp_loc + wflxac)
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !      CO2 exchange with biosphere is currently not computed in LEAF-3... Feel       !
      ! free to include here: cflxgc would be the root respiration, and cflxvc R-GPP.      !
      !------------------------------------------------------------------------------------!
      cflxgc   = 0.
      cflxvc   = 0.
      cflxcv   = 0.
      can_co2  = can_co2  + dtllowcc * mmdry * (cflxgc + cflxvc + cflxcv + cflxac)
      !------------------------------------------------------------------------------------!



      !----- Update the fluxes. -----------------------------------------------------------!
      sensible = sensible + (hflxvc  + hflxgc ) * dtll_factor
      evap     = evap     + (qwflxvc + qwflxgc) * dtll_factor
      transp   = transp   + (qtransp_loc      ) * dtll_factor

      !----- Find the vegetation diagnostic variables for next time step. -----------------!
      call qwtk(veg_energy,veg_water,veg_hcap,veg_temp,veg_fliq)

      !----- Find the canopy air diagnostic variables for next time step. -----------------!
      can_theiv = exp(can_lntheiv)
      can_shv   = can_rvap / (can_rvap + 1.)
      can_theta = thetaeiv2thil(can_theiv,can_prss,can_rvap)
      can_temp  = cpi * can_theta * can_exner
      can_rhos  = idealdenssh(can_prss,can_temp,can_shv)
      !------------------------------------------------------------------------------------!


      if (can_lntheiv /= can_lntheiv .or. can_theiv /= can_theiv) then
         nsoil = nint(soil_text(mzg))
         write (unit=*,fmt=fmtc) '------------ Canopy theta_Eiv is screwed. ------------'
         write (unit=*,fmt=fmtc) ' - Vegetated patch. '
         write (unit=*,fmt=fmti) ' - SOIL_TEXT        = ',nsoil
         write (unit=*,fmt=fmti) ' - LEAF_CLASS       = ',nveg
         write (unit=*,fmt=fmti) ' - KSN              = ',ksn
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - CAN_THEIV        = ',can_theiv
         write (unit=*,fmt=fmtf) ' - OLD_CAN_THEIV    = ',old_can_theiv
         write (unit=*,fmt=fmtf) ' - CAN_RVAP         = ',can_rvap
         write (unit=*,fmt=fmtf) ' - OLD_CAN_RVAP     = ',old_can_rvap
         write (unit=*,fmt=fmtf) ' - CAN_CO2          = ',can_co2
         write (unit=*,fmt=fmtf) ' - CAN_PRSS         = ',can_prss
         write (unit=*,fmt=fmtf) ' - OLD_CAN_THETA    = ',old_can_theta
         write (unit=*,fmt=fmtf) ' - OLD_CAN_TEMP     = ',old_can_temp
         write (unit=*,fmt=fmtf) ' - OLD_CAN_RHOS     = ',old_can_rhos
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - VEG_ENERGY       = ',veg_energy
         write (unit=*,fmt=fmtf) ' - OLD_VEG_ENERGY   = ',old_veg_energy
         write (unit=*,fmt=fmtf) ' - VEG_WATER        = ',veg_water
         write (unit=*,fmt=fmtf) ' - OLD_VEG_WATER    = ',old_veg_water
         write (unit=*,fmt=fmtf) ' - VEG_HCAP         = ',veg_hcap
         write (unit=*,fmt=fmtf) ' - VEG_TEMP         = ',veg_temp
         write (unit=*,fmt=fmtf) ' - VEG_FLIQ         = ',veg_fliq
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - GROUND_RVAP      = ',ground_rvap
         write (unit=*,fmt=fmtf) ' - GROUND_RSAT      = ',ground_rsat
         write (unit=*,fmt=fmtf) ' - GROUND_TEMP      = ',ground_temp
         write (unit=*,fmt=fmtf) ' - GROUND_FLIQ      = ',ground_fliq
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - HFLXGC           = ',hflxgc
         write (unit=*,fmt=fmtf) ' - HFLXVC           = ',hflxvc
         write (unit=*,fmt=fmtf) ' - EFLXAC           = ',eflxac
         write (unit=*,fmt=fmtf) ' - QWFLXGC          = ',qwflxgc
         write (unit=*,fmt=fmtf) ' - QDEWGNDFLX       = ',qdewgndflx
         write (unit=*,fmt=fmtf) ' - QWFLXVC          = ',qwflxvc
         write (unit=*,fmt=fmtf) ' - QTRANSP_LOC      = ',qtransp_loc
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - RSHORT_V         = ',rshort_v
         write (unit=*,fmt=fmtf) ' - RLONGA_V         = ',rlonga_v
         write (unit=*,fmt=fmtf) ' - RLONGGS_V        = ',rlonggs_v
         write (unit=*,fmt=fmtf) ' - RLONGV_GS        = ',rlongv_gs
         write (unit=*,fmt=fmtf) ' - RLONGV_A         = ',rlongv_a
         write (unit=*,fmt=fmtf) ' - QPCPGL           = ',qpcpgl
         write (unit=*,fmt=fmtf) ' - VF               = ',vf
         write (unit=*,fmt=fmtf) ' - QWSHED_LOC       = ',qwshed_loc
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - WFLXGC           = ',wflxgc
         write (unit=*,fmt=fmtf) ' - DEWGNDFLX        = ',dewgndflx
         write (unit=*,fmt=fmtf) ' - TRANSP_LOC       = ',transp_loc
         write (unit=*,fmt=fmtf) ' - WFLXAC           = ',wflxac
         write (unit=*,fmt=fmtf) ' - WFLXVC           = ',wflxvc
         write (unit=*,fmt=fmtf) ' - PCPGL            = ',pcpgl
         write (unit=*,fmt=fmtf) ' - WSHED_LOC        = ',wshed_loc
         write (unit=*,fmt=fmtf) ' - VF               = ',vf
         write (unit=*,fmt=fmtc) '-----------------------------------------------------'
      end if

   else
      !------------------------------------------------------------------------------------!
      !    No vegetation, or vegetation buried in snow...                                  !
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the atmosphere -> canopy fluxes.                                          !
      !------------------------------------------------------------------------------------!
      rho_ustar  = can_rhos  * ustar
      eflxac     = rho_ustar * estar * cp * can_temp ! Enthalpy exchange
      wflxac     = rho_ustar * rstar                 ! Water vapour exchange
      cflxac     = rho_ustar * cstar * mmdryi        ! CO2 exchange
      !------------------------------------------------------------------------------------!


      !----- Update the canopy prognostic variables. --------------------------------------!
      can_lntheiv  = can_lntheiv  + dtllohcc * (hflxgc + qwflxgc - qdewgndflx + eflxac)
      can_rvap     = can_rvap     + dtllowcc * (wflxgc - dewgndflx + wflxac)

      !------------------------------------------------------------------------------------!
      !      CO2 exchange with biosphere is currently not computed in LEAF-3... Feel       !
      ! free to include here: cflxgc would be the root respiration, and cflxvc R-GPP.      !
      !------------------------------------------------------------------------------------!
      cflxgc       = 0.
      can_co2      = can_co2      + dtlloccc * (cflxgc + cflxac)
      !------------------------------------------------------------------------------------!

      !----- Update the fluxes. -----------------------------------------------------------!
      sensible = sensible + hflxgc  * dtll_factor
      evap     = evap     + qwflxgc * dtll_factor

      !----- Find the diagnostic canopy air variables for next time step. -----------------!
      can_theiv = exp(can_lntheiv)
      can_shv   = can_rvap / (can_rvap + 1.)
      can_theta = thetaeiv2thil(can_theiv,can_prss,can_rvap)
      can_temp  = cpi * can_theta * can_exner
      can_rhos  = idealdenssh(can_prss,can_temp,can_shv)
      !------------------------------------------------------------------------------------!

      if (can_lntheiv /= can_lntheiv .or. can_theiv /= can_theiv) then
         write (unit=*,fmt=fmtc) '------------ Canopy theta_Eiv is screwed. ------------'
         write (unit=*,fmt=fmtc) ' - Non-vegetated patch. '
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - CAN_THEIV        = ',can_theiv
         write (unit=*,fmt=fmtf) ' - OLD_CAN_THEIV    = ',old_can_theiv
         write (unit=*,fmt=fmtf) ' - CAN_RVAP         = ',can_rvap
         write (unit=*,fmt=fmtf) ' - OLD_CAN_RVAP     = ',old_can_rvap
         write (unit=*,fmt=fmtf) ' - CAN_CO2          = ',can_co2
         write (unit=*,fmt=fmtf) ' - CAN_PRSS         = ',can_prss
         write (unit=*,fmt=fmtf) ' - OLD_CAN_THETA    = ',old_can_theta
         write (unit=*,fmt=fmtf) ' - OLD_CAN_TEMP     = ',old_can_temp
         write (unit=*,fmt=fmtf) ' - OLD_CAN_RHOS     = ',old_can_rhos
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - HFLXGC           = ',hflxgc
         write (unit=*,fmt=fmtf) ' - EFLXAC           = ',eflxac
         write (unit=*,fmt=fmtf) ' - QWFLXGC          = ',qwflxgc
         write (unit=*,fmt=fmtf) ' - QDEWGNDFLX       = ',qdewgndflx
         write (unit=*,fmt=fmtc) ' '
         write (unit=*,fmt=fmtf) ' - WFLXGC           = ',wflxgc
         write (unit=*,fmt=fmtf) ' - DEWGNDFLX        = ',dewgndflx
         write (unit=*,fmt=fmtf) ' - WFLXAC           = ',wflxac
         write (unit=*,fmt=fmtc) '-----------------------------------------------------'
      end if
   end if

   return
end subroutine leaf_canopy
!==========================================================================================!
!==========================================================================================!
