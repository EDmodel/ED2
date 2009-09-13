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
   real, dimension(mxp,myp)             :: l_ths2,l_rvs2,l_pis2,l_dens2,l_ups2
   real, dimension(mxp,myp)             :: l_vps2,l_zts2,l_co2s2
   !---------------------------------------------------------------------------------------!
   
   if (nstbot == 0) return

   ng=ngrid


   !----- Calling LEAF-3 main driver ------------------------------------------------------!
   call leaf3(mzp,mxp,myp,nzg,nzs,npatch,ia,iz,ja,jz,leaf_g (ng),basic_g (ng),turb_g (ng)  &
             ,radiate_g(ng),grid_g (ng),cuparm_g(ng),micro_g(ng),l_ths2,l_rvs2,l_co2s2     &
             ,l_pis2,l_dens2,l_ups2,l_vps2,l_zts2,teb_g(ng),tebc_g(ng))

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
                  , leaf_g(ng)%cstar               , leaf_g(ng)%veg_albedo                 &
                  , leaf_g(ng)%veg_fracarea        , leaf_g(ng)%veg_lai                    &
                  , leaf_g(ng)%veg_tai             , leaf_g(ng)%veg_rough                  &
                  , leaf_g(ng)%veg_height          , leaf_g(ng)%patch_area                 &
                  , leaf_g(ng)%patch_rough         , leaf_g(ng)%patch_wetind               &
                  , leaf_g(ng)%leaf_class          , leaf_g(ng)%soil_rough                 &
                  , leaf_g(ng)%sfcwater_nlev       , leaf_g(ng)%stom_resist                &
                  , leaf_g(ng)%ground_rsat         , leaf_g(ng)%ground_rvap                &
                  , leaf_g(ng)%veg_water           , leaf_g(ng)%veg_hcap                   &
                  , leaf_g(ng)%veg_energy          , leaf_g(ng)%can_prss                   &
                  , leaf_g(ng)%can_theta           , leaf_g(ng)%can_rvap                   &
                  , leaf_g(ng)%can_co2             , leaf_g(ng)%sensible                   &
                  , leaf_g(ng)%evap                , leaf_g(ng)%transp                     &
                  , leaf_g(ng)%gpp                 , leaf_g(ng)%plresp                     &
                  , leaf_g(ng)%resphet             , leaf_g(ng)%veg_ndvip                  &
                  , leaf_g(ng)%veg_ndvic           , leaf_g(ng)%veg_ndvif                  )
   
   return
end subroutine sfclyr
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This is the LEAF-3 main driver.                                                       !
!------------------------------------------------------------------------------------------!
subroutine leaf3(m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz,leaf,basic,turb,radiate,grid,cuparm,micro &
                ,ths2,rvs2,co2s2,pis2,dens2,ups2,vps2,zts2,teb,tebc)

   use mem_all
   use leaf_coms
   use rconstants
   use node_mod       , only : mynum         ! ! intent(in)
   use catt_start     , only : CATT          ! ! intent(in)
   use teb_spm_start  , only : TEB_SPM       ! ! intent(in)
   use mem_teb        , only : teb_vars      ! ! type
   use mem_teb_common , only : teb_common    ! ! type
   use therm_lib      , only : rslif         & ! function
                             , reducedpress  & ! function
                             , ptqz2enthalpy & ! function
                             , hpqz2temp     & ! function
                             , idealdenssh   & ! function
                             , qwtk          & ! subroutine
                             , level         ! ! intent(in)
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
   real, dimension(m2,m3), intent(inout) :: ths2,rvs2,co2s2,pis2,dens2,ups2,vps2,zts2
   !----- Local variables -----------------------------------------------------------------!
   type(teb_vars)                      :: teb
   type(teb_common)                    :: tebc
   integer                             :: i,j,k,ksn,ip,iter_leaf
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
   !    Define LEAF-3 and canopy time-split timesteps here.  This ensures that LEAF-3 will !
   ! not use a timestep longer than about 40 seconds, and canopy will not use a timestep   !
   ! longer than about 15 seconds.  This allows values of can_depth = 20. and              !
   ! hcapveg = 3.e4 as are now defined below.                                              !
   !---------------------------------------------------------------------------------------!
   niter_leaf  = max(1,nint(dtlt/40.+.4))
   niter_can   = max(1,nint(dtll/15.+.4))
   dtll_factor = 1. / float(niter_leaf)
   dtll        = dtlt * dtll_factor
   dtlc_factor = 1. / float(niter_can)
   dtlc        = dtll * dtlc_factor

   !----- Check whether we have CO2, and copy to an scratch array. ------------------------!
   if (co2_on) then
      call atob(m1*m2*m3,basic%co2p,scratch%vt3do)
   else
      call ae0(m1*m2*m3,scratch%vt3do,co2con(1))
   end if

   !----- Copy surface atmospheric variables into 2d arrays for input to LEAF. ------------!
   select case (if_adap)
   case (0)
      call sfc_fields(m1,m2,m3,ia,iz,ja,jz,jdim           , basic%theta   , basic%rv       &
                          , scratch%vt3do , basic%up      , basic%vp      , basic%dn0      &
                          ,  basic%pp     , basic%pi0     , grid%rtgt                      &
                          , zt,zm,ths2,rvs2,co2s2,ups2,vps2,pis2,dens2,zts2)
   case (1)
      call sfc_fields_adap(m1,m2,m3,ia,iz,ja,jz,jdim      , grid%flpu     , grid%flpv      &
                          , grid%flpw     , grid%topma    , grid%aru      , grid%arv       &
                          , basic%theta   , basic%rv      , scratch%vt3do , basic%up       &
                          , basic%vp      , basic%dn0     , basic%pp      , basic%pi0      &
                          , zt,zm,dzt,ths2,rvs2,co2s2,ups2,vps2,pis2,dens2,zts2)
   end select

   !----- Big domain loop -----------------------------------------------------------------!
   jloop1: do j = ja,jz
      iloop1: do i = ia,iz

         !----- Copy surface variables to single-column values ----------------------------!
         atm_up       = ups2(i,j)
         atm_vp       = vps2(i,j)
         atm_theta    = ths2(i,j)
         atm_rvap     = rvs2(i,j)
         atm_shv      = atm_rvap / (1. + atm_rvap)
         geoht        = zts2(i,j)
         atm_exner    = pis2(i,j)
         atm_co2      = co2s2(i,j)
         atm_prss     = atm_exner ** cpor * p00
         vels         = sqrt(atm_up ** 2 + atm_vp ** 2)
         gzotheta     = grav * geoht / atm_theta
         atm_temp     = atm_theta * atm_exner
         atm_enthalpy = ptqz2enthalpy(atm_prss,atm_temp,atm_shv,geoht)

         !----- Update water internal energy from time-dependent SST. ---------------------!
         leaf%soil_energy(mzg,i,j,1) =  cliq * (leaf%seatp(i,j)                            &
                                     + (leaf%seatf(i,j) - leaf%seatp(i,j)) * timefac_sst   &
                                     - tsupercool)

         !----- Fill surface precipitation arrays for input to LEAF-3 ---------------------!
         call sfc_pcp(nnqparm(ngrid),level,i,j,cuparm,micro)

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
            leaf%patch_area(i,j,1) = 1. - pctlcon
            leaf%patch_area(i,j,2) = pctlcon
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
            if (ip >= 2 .and. leaf%patch_area(i,j,ip) >= .009) then

               if (ip >= 2 .and. isfcl >= 1) then
                  call vegndvi( ngrid                      , leaf%patch_area  (i,j,ip)     &
                              , leaf%leaf_class(i,j,ip)    , leaf%veg_fracarea(i,j,ip)     &
                              , leaf%veg_lai   (i,j,ip)    , leaf%veg_tai     (i,j,ip)     &
                              , leaf%veg_rough (i,j,ip)    , leaf%veg_height  (i,j,ip)     &
                              , leaf%veg_albedo(i,j,ip)    , leaf%veg_ndvip   (i,j,ip)     &
                              , leaf%veg_ndvic (i,j,ip)    , leaf%veg_ndvif   (i,j,ip)     )
               end if
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

            !----- Set initial value of specific humidity. --------------------------------!
            can_shv = leaf%can_rvap(i,j,ip) / (leaf%can_rvap(i,j,ip) + 1.)

            !------------------------------------------------------------------------------!
            !    Update canopy air pressure here.  Canopy air pressure is assumed to       !
            ! remain constant during one LEAF full timestep, which means that heat equals  !
            ! to change in enthalpy.  Between time steps, atmospheric pressure changes so  !
            ! enthalpy is no longer a conserved variable.  Potential temperature, on the   !
            ! other hand, is conserved if no heat happens, which is the case between       !
            ! successive calls.  Therefore, we track theta outside LEAF, but enthalpy      !
            ! inside LEAF.                                                                 !
            !------------------------------------------------------------------------------!
            leaf%can_prss(i,j,ip) = reducedpress(atm_prss,atm_theta,atm_shv,geoht          &
                                                ,leaf%can_theta(i,j,ip),can_shv,can_depth)

            !----- Update canopy air temperature, density, and enthalpy. ------------------!
            can_temp     = leaf%can_theta(i,j,ip)                                          &
                         * (p00i * leaf%can_prss(i,j,ip)) ** rocp
            can_rhos     = idealdenssh (leaf%can_prss(i,j,ip),can_temp,can_shv)
            can_enthalpy = ptqz2enthalpy(leaf%can_prss(i,j,ip),can_temp,can_shv,can_depth)


            !----- Begin LEAF-3 small timestep here. --------------------------------------!
            tloop: do iter_leaf = 1,niter_leaf

               !----- Finding vegetation temperature and surface liquid water fraction. ---!
               call qwtk(leaf%veg_energy(i,j,ip),leaf%veg_water(i,j,ip)                    &
                        ,leaf%veg_hcap(i,j,ip),veg_temp,veg_fliq)

               !---------------------------------------------------------------------------!
               !    Calculate radiative fluxes between atmosphere, vegetation, and ground/ !
               ! /snow based on already-computed downward shortwave and longwave fluxes    !
               ! from the atmosphere.  Fill tempk array with soil and snow temperature (C) !
               ! and fracliq array with liquid fraction of water content in soil and snow. !
               ! Other snowcover properties are also computed here.                        !
               !---------------------------------------------------------------------------!

               if (iswrtyp > 0 .or. ilwrtyp > 0) then

                  if (ip == 1 .or. leaf%patch_area(i,j,ip) >= .009) then

                     ! If TEB is on, copy the values to single variables ------------------!
                     if (teb_spm == 1) then
                        g_urban   = leaf%g_urban(i,j,ip)
                        emis_town = tebc%emis_town(i,j)
                        alb_town  = tebc%alb_town(i,j)
                        ts_town   = tebc%ts_town(i,j)
                     end if

                     call sfcrad(mzg,mzs,ip              , leaf%soil_energy     (:,i,j,ip) &
                       , leaf%soil_water      (:,i,j,ip) , leaf%soil_text       (:,i,j,ip) &
                       , leaf%sfcwater_energy (:,i,j,ip) , leaf%sfcwater_mass   (:,i,j,ip) &
                       , leaf%sfcwater_depth  (:,i,j,ip) , leaf%patch_area      (  i,j,ip) &
                       , leaf%leaf_class      (  i,j,ip) , leaf%veg_height      (  i,j,ip) &
                       , leaf%veg_fracarea    (  i,j,ip) , leaf%veg_albedo      (  i,j,ip) &
                       , leaf%sfcwater_nlev   (  i,j,ip) , leaf%veg_energy      (  i,j,ip) &
                       , leaf%veg_water       (  i,j,ip) , leaf%veg_hcap        (  i,j,ip) &
                       , leaf%can_prss        (  i,j,ip) , leaf%can_theta       (  i,j,ip) &
                       , leaf%can_rvap        (  i,j,ip) , radiate%rshort       (  i,j   ) &
                       , radiate%rlong        (  i,j   ) , radiate%albedt       (  i,j   ) &
                       , radiate%rlongup      (  i,j   ) , radiate%cosz         (  i,j   ) &
                       , g_urban                         , emis_town                       &
                       , alb_town                        , ts_town                         )
                  end if
               end if

               !---------------------------------------------------------------------------!
               !    For water surface (patch 1), compute surface saturation mixing ratio   !
               ! and roughness length based on previous ustar.  For soil patches, compute  !
               ! roughness length based on vegetation and snow.                            !
               !---------------------------------------------------------------------------!
               if (ip == 1) then

                  leaf%ground_rsat(i,j,ip) = rslif(leaf%can_prss(i,j,ip),tempk(mzg))
                  leaf%patch_rough(i,j,ip) = max( z0fac_water * leaf%ustar(i,j,ip) ** 2    &
                                                , min_waterrough)

               elseif (isfcl >= 1) then

                  if (leaf%patch_area(i,j,ip) >= .009) then
                     leaf%patch_rough(i,j,ip) = max( grid%topzo(i,j)                       &
                                                   , leaf%soil_rough(i,j,ip)               &
                                                   , leaf%veg_rough(i,j,ip))               &
                                              * (1. - snowfac) + snowrough * snowfac
                  end if
               end if

               !---------------------------------------------------------------------------!
               !    Deciding the minimum velocity based on the stability.                  !
               !---------------------------------------------------------------------------!
               if (leaf%can_theta(i,j,ip) < atm_theta) then
                  !----- Stable case. -----------------------------------------------------!
                  vels_pat = max(ubmin_stab,vels)
               else
                  !----- Unstable case. ---------------------------------------------------!
                  vels_pat = max(ubmin_unstab,vels)
               end if

               !----- Compute the characteristic scales. ----------------------------------!
               call leaf_stars(atm_theta,atm_enthalpy,atm_shv,atm_rvap,atm_co2             &
                              ,leaf%can_theta(i,j,ip),can_enthalpy,can_shv                 &
                              ,leaf%can_rvap(i,j,ip),leaf%can_co2(i,j,ip)                  &
                              ,geoht,vels_pat,dtll                                         &
                              ,leaf%patch_rough(i,j,ip),leaf%ustar(i,j,ip)                 &
                              ,leaf%tstar(i,j,ip),estar,qstar,leaf%rstar(i,j,ip)           &
                              ,leaf%cstar(i,j,ip),leaf%R_aer(i,j,ip))


               if (teb_spm==1) g_urban = leaf%g_urban(i,j,ip)

               call sfclmcv( leaf%ustar        (i,j,ip) , leaf%tstar              (i,j,ip) &
                           , leaf%rstar        (i,j,ip) , leaf%cstar              (i,j,ip) &
                           , vels                       , vels_pat                         &
                           , atm_up                     , atm_vp                           &
                           , gzotheta                   , leaf%patch_area   (i,j,ip)       &
                           , turb%sflux_u      (i,j   ) , turb%sflux_v      (i,j   )       &
                           , turb%sflux_w      (i,j   ) , turb%sflux_t      (i,j   )       &
                           , turb%sflux_r      (i,j   ) , turb%sflux_c      (i,j   )       &
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
                  dtlloccc  = mmdry * dtllowcc

                  !----- Compute the fluxes from water body to canopy. --------------------!
                  hflxgc  = rdi * cp * (tempk(mzg) - can_temp)
                  wflxgc  = rdi * (leaf%ground_rsat(i,j,ip) - leaf%can_rvap(i,j,ip))
                  qwflxgc = wflxgc * alvl
                  cflxgc  = 0. !----- No water carbon emission model available...

                  !----- Compute the fluxes from atmosphere to canopy air space. ----------!
                  eflxac  = rho_ustar * estar
                  wflxac  = rho_ustar * leaf%rstar(i,j,ip)
                  cflxac  = rho_ustar * leaf%cstar(i,j,ip) * mmdryi

                  !----- Integrate the state variables. -----------------------------------!
                  can_enthalpy         = can_enthalpy                                      &
                                       + dtllowcc * (hflxgc + qwflxgc + eflxac)

                  leaf%can_rvap(i,j,ip) = leaf%can_rvap(i,j,ip)                            &
                                        + dtllowcc * (wflxgc + wflxac)

                  leaf%can_co2(i,j,ip)  = leaf%can_co2(i,j,ip)                             &
                                        + dtlloccc * (cflxgc + cflxac)

                  !----- Integrate the fluxes. --------------------------------------------!
                  leaf%sensible(i,j,ip) = leaf%sensible(i,j,ip) + hflxgc * dtll_factor
                  leaf%evap    (i,j,ip) = leaf%evap    (i,j,ip) + wflxgc * dtll_factor
                  leaf%plresp  (i,j,ip) = leaf%plresp  (i,j,ip) + cflxgc * dtll_factor

                  !------------------------------------------------------------------------!
                  !    Updating canopy air properties.  This is done inside the loop to    !
                  ! ensure that we leave here with the right canopy air potential temper-  !
                  ! ature.                                                                 !
                  !------------------------------------------------------------------------!
                  can_shv  = leaf%can_rvap(i,j,ip) / (leaf%can_rvap(i,j,ip) + 1.)
                  can_temp = hpqz2temp(can_enthalpy,leaf%can_prss(i,j,ip),can_shv,can_depth)
                  can_rhos = idealdenssh(leaf%can_prss(i,j,ip),can_temp,can_shv)
                  leaf%can_theta(i,j,ip) = can_temp * (p00 / leaf%can_prss(i,j,ip)) ** rocp

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
                     if (leaf%patch_area(i,j,ip) >= .009) then
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
                          , leaf%veg_water        (i,j,ip), leaf%veg_hcap         (i,j,ip) &
                          , leaf%veg_energy       (i,j,ip), leaf%can_prss         (i,j,ip) &
                          , leaf%can_theta        (i,j,ip), leaf%can_rvap         (i,j,ip) &
                          , leaf%can_co2          (i,j,ip), leaf%sensible         (i,j,ip) &
                          , leaf%evap             (i,j,ip), leaf%transp           (i,j,ip) &
                          , leaf%gpp              (i,j,ip), leaf%plresp           (i,j,ip) &
                          , leaf%resphet          (i,j,ip), leaf%veg_ndvip        (i,j,ip) &
                          , leaf%veg_ndvic        (i,j,ip), leaf%veg_ndvif        (i,j,ip) &
                          , radiate%rshort        (i,j)   , radiate%cosz          (i,j)    &
                          , ip,i,j)
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
                 ,veg_water,veg_hcap,veg_energy,can_prss,can_theta,can_rvap,can_co2        &
                 ,sensible,evap,transp,gpp,plresp,resphet,veg_ndvip,veg_ndvic,veg_ndvif    &
                 ,rshort,cosz,ip,i,j)

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
   integer                , intent(in)    :: mzg,mzs,np
   real   , dimension(mzg), intent(inout) :: soil_water,soil_energy,soil_text
   real   , dimension(mzs), intent(inout) :: sfcwater_mass,sfcwater_depth
   real                   , intent(in)    :: ustar,tstar,rstar,cstar
   real                   , intent(in)    :: veg_albedo,veg_fracarea,veg_lai,veg_tai
   real                   , intent(in)    :: veg_rough,veg_height,patch_area,patch_rough
   real                   , intent(in)    :: patch_wetind,leaf_class,soil_rough
   real                   , intent(inout) :: sfcwater_nlev, stom_resist
   real                   , intent(inout) :: ground_rsat,ground_rvap
   real                   , intent(inout) :: veg_water,veg_hcap,veg_energy
   real                   , intent(inout) :: can_prss,can_theta,can_rvap,can_co2
   real                   , intent(inout) :: sensible,evap,transp
   real                   , intent(inout) :: gpp,plresp,resphet
   real                   , intent(in)    :: veg_ndvip,veg_ndvic,veg_ndvif
   real                   , intent(in)    :: rshort,cosz
   !----- Local variables. ----------------------------------------------------------------!
   integer :: ip,k,nveg,ksn,nsoil,ksnnew,newlayers,nlayers,kold,ktrans,nsl,kzs
   integer :: i,j
   real    :: sfcw_energy_int,stretch,thik,pf,snden,vegfracc,qwfree,wfree
   real    :: depthgain,totsnow,qw,w,qwt,wt,soilhcap,fac,wfreeb,depthloss,soilcap
   real    :: sndenmax,sndenmin,snowmin,hold,wtnew,wtold,wdiff,watermid,availwat,wg
   real    :: wloss,soilcond,waterfrac
   !----- Locally saved variables. --------------------------------------------------------!
   logical               , save :: ncall=.true.
   real, dimension(20)   , save :: thicknet
   real, dimension(20,20), save :: thick
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
   !      Evaluate any exchanges of heat and moisture to or from vegetation, apply moist-  !
   ! ure and heat changes to vegetation, and evaluate the resistance parameter rd between  !
   ! canopy air and the top soil or snow surface.                                          !
   !---------------------------------------------------------------------------------------!
   call canopy(mzg,mzs,ksn,nveg,soil_energy,soil_water,soil_text,sfcwater_mass,ustar,tstar &
              ,rstar,cstar,soil_rough,veg_rough,veg_height,veg_lai,veg_tai,veg_water       &
              ,veg_hcap,veg_energy,leaf_class,veg_fracarea,stom_resist,can_prss,can_theta  &
              ,can_rvap,can_co2,sensible,evap,transp,gpp,plresp,resphet,ground_rsat        &
              ,ground_rvap,rshort,i,j,ip)

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
         sndenmin          = max(30.,200. * (wfree + wfreeb) / max(min_sfcwater_mass,sfcwater_mass(k)))
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
      soil_liq(k) = dslz(k) * min(soil_water(k) - soilcp(nsoil),soil_water(k)*fracliq(k))
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
      wflux(k) = dslztidt(k) * slcons1(k,nsoil)  &
               * (watermid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)  &
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
   !    Remove water only from moistest level in root zone for transpiration. Bottom       !
   ! k-index in root zone is kroot(nveg).  Limit soil moisture to values above soilcp.     !
   ! More sophisticated use of root profile function and water extractibility function,    !
   ! and improved minimum value are being considered. Units of wloss are m3/m3, of transp  !
   ! are kg/m2/s.                                                                          !
   !---------------------------------------------------------------------------------------!
   wg  = 0.
   nsl = nint(soil_text(mzg))
   ktrans = mzg
   do k = kroot(nveg),mzg
      nsoil = nint(soil_text(k))
      availwat = soil_water(k) * fracliq(k)
      if (wg < availwat) then
         wg = availwat
         ktrans = k
         nsl = nsoil
      end if
   end do

   wloss = min(transp_tot * dslzidt(ktrans) * wdnsi, wg, soil_water(ktrans) - soilcp(nsl))
   soil_water(ktrans)  = soil_water(ktrans)  - wloss
   soil_energy(ktrans) = soil_energy(ktrans)                                               &
                       - wloss * cliqvlme * (tempk(ktrans) - tsupercool)

   !---------------------------------------------------------------------------------------!
   ! Compute ground vap mxrat for availability on next timestep; put into ground_rsat.     !
   !---------------------------------------------------------------------------------------!
   newlayers = max(1,nint(sfcwater_nlev))
   if (sfcwater_mass(newlayers) > min_sfcwater_mass) then
      sfcw_energy_int = sfcw_energy_ext(newlayers) / sfcwater_mass(newlayers)
   else
      sfcw_energy_int = 0.
   end if
   call grndvap(soil_energy(mzg),soil_water(mzg),soil_text(mzg),sfcw_energy_int            &
               ,sfcwater_nlev,ground_rsat,ground_rvap,can_rvap,can_prss)

   return
end subroutine leaftw
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the variables that depend on heat, water, and (eventually)  !
! carbon exchanges with the canopy air space and vegetation.                               !
!------------------------------------------------------------------------------------------!
subroutine canopy(mzg,mzs,ksn,nveg,soil_energy,soil_water,soil_text,sfcwater_mass,ustar    &
                 ,tstar,rstar,cstar,soil_rough,veg_rough,veg_height,veg_lai,veg_tai        &
                 ,veg_water,veg_hcap,veg_energy,leaf_class,veg_fracarea,stom_resist        &
                 ,can_prss,can_theta,can_rvap,can_co2,sensible,evap,transp,gpp,plresp      &
                 ,resphet,ground_rsat,ground_rvap,rshort,i,j,ip)
   use leaf_coms
   use rconstants
   use therm_lib , only : eslif       & ! function
                        , rslif       & ! function
                        , hpqz2temp   & ! function
                        , idealdenssh & ! function
                        , qwtk        ! ! subroutine
   use catt_start, only : CATT
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer             , intent(in)    :: mzg,mzs,ksn
   integer             , intent(inout) :: nveg
   real, dimension(mzg), intent(in)    :: soil_energy,soil_water,soil_text
   real, dimension(mzs), intent(in)    :: sfcwater_mass
   real                , intent(in)    :: ustar,tstar,rstar,cstar
   real                , intent(in)    :: soil_rough,veg_rough,veg_height,veg_lai,veg_tai
   real                , intent(in)    :: leaf_class,veg_fracarea
   real                , intent(in)    :: ground_rsat,ground_rvap,rshort
   real                , intent(inout) :: veg_water,veg_hcap,veg_energy
   real                , intent(inout) :: stom_resist,can_prss,can_theta,can_rvap,can_co2
   real                , intent(inout) :: sensible,evap,transp,gpp,plresp,resphet
   !----- Local arguments. ----------------------------------------------------------------!
   integer                             :: k,kk,nsoil,iter_can
   integer                             :: i,j,ip
   real                                :: aux,c2,c3,c4,dsm
   real                                :: es,fac,factv
   real                                :: fthi,ftlo,frad,fswp,fvpd,qwtot
   real                                :: rasgnd,rasveg,rleaf,rsatveg
   real                                :: sigmaw ,slai,stai,slpotv,swp
   real                                :: vpd,wtemp,wtroot,x,zognd,zoveg,zdisp,zveg
   real                                :: wflx,dewgndflx,qdewgndflx
   real                                :: transp_test,rc,rc_inf
   real                                :: wshed_loc,qwshed_loc


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

   hflxgc     = cp * rdi * (tempk(mzg+ksn) - can_temp)
   wflx       =      rdi * (ground_rsat    - can_rvap)
   dewgndflx  = max(0.,-wflx)
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
      rb  = (1. + .5 * veg_tai) / (.01 * sqrt(ustar * c1))

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
      !      Begin canopy time-split iterations.                                           !
      !------------------------------------------------------------------------------------!
      do iter_can = 1,niter_can

         !----- Calculate the saturation vapor pressure at TVEG. --------------------------!
         es      = eslif(veg_temp)
         rsatveg = rslif(can_prss,veg_temp)
         
         !----- Compute mixing ratio at leaf surface using previous rc --------------------!
         rc    = stom_resist
         rleaf = (rb * rsatveg + rc * can_rvap) / (rb + rc)
         vpd   = max((es - rleaf * can_prss / (ep + rleaf)),0.)

         !----- Evaluate 5 environmental factors and new rc -------------------------------!
         ftlo = 1. + exp(-stlo * (veg_temp - btlo))
         fthi = 1. + exp(-sthi * (veg_temp - bthi))
         frad = 1. + exp(-srad * (rshort   - brad))
         fswp = 1. + exp(-sswp * (swp      - bswp))
         fvpd = 1. + exp(-svpd * (vpd      - bvpd))

         !----- 15-minute response time for stomatal conductance (must have dtlc <= 900.). !
         rc_inf = ftlo * fthi * frad * fvpd * fswp * rcmin(nveg)
         rc     = 1. / (1. / rc + .0011 * dtlc * (1. / rc_inf - 1. / rc))

         !---------------------------------------------------------------------------------!
         !    Limit maximum transpiration to be <= 400 w m-2 by increasing rc if           !
         ! necessary.                                                                      !
         !---------------------------------------------------------------------------------!
         transp_test = alvl * can_rhos * veg_lai * (rsatveg - can_rvap) / (rb + rc)
         if (transp_test > 400.) then
            rc = (rb + rc) * transp_test * .0025 - rb
         endif
         stom_resist = rc

         !----- Flux of heat and moisture from vegetation to canopy. ----------------------!
         stai = veg_tai * (1. - snowfac)
         slai = veg_lai * (1. - snowfac)

         !----- Sensible heat flux from leaf to canopy. -----------------------------------!
         hflxvc = 2. * stai * cp * can_rhos * (veg_temp - can_temp) / rb

         c3     = can_rhos * (rsatveg - can_rvap)

         if (c3 >= 0.) then
            !----- Flow will go towards the canopy, allow evaporation and transpiration. --!
            sigmaw     = min(1.,(veg_water / (.2 * stai)) ** .66667)
            wflxvc     = min(c3 * 2. * stai * sigmaw / rb,veg_water/dtlc)
            transp_loc = max(0.,c3 * slai * (1.-sigmaw) / (rb + rc))
         else
            !----- Flow is towards the leaf, dew/frost formation and no transpiration. ----!
            wflxvc     = c3 * 2. * stai / rb
            transp_loc = 0.
         end if
         transp_tot = transp_tot + transp_loc

         !---------------------------------------------------------------------------------!
         !     Update vegetation moisture from precipitation, vapor flux with canopy,      !
         ! and shedding of excess moisture.                                                !
         !---------------------------------------------------------------------------------!
         wtemp     = veg_water - wflxvc * dtlc + pcpgc * vf
         !----- If this result will lead to negative solution, then reduce evaporation. ---!
         if (wtemp < 0.) then
            wtemp     = 0.
            wshed_loc = 0.
            wflxvc    = (veg_water + pcpgc * vf) / dtlc
         else
            wshed_loc = max(0.,wtemp - 0.2 * stai)
            wshed_tot = wshed_tot + wshed_loc
            veg_water = wtemp - wshed_loc
         end if
         !---------------------------------------------------------------------------------!

         !------ Find the associated latent heat flux from vegetation to canopy. ----------!
         qwflxvc     = wflxvc     * (alvi - alli * veg_fliq)
         qtransp_loc = transp_loc *  alvl !----- Liquid phase only in transpiration. ------!
         !---------------------------------------------------------------------------------!


         !----- Exchange heat between vegetation and precipitation in implicit scheme -----!
         qwshed_loc = wshed_loc * (     veg_fliq  * cliq * (veg_temp - tsupercool)         &
                                  + (1.-veg_fliq) * cice *  veg_temp )
         qwshed_tot = qwshed_tot + qwshed_loc

         !---------------------------------------------------------------------------------!
         !      Update vegetation internal energy from radiative fluxes and sensible and   !
         ! latent heat transfer with canopy air.                                           !
         !---------------------------------------------------------------------------------!
         veg_energy = veg_energy                                                           &
                    + dtlc * ( rshort_v + rlonga_v + rlonggs_v - rlongv_gs - rlongv_a      &
                             - hflxvc   - qwflxvc   - qtransp_loc)                         &
                    + qpcpgc * vf - qwshed_loc

         !---------------------------------------------------------------------------------!
         !     Find the atmosphere -> canopy fluxes.                                       !
         !---------------------------------------------------------------------------------!
         rho_ustar  = can_rhos  * ustar
         eflxac     = rho_ustar * estar          ! Enthalpy exchange
         wflxac     = rho_ustar * rstar          ! Water vapour exchange
         cflxac     = rho_ustar * mmdryi * cstar ! CO2 exchange
         
         dtlcowcc = dtlc / (can_depth * can_rhos)
         dtlcoccc = dtlcowcc * mmdry

         !---------------------------------------------------------------------------------!
         !      Update enthalpy, CO2, and canopy mixing ratio.                             !
         !---------------------------------------------------------------------------------!
         can_enthalpy = can_enthalpy                                                       &
                      + dtlcowcc * ( hflxgc + hflxvc + eflxac                              &
                                   + qwflxgc - qdewgndflx + qwflxvc + qtransp_loc)

         can_rvap     = can_rvap                                                           &
                      + dtlcowcc * (wflxgc - dewgndflx + wflxvc + transp_loc + wflxac)
         !---------------------------------------------------------------------------------!
         !      CO2 exchange with biosphere is currently not computed in LEAF-3... Feel    !
         ! free to include here: cflxgc would be the root respiration, and cflxvc R-GPP.   !
         !---------------------------------------------------------------------------------!
         cflxgc   = 0.
         cflxvc   = 0.
         cflxcv   = 0.
         can_co2  = can_co2  + dtlcowcc * mmdry * (cflxgc + cflxvc + cflxcv + cflxac)

         !----- Update the fluxes. --------------------------------------------------------!
         sensible = sensible + (hflxvc  + hflxgc ) * dtlc_factor
         evap     = evap     + (qwflxvc + qwflxgc) * dtlc_factor
         transp   = transp   + (qtransp_loc      ) * dtlc_factor
         ! gpp      = gpp     + cflxcv * dtlc_factor
         ! plresp   = plresp  + cflxvc * dtlc_factor
         ! resphet  = resphet + cflxgc * dtlc_factor

         !----- Find the vegetation diagnostic variables for next time step. --------------!
         call qwtk(veg_energy,veg_water,veg_hcap,veg_temp,veg_fliq)

         !----- Find the canopy air diagnostic variables for next time step. --------------!
         can_shv   = can_rvap / (can_rvap + 1.)
         can_temp  = hpqz2temp(can_enthalpy,can_prss,can_shv,can_depth)
         can_rhos  = idealdenssh(can_prss,can_temp,can_shv)
         can_theta = can_temp * (p00 / can_prss) ** rocp

      end do

   else
      !------------------------------------------------------------------------------------!
      !    No vegetation, or vegetation buried in snow... Update the canopy properties and !
      ! move to the next time step without nesting.                                        !
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Find the atmosphere -> canopy fluxes.                                          !
      !------------------------------------------------------------------------------------!
      rho_ustar  = can_rhos  * ustar
      eflxac     = rho_ustar * estar          ! Enthalpy exchange
      wflxac     = rho_ustar * rstar          ! Water vapour exchange
      cflxac     = rho_ustar * mmdryi * cstar ! CO2 exchange
      !------------------------------------------------------------------------------------!

      dtllowcc = dtll / (can_depth * can_rhos)
      dtlloccc = dtllowcc * mmdry

      !----- Update the canopy prognostic variables. --------------------------------------!
      can_enthalpy = can_enthalpy + dtllowcc * (hflxgc + qwflxgc - qdewgndflx + eflxac)
      can_rvap     = can_rvap     + dtllowcc * (wflxgc - dewgndflx + wflxac)
      !------------------------------------------------------------------------------------!
      !      CO2 exchange with biosphere is currently not computed in LEAF-3... Feel       !
      ! free to include here: cflxgc would be the root respiration, and cflxvc R-GPP.      !
      !------------------------------------------------------------------------------------!
      cflxgc       = 0.
      can_co2      = can_co2      + dtlcoccc * (cflxgc + cflxac)
      !------------------------------------------------------------------------------------!

      !----- Update the fluxes. -----------------------------------------------------------!
      sensible = sensible + hflxgc  * dtll_factor
      evap     = evap     + qwflxgc * dtll_factor
      ! resphet  = resphet  + cflxgc  * dtll_factor

      !----- Find the canopy air diagnostic variables for next time step. -----------------!
      can_shv   = can_rvap / (can_rvap + 1.)
      can_temp  = hpqz2temp(can_enthalpy,can_prss,can_shv,can_depth)
      can_rhos  = idealdenssh(can_prss,can_temp,can_shv)
      can_theta = can_temp * (p00 / can_prss) ** rocp
      !------------------------------------------------------------------------------------!
   end if

   return
end subroutine canopy
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the characteristic scales, using the surface layer          !
! parameterization proposed by:                                                            !
! LOUIS, J.F., Parametric Model of vertical eddy fluxes in the atmosphere. Boundary-Layer  !
!     Meteor., 17, 187-202, 1979.                                                          !
!------------------------------------------------------------------------------------------!
subroutine leaf_stars(theta_atm,enthalpy_atm,shv_atm,rvap_atm,co2_atm,theta_can            &
                     ,enthalpy_can,shv_can,rvap_can,co2_can,zref,uref,dtll,rough           &
                     ,ustar,tstar,estar,qstar,rstar,cstar,r_aer)
   use rconstants, only : grav          & ! intent(in)
                        , vonk          ! ! intent(in)
   use leaf_coms , only : ustmin_stab   & ! intent(in)
                        , ustmin_unstab ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   real, intent(in)  :: theta_atm    ! Above canopy air pot. temperature        [        K]
   real, intent(in)  :: enthalpy_atm ! Above canopy air enthalpy                [     J/kg]
   real, intent(in)  :: shv_atm      ! Above canopy vapour spec. hum.           [kg/kg_air]
   real, intent(in)  :: rvap_atm     ! Above canopy vapour mixing ratio         [kg/kg_air]
   real, intent(in)  :: co2_atm      ! Above canopy CO2 mixing ratio            [ µmol/mol]
   real, intent(in)  :: theta_can    ! Canopy air potential temperature         [        K]
   real, intent(in)  :: enthalpy_can ! Canopy air enthalpy                      [     J/kg]
   real, intent(in)  :: shv_can      ! Canopy air vapour spec. humidity         [kg/kg_air]
   real, intent(in)  :: rvap_can     ! Canopy air vapour mixing ratio           [kg/kg_air]
   real, intent(in)  :: co2_can      ! Canopy air CO2 mixing ratio              [ µmol/mol]
   real, intent(in)  :: zref         ! Height at reference point                [        m]
   real, intent(in)  :: uref         ! Wind speed at reference height           [      m/s]
   real, intent(in)  :: dtll         ! Time step                                [        m]
   real, intent(in)  :: rough        ! Roughness                                [        m]
   real, intent(out) :: ustar        ! U*, friction velocity                    [      m/s]
   real, intent(out) :: tstar        ! Temperature friction scale               [        K]
   real, intent(out) :: estar        ! Enthalpy friction scale                  [     J/kg]
   real, intent(out) :: qstar        ! Specific humidity friction scale         [kg/kg_air]
   real, intent(out) :: rstar        ! Vapour mixing ratio friction scale       [kg/kg_air]
   real, intent(out) :: cstar        ! CO2 mixing ratio                         [ µmol/mol]
   real, intent(out) :: r_aer        ! Aerodynamic resistance                   [      s/m]
   !----- Local variables. ----------------------------------------------------------------!
   real              :: a2           ! Drag coefficient in neutral conditions, 
   real              :: c1           !
   real              :: ri           ! Bulk richardson numer, eq. 3.45 in Garratt
   real              :: fh           ! Stability parameter for heat
   real              :: fm           ! Stability parameter for momentum
   real              :: c2           !
   real              :: cm           !
   real              :: ch           !
   real              :: c3           !
   real              :: delz         !
   real              :: d_vel        !
   real              :: vel_new      !
   !----- Local constants. ----------------------------------------------------------------!
   real, parameter   :: b   = 5.0    !
   real, parameter   :: csm = 7.5    !
   real, parameter   :: csh = 5.0    !
   real, parameter   :: d   = 5.0    !
   !---------------------------------------------------------------------------------------!


   !----- Make the log profile (constant shear assumption). -------------------------------!
   a2 = (vonk / log(zref/rough)) ** 2.
   c1 = a2 * uref
   ri = grav * zref * (theta_atm-theta_can) / (0.5 * (theta_atm+theta_can) * uref * uref )

   if (theta_atm - theta_can > 0.0) then
      !----- Stable case ------------------------------------------------------------------!
      fm = 1.0 / (1.0 + (2.0 * b * ri / sqrt(1.0 + d * ri)))
      fh = 1.0 / (1.0 + (3.0 * b * ri * sqrt(1.0 + d * ri)))

      !----- Making sure ustar is not too small. ------------------------------------------!
      ustar = max(ustmin_stab,sqrt(c1 * uref * fm))

   else
      !----- Unstable case ----------------------------------------------------------------!

      c2 = b * a2 * sqrt(zref/rough * (abs(ri)))
      cm = csm * c2
      ch = csh * c2
      fm = (1.0 - 2.0 * b * ri / (1.0 + 2.0 * cm))
      fh = (1.0 - 3.0 * b * ri / (1.0 + 3.0 * ch))

      !----- Making sure ustar is not too small. ------------------------------------------!
      ustar = max(ustmin_unstab,sqrt(c1 * uref * fm))
   end if

   c3 = c1 * fh / ustar
   tstar = c3 * (theta_atm    - theta_can   )
   estar = c3 * (enthalpy_atm - enthalpy_can)
   qstar = c3 * (shv_atm      - shv_can     )
   rstar = c3 * (rvap_atm     - rvap_can    )
   cstar = c3 * (co2_atm      - co2_can     )
   
   if (abs(tstar) < 1.e-7) tstar = 0.
   if (abs(estar) < 1.e-7) estar = 0.
   if (abs(qstar) < 1.e-7) qstar = 0.
   if (abs(rstar) < 1.e-7) rstar = 0.
   if (abs(cstar) < 1.e-7) cstar = 0.

   r_aer = 1. / (a2 * uref * fh)

   !----- Limit ustar so that the flux cannot take more than 1/2 velocity in a timestep ---!
   delz  = 2. * zref
   d_vel = - ustar * ustar * dtll / delz
   vel_new = uref + d_vel
   if (vel_new < .5 * uref) then
      d_vel = .5 * uref
      ustar = sqrt(d_vel * delz / dtll)
   end if

   return
end subroutine leaf_stars
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine computes the turbulent fluxes of momentum, heat and moisture from the    !
! surface layer using the  Manton-Cotton algebraic surface layer equations.                !
!------------------------------------------------------------------------------------------!
subroutine sfclmcv(ustar,tstar,rstar,cstar,vels,vels_pat,ups,vps,gzotheta,patch_area       &
                  ,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r,sflux_c,g_urban)
   use rconstants
   use teb_spm_start , only : teb_spm ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real , intent(in)    :: ustar,tstar,rstar,cstar
   real , intent(in)    :: vels,vels_pat,ups,vps,gzotheta,patch_area
   real , intent(inout) :: sflux_u,sflux_v,sflux_w,sflux_t,sflux_r,sflux_c
   real , intent(inout) :: g_urban
   !----- Local variables. ----------------------------------------------------------------!
   real                 :: zoverl
   real                 :: cosine1,sine1,vtscr,cx,psin
   !----- Local constants. ----------------------------------------------------------------!
   real , parameter     :: wtol = 1.e-20
   !---------------------------------------------------------------------------------------!

   cosine1 = ups / vels_pat
   sine1   = vps / vels_pat

   vtscr = ustar * patch_area

   !----- Check whether TEB is used. If it is, skip the lines unless g_urban is zero. -----!
   if (teb_spm /= 1 .or. nint(g_urban) == 0) then
      sflux_u = sflux_u - ustar * vtscr * cosine1
      sflux_v = sflux_v - ustar * vtscr * sine1
      sflux_t = sflux_t - tstar * vtscr
      sflux_r = sflux_r - rstar * vtscr
   end if

   !----- TEB currently doesn't save CO2, so compute sflux_c outside the if statement. ----!
   sflux_c = sflux_c - cstar * vtscr


   !----- Compute the vertical flux. ------------------------------------------------------!
   zoverl = gzotheta * vonk * tstar / (ustar * ustar)

   if (zoverl < 0.)then
      cx = zoverl * sqrt(sqrt(1. - 15. * zoverl))
   else
      cx = zoverl / (1.0 + 4.7 * zoverl)
   end if

   psin    = sqrt((1.-2.86 * cx) / (1. + cx * (-5.39 + cx * 6.998 )))
   sflux_w = sflux_w + (0.27 * max(6.25 * (1. - cx) * psin,wtol) - 1.18 * cx * psin)       &
                     * ustar * vtscr

   return
end subroutine sfclmcv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine grndvap(soil_energy,soil_water,soil_text,sfcw_energy_int  &
   ,sfcwater_nlev,ground_rsat,ground_rvap,can_rvap,prsg)

  use leaf_coms
  use rconstants
  use therm_lib, only: rslif,qwtk,qtk

  implicit none

  real :: soil_energy,soil_water,soil_text,sfcw_energy_int
  real :: sfcwater_nlev,ground_rsat,ground_rvap,can_rvap,prsg

  integer :: ksn,nsoil

  real :: slpotvn,alpha,beta


  ksn = nint(sfcwater_nlev)

  ! ground_rsat(i,j) is the saturation mixing ratio of the top soil/snow surface
  ! and is used for dew formation and snow evaporation.

  if (ksn > 0 .and. sfcw_energy_int > 0.) then
     call qtk(sfcw_energy_int,ground_temp,ground_fliq)
     ground_rsat = rslif(prsg,ground_temp)
  else

  ! Without snowcover, ground_rvap is the effective saturation mixing
  ! ratio of soil and is used for soil evaporation.  First, compute the
  ! "alpha" term or soil "relative humidity" and the "beta" term.

     nsoil = nint(soil_text)

     call qwtk(soil_energy,soil_water*wdns,slcpd(nsoil),ground_temp,ground_fliq)
     ground_rsat = rslif(prsg,ground_temp)

     slpotvn = slpots(nsoil) * (slmsts(nsoil) / soil_water) ** slbs(nsoil)
     alpha = exp(gorh2o * slpotvn / ground_temp)
     beta = .25 * (1. - cos (min(1.,soil_water / sfldcap(nsoil)) * pi1)) ** 2
     ground_rvap = ground_rsat * alpha * beta + (1. - beta) * can_rvap
  end if

  return
end subroutine grndvap
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will aplly the boundary condition to all leaf variables at the        !
! absolute domain edges.                                                                   !
!------------------------------------------------------------------------------------------!
subroutine leaf_bcond(m2,m3,mzg,mzs,npat,jdim,soil_water,sfcwater_mass,soil_energy         &
                     ,sfcwater_energy,soil_text,sfcwater_depth,ustar,tstar,rstar,cstar     &
                     ,veg_albedo,veg_fracarea,veg_lai,veg_tai,veg_rough,veg_height         &
                     ,patch_area,patch_rough,patch_wetind,leaf_class,soil_rough            &
                     ,sfcwater_nlev,stom_resist,ground_rsat,ground_rvap,veg_water          &
                     ,veg_hcap,veg_energy,can_prss,can_theta,can_rvap,can_co2,sensible     &
                     ,evap,transp,gpp,plresp,resphet,veg_ndvip,veg_ndvic,veg_ndvif)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                        , intent(in)    :: m2,m3,mzg,mzs,npat,jdim
   real, dimension(mzg,m2,m3,npat), intent(inout) :: soil_water,soil_energy,soil_text
   real, dimension(mzs,m2,m3,npat), intent(inout) :: sfcwater_mass,sfcwater_energy
   real, dimension(mzs,m2,m3,npat), intent(inout) :: sfcwater_depth
   real, dimension(m2,m3,npat)    , intent(inout) :: ustar,tstar,rstar,cstar
   real, dimension(m2,m3,npat)    , intent(inout) :: veg_albedo,veg_fracarea
   real, dimension(m2,m3,npat)    , intent(inout) :: veg_lai,veg_tai,veg_rough,veg_height
   real, dimension(m2,m3,npat)    , intent(inout) :: patch_area,patch_rough,patch_wetind
   real, dimension(m2,m3,npat)    , intent(inout) :: leaf_class,soil_rough,sfcwater_nlev
   real, dimension(m2,m3,npat)    , intent(inout) :: stom_resist,ground_rsat,ground_rvap
   real, dimension(m2,m3,npat)    , intent(inout) :: veg_water,veg_energy,veg_hcap
   real, dimension(m2,m3,npat)    , intent(inout) :: can_prss,can_theta,can_rvap,can_co2
   real, dimension(m2,m3,npat)    , intent(inout) :: sensible,evap,transp
   real, dimension(m2,m3,npat)    , intent(inout) :: gpp,plresp,resphet
   real, dimension(m2,m3,npat)    , intent(inout) :: veg_ndvip,veg_ndvic,veg_ndvif
   !----- Local variables. ----------------------------------------------------------------!
   integer                                        :: i,j,k,ipat
   !---------------------------------------------------------------------------------------!
   do ipat = 1,npat
      do j = 1,m3

         ustar          (1,j,ipat) = ustar            (2,j,ipat)
         tstar          (1,j,ipat) = tstar            (2,j,ipat)
         rstar          (1,j,ipat) = rstar            (2,j,ipat)
         cstar          (1,j,ipat) = cstar            (2,j,ipat)
         veg_fracarea   (1,j,ipat) = veg_fracarea     (2,j,ipat)
         veg_lai        (1,j,ipat) = veg_lai          (2,j,ipat)
         veg_tai        (1,j,ipat) = veg_tai          (2,j,ipat)
         veg_rough      (1,j,ipat) = veg_rough        (2,j,ipat)
         veg_height     (1,j,ipat) = veg_height       (2,j,ipat)
         patch_area     (1,j,ipat) = patch_area       (2,j,ipat)
         patch_rough    (1,j,ipat) = patch_rough      (2,j,ipat)
         patch_wetind   (1,j,ipat) = patch_wetind     (2,j,ipat)
         leaf_class     (1,j,ipat) = leaf_class       (2,j,ipat)
         soil_rough     (1,j,ipat) = soil_rough       (2,j,ipat)
         sfcwater_nlev  (1,j,ipat) = sfcwater_nlev    (2,j,ipat)
         stom_resist    (1,j,ipat) = stom_resist      (2,j,ipat)
         ground_rsat    (1,j,ipat) = ground_rsat      (2,j,ipat)
         ground_rvap    (1,j,ipat) = ground_rvap      (2,j,ipat)
         veg_water      (1,j,ipat) = veg_water        (2,j,ipat)
         veg_hcap       (1,j,ipat) = veg_hcap         (2,j,ipat)
         veg_energy     (1,j,ipat) = veg_energy       (2,j,ipat)
         can_prss       (1,j,ipat) = can_prss         (2,j,ipat)
         can_theta      (1,j,ipat) = can_theta        (2,j,ipat)
         can_rvap       (1,j,ipat) = can_rvap         (2,j,ipat)
         can_co2        (1,j,ipat) = can_co2          (2,j,ipat)
         sensible       (1,j,ipat) = sensible         (2,j,ipat)
         evap           (1,j,ipat) = evap             (2,j,ipat)
         transp         (1,j,ipat) = transp           (2,j,ipat)
         gpp            (1,j,ipat) = gpp              (2,j,ipat)
         plresp         (1,j,ipat) = plresp           (2,j,ipat)
         resphet        (1,j,ipat) = resphet          (2,j,ipat)
         veg_ndvip      (1,j,ipat) = veg_ndvip        (2,j,ipat)
         veg_ndvic      (1,j,ipat) = veg_ndvic        (2,j,ipat)
         veg_ndvif      (1,j,ipat) = veg_ndvif        (2,j,ipat)

         ustar         (m2,j,ipat) = ustar         (m2-1,j,ipat)
         tstar         (m2,j,ipat) = tstar         (m2-1,j,ipat)
         rstar         (m2,j,ipat) = rstar         (m2-1,j,ipat)
         cstar         (m2,j,ipat) = cstar         (m2-1,j,ipat)
         veg_albedo    (m2,j,ipat) = veg_albedo    (m2-1,j,ipat)
         veg_fracarea  (m2,j,ipat) = veg_fracarea  (m2-1,j,ipat)
         veg_lai       (m2,j,ipat) = veg_lai       (m2-1,j,ipat)
         veg_tai       (m2,j,ipat) = veg_tai       (m2-1,j,ipat)
         veg_rough     (m2,j,ipat) = veg_rough     (m2-1,j,ipat)
         veg_height    (m2,j,ipat) = veg_height    (m2-1,j,ipat)
         patch_area    (m2,j,ipat) = patch_area    (m2-1,j,ipat)
         patch_rough   (m2,j,ipat) = patch_rough   (m2-1,j,ipat)
         patch_wetind  (m2,j,ipat) = patch_wetind  (m2-1,j,ipat)
         leaf_class    (m2,j,ipat) = leaf_class    (m2-1,j,ipat)
         soil_rough    (m2,j,ipat) = soil_rough    (m2-1,j,ipat)
         sfcwater_nlev (m2,j,ipat) = sfcwater_nlev (m2-1,j,ipat)
         stom_resist   (m2,j,ipat) = stom_resist   (m2-1,j,ipat)
         ground_rsat   (m2,j,ipat) = ground_rsat   (m2-1,j,ipat)
         ground_rvap   (m2,j,ipat) = ground_rvap   (m2-1,j,ipat)
         veg_water     (m2,j,ipat) = veg_water     (m2-1,j,ipat)
         veg_hcap      (m2,j,ipat) = veg_hcap      (m2-1,j,ipat)
         veg_energy    (m2,j,ipat) = veg_energy    (m2-1,j,ipat)
         can_prss      (m2,j,ipat) = can_prss      (m2-1,j,ipat)
         can_theta     (m2,j,ipat) = can_theta     (m2-1,j,ipat)
         can_rvap      (m2,j,ipat) = can_rvap      (m2-1,j,ipat)
         can_co2       (m2,j,ipat) = can_co2       (m2-1,j,ipat)
         sensible      (m2,j,ipat) = sensible      (m2-1,j,ipat)
         evap          (m2,j,ipat) = evap          (m2-1,j,ipat)
         transp        (m2,j,ipat) = transp        (m2-1,j,ipat)
         gpp           (m2,j,ipat) = gpp           (m2-1,j,ipat)
         plresp        (m2,j,ipat) = plresp        (m2-1,j,ipat)
         resphet       (m2,j,ipat) = resphet       (m2-1,j,ipat)
         veg_ndvip     (m2,j,ipat) = veg_ndvip     (m2-1,j,ipat)
         veg_ndvic     (m2,j,ipat) = veg_ndvic     (m2-1,j,ipat)
         veg_ndvif     (m2,j,ipat) = veg_ndvif     (m2-1,j,ipat)

         do k = 1,mzg
            soil_water       (k,1,j,ipat) = soil_water         (k,2,j,ipat)
            soil_energy      (k,1,j,ipat) = soil_energy        (k,2,j,ipat)
            soil_text        (k,1,j,ipat) = soil_text          (k,2,j,ipat)

            soil_water      (k,m2,j,ipat) = soil_water      (k,m2-1,j,ipat)
            soil_energy     (k,m2,j,ipat) = soil_energy     (k,m2-1,j,ipat)
            soil_text       (k,m2,j,ipat) = soil_text       (k,m2-1,j,ipat)
         end do

         do k = 1,mzs
            sfcwater_mass    (k,1,j,ipat) = sfcwater_mass      (k,2,j,ipat)
            sfcwater_energy  (k,1,j,ipat) = sfcwater_energy    (k,2,j,ipat)
            sfcwater_depth   (k,1,j,ipat) = sfcwater_depth     (k,2,j,ipat)

            sfcwater_mass   (k,m2,j,ipat) = sfcwater_mass   (k,m2-1,j,ipat)
            sfcwater_energy (k,m2,j,ipat) = sfcwater_energy (k,m2-1,j,ipat)
            sfcwater_depth  (k,m2,j,ipat) = sfcwater_depth  (k,m2-1,j,ipat)
         end do
      end do


      if (jdim == 1) then
         do i = 1,m2
            ustar          (i,1,ipat) = ustar            (i,2,ipat)
            tstar          (i,1,ipat) = tstar            (i,2,ipat)
            rstar          (i,1,ipat) = rstar            (i,2,ipat)
            cstar          (i,1,ipat) = cstar            (i,2,ipat)
            veg_albedo     (i,1,ipat) = veg_albedo       (i,2,ipat)
            veg_fracarea   (i,1,ipat) = veg_fracarea     (i,2,ipat)
            veg_lai        (i,1,ipat) = veg_lai          (i,2,ipat)
            veg_tai        (i,1,ipat) = veg_tai          (i,2,ipat)
            veg_rough      (i,1,ipat) = veg_rough        (i,2,ipat)
            veg_height     (i,1,ipat) = veg_height       (i,2,ipat)
            patch_area     (i,1,ipat) = patch_area       (i,2,ipat)
            patch_rough    (i,1,ipat) = patch_rough      (i,2,ipat)
            patch_wetind   (i,1,ipat) = patch_wetind     (i,2,ipat)
            leaf_class     (i,1,ipat) = leaf_class       (i,2,ipat)
            soil_rough     (i,1,ipat) = soil_rough       (i,2,ipat)
            sfcwater_nlev  (i,1,ipat) = sfcwater_nlev    (i,2,ipat)
            stom_resist    (i,1,ipat) = stom_resist      (i,2,ipat)
            ground_rsat    (i,1,ipat) = ground_rsat      (i,2,ipat)
            ground_rvap    (i,1,ipat) = ground_rvap      (i,2,ipat)
            veg_water      (i,1,ipat) = veg_water        (i,2,ipat)
            veg_hcap       (i,1,ipat) = veg_hcap         (i,2,ipat)
            veg_energy     (i,1,ipat) = veg_energy       (i,2,ipat)
            can_prss       (i,1,ipat) = can_prss         (i,2,ipat)
            can_theta      (i,1,ipat) = can_theta        (i,2,ipat)
            can_rvap       (i,1,ipat) = can_rvap         (i,2,ipat)
            can_co2        (i,1,ipat) = can_co2          (i,2,ipat)
            sensible       (i,1,ipat) = sensible         (i,2,ipat)
            evap           (i,1,ipat) = evap             (i,2,ipat)
            transp         (i,1,ipat) = transp           (i,2,ipat)
            gpp            (i,1,ipat) = gpp              (i,2,ipat)
            plresp         (i,1,ipat) = plresp           (i,2,ipat)
            resphet        (i,1,ipat) = resphet          (i,2,ipat)
            veg_ndvip      (i,1,ipat) = veg_ndvip        (i,2,ipat)
            veg_ndvic      (i,1,ipat) = veg_ndvic        (i,2,ipat)
            veg_ndvif      (i,1,ipat) = veg_ndvif        (i,2,ipat)

            ustar         (i,m3,ipat) = ustar         (i,m3-1,ipat)
            tstar         (i,m3,ipat) = tstar         (i,m3-1,ipat)
            rstar         (i,m3,ipat) = rstar         (i,m3-1,ipat)
            cstar         (i,m3,ipat) = cstar         (i,m3-1,ipat)
            veg_albedo    (i,m3,ipat) = veg_albedo    (i,m3-1,ipat)
            veg_fracarea  (i,m3,ipat) = veg_fracarea  (i,m3-1,ipat)
            veg_lai       (i,m3,ipat) = veg_lai       (i,m3-1,ipat)
            veg_tai       (i,m3,ipat) = veg_tai       (i,m3-1,ipat)
            veg_rough     (i,m3,ipat) = veg_rough     (i,m3-1,ipat)
            veg_height    (i,m3,ipat) = veg_height    (i,m3-1,ipat)
            patch_area    (i,m3,ipat) = patch_area    (i,m3-1,ipat)
            patch_rough   (i,m3,ipat) = patch_rough   (i,m3-1,ipat)
            patch_wetind  (i,m3,ipat) = patch_wetind  (i,m3-1,ipat)
            leaf_class    (i,m3,ipat) = leaf_class    (i,m3-1,ipat)
            soil_rough    (i,m3,ipat) = soil_rough    (i,m3-1,ipat)
            sfcwater_nlev (i,m3,ipat) = sfcwater_nlev (i,m3-1,ipat)
            stom_resist   (i,m3,ipat) = stom_resist   (i,m3-1,ipat)
            ground_rsat   (i,m3,ipat) = ground_rsat   (i,m3-1,ipat)
            ground_rvap   (i,m3,ipat) = ground_rvap   (i,m3-1,ipat)
            veg_water     (i,m3,ipat) = veg_water     (i,m3-1,ipat)
            veg_hcap      (i,m3,ipat) = veg_hcap      (i,m3-1,ipat)
            veg_energy    (i,m3,ipat) = veg_energy    (i,m3-1,ipat)
            can_prss      (i,m3,ipat) = can_prss      (i,m3-1,ipat)
            can_theta     (i,m3,ipat) = can_theta     (i,m3-1,ipat)
            can_rvap      (i,m3,ipat) = can_rvap      (i,m3-1,ipat)
            can_co2       (i,m3,ipat) = can_co2       (i,m3-1,ipat)
            sensible      (i,m3,ipat) = can_co2       (i,m3-1,ipat)
            evap          (i,m3,ipat) = can_co2       (i,m3-1,ipat)
            transp        (i,m3,ipat) = can_co2       (i,m3-1,ipat)
            gpp           (i,m3,ipat) = gpp           (i,m3-1,ipat)
            plresp        (i,m3,ipat) = plresp        (i,m3-1,ipat)
            resphet       (i,m3,ipat) = resphet       (i,m3-1,ipat)
            veg_ndvip     (i,m3,ipat) = veg_ndvip     (i,m3-1,ipat)
            veg_ndvic     (i,m3,ipat) = veg_ndvic     (i,m3-1,ipat)
            veg_ndvif     (i,m3,ipat) = veg_ndvif     (i,m3-1,ipat)

            do k = 1,mzg
               soil_water       (k,i,1,ipat) = soil_water         (k,i,2,ipat)
               soil_energy      (k,i,1,ipat) = soil_energy        (k,i,2,ipat)
               soil_text        (k,i,1,ipat) = soil_text          (k,i,2,ipat)

               soil_water      (k,i,m3,ipat) = soil_water      (k,i,m3-1,ipat)
               soil_energy     (k,i,m3,ipat) = soil_energy     (k,i,m3-1,ipat)
               soil_text       (k,i,m3,ipat) = soil_text       (k,i,m3-1,ipat)
            end do

            do k = 1,mzs
               sfcwater_mass    (k,i,1,ipat) = sfcwater_mass      (k,i,2,ipat)
               sfcwater_energy  (k,i,1,ipat) = sfcwater_energy    (k,i,2,ipat)
               sfcwater_depth   (k,i,1,ipat) = sfcwater_depth     (k,i,2,ipat)

               sfcwater_mass   (k,i,m3,ipat) = sfcwater_mass   (k,i,m3-1,ipat)
               sfcwater_energy (k,i,m3,ipat) = sfcwater_energy (k,i,m3-1,ipat)
               sfcwater_depth  (k,i,m3,ipat) = sfcwater_depth  (k,i,m3-1,ipat)
            end do
         end do
      end if

   end do
   return
end subroutine leaf_bcond
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine sfc_fields(m1,m2,m3,ia,iz,ja,jz,jd  &
   ,theta,rv,co2p,up,vp,dn0,pp,pi0,rtgt,zt,zm,ths2,rvs2,co2s2,ups2,vps2,pis2,dens2,zts2)

   use leaf_coms
   use rconstants

   implicit none

   integer :: m1,m2,m3,ia,iz,ja,jz,jd
   real, dimension(m1,m2,m3) :: theta,rv,co2p,up,vp,dn0,pp,pi0
   real, dimension(m2,m3) :: rtgt,ths2,rvs2,co2s2,ups2,vps2,pis2,dens2,zts2
   real, dimension(m1) :: zt,zm

   integer :: i,j
   real :: hcpi

   ! Compute surface atmospheric conditions

   hcpi = .5 * cpi

   do j = ja,jz
      do i = ia,iz
         ths2(i,j) = theta(2,i,j)
         rvs2(i,j) = rv(2,i,j)
         co2s2(i,j) = co2p(2,i,j)
         ups2(i,j) = (up(2,i-1,j) + up(2,i,j)) * .5
         vps2(i,j) = (vp(2,i,j-jd) + vp(2,i,j)) * .5
         zts2(i,j) = zt(2) * rtgt(i,j)
         pis2(i,j) = (pp(1,i,j) + pi0(1,i,j) + pp(2,i,j) + pi0(2,i,j)) * hcpi
         dens2(i,j) = (dn0(1,i,j) + dn0(2,i,j)) * .5
      enddo
   enddo

   return
end subroutine sfc_fields

!****************************************************************************

subroutine sfc_fields_adap(m1,m2,m3,ia,iz,ja,jz,jd,flpu,flpv,flpw  &
   ,topma,aru,arv,theta,rv,co2p,up,vp,dn0,pp,pi0,zt,zm,dzt       &
   ,ths2,rvs2,co2s2,ups2,vps2,pis2,dens2,zts2)

   use leaf_coms
   use rconstants

   implicit none

   integer :: m1,m2,m3,ia,iz,ja,jz,jd
   real, dimension(m2,m3) :: flpu,flpv,flpw
   real, dimension(m1,m2,m3) :: aru,arv,theta,rv,co2p,up,vp,dn0,pp,pi0
   real, dimension(m2,m3) :: topma,ths2,rvs2,co2s2,ups2,vps2,pis2,dens2,zts2
   real, dimension(m1) :: zt,zm,dzt

   integer :: i,j,k1,k2,k3
   real :: topma_t,wtw,wtu1,wtu2,wtv1,wtv2

   ! Compute surface atmospheric conditions

   do j = ja,jz
      do i = ia,iz
         k2 = nint(flpw(i,j))
         k1 = k2 - 1
         k3 = k2 + 1

         topma_t = .25 * (topma(i,j) + topma(i-1,j)  &
                 + topma(i,j-jd) + topma(i-1,j-jd))

   ! weights for lowest predicted points, relative to points above them

         wtw = (zm(k2) - topma_t) * dzt(k2)
         wtu1 = aru(nint(flpu(i-1,j)),i-1,j)   / aru(nint(flpu(i-1,j))+1,i-1,j)
         wtu2 = aru(nint(flpu(i,j)),i,j)       / aru(nint(flpu(i,j))+1,i,j)
         wtv1 = arv(nint(flpv(i,j-jd)),i,j-jd) / arv(nint(flpv(i,j-jd))+1,i,j-jd)
         wtv2 = arv(nint(flpv(i,j)),i,j)       / arv(nint(flpv(i,j))+1,i,j)

         ths2(i,j)  =  wtw * theta(k2,i,j) + (1. - wtw)  * theta(k3,i,j)

         rvs2(i,j)  =  wtw * rv(k2,i,j)    + (1. - wtw)  * rv(k3,i,j)
         co2s2(i,j) =  wtw * co2p(k2,i,j)  + (1. - wtw)  * co2p(k3,i,j)

         ups2(i,j) = (wtu1        * up(nint(flpu(i-1,j)),i-1,j)    &
                   +  (1. - wtu1) * up(nint(flpu(i-1,j))+1,i-1,j)  &
                   +  wtu2        * up(nint(flpu(i,j)),i,j)        &
                   +  (1. - wtu2) * up(nint(flpu(i,j))+1,i,j)) * .5

         vps2(i,j) = (wtv1        * vp(nint(flpv(i,j-jd)),i,j-jd)    &
                   +  (1. - wtv1) * vp(nint(flpv(i,j-jd))+1,i,j-jd)  &
                   +  wtv2        * vp(nint(flpv(i,j)),i,j)          &
                   +  (1. - wtv2) * vp(nint(flpv(i,j))+1,i,j)) * .5

         zts2(i,j) = (wtw * (zt(k2) - zm(k1))  &
                   + (1. - wtw) * (zt(k3) - zm(k2)))

         if (wtw >= .5) then
            pis2(i,j)  = ((wtw - .5) * (pp(k1,i,j) + pi0(k1,i,j))  &
                       + (1.5 - wtw) * (pp(k2,i,j) + pi0(k2,i,j))) * cpi
            dens2(i,j) = (wtw - .5)  * dn0(k1,i,j)  &
                       + (1.5 - wtw) * dn0(k2,i,j)
         else
            pis2(i,j)  = ((wtw + .5) * (pp(k2,i,j) + pi0(k2,i,j))  &
                       + (.5 - wtw) * (pp(k3,i,j) + pi0(k3,i,j))) * cpi
            dens2(i,j) = (wtw + .5) * dn0(k2,i,j)  &
                       + (.5 - wtw) * dn0(k3,i,j)
         endif

      enddo
   enddo

   return
end subroutine sfc_fields_adap

!****************************************************************************

subroutine sfc_pcp(nqparm,level,i,j,cuparm,micro)

   use mem_basic
   use mem_grid
   use mem_micro
   use mem_cuparm
   use leaf_coms
   use rconstants, only : cliq,tsupercool,wdnsi

   implicit none

   integer :: nqparm,level,i,j,icld
   type (cuparm_vars)  cuparm
   type (micro_vars)   micro

   if (nqparm > 0) then
      pcpgl = 0.
      do icld=1,nclouds
         pcpgl  = pcpgl + cuparm%conprr(i,j,icld)
      end do
      pcpgl  = pcpgl  * dtll
      qpcpgl = pcpgl  * cliq * (atm_temp - tsupercool)
      dpcpgl = pcpgl  * wdnsi
      pcpgc  = dtlc_factor * pcpgl
      qpcpgc = dtlc_factor * qpcpgl

   else

      pcpgl  = 0.
      qpcpgl = 0.
      dpcpgl = 0.
      pcpgc  = 0.
      qpcpgc = 0.

   endif

   if (level >= 3) then

      pcpgl  = pcpgl  + dtll_factor * micro%pcpg(i,j)
      qpcpgl = qpcpgl + dtll_factor * micro%qpcpg(i,j)
      dpcpgl = dpcpgl + dtll_factor * micro%dpcpg(i,j)
      pcpgc  = pcpgc  + dtlc_factor * dtll_factor * micro%pcpg(i,j)
      qpcpgc = qpcpgc + dtlc_factor * dtll_factor * micro%qpcpg(i,j)

   endif

   return
end subroutine sfc_pcp

!****************************************************************************

subroutine vegndvi(ifm    &
   ,patch_area,leaf_class,veg_fracarea,veg_lai,veg_tai,veg_rough  &
   ,veg_height,veg_albedo,veg_ndvip,veg_ndvic,veg_ndvif)

   use leaf_coms
   use rconstants
   use io_params
   use mem_grid

   !CATT
   use catt_start, only: CATT  ! INTENT(IN)

   implicit none

   integer :: ifm
   real :: patch_area,leaf_class,veg_fracarea,veg_lai,veg_tai,veg_rough  &
      ,veg_height,veg_albedo,veg_ndvip,veg_ndvic,veg_ndvif

   integer :: nveg
   integer, save :: nvcall = 0

   real :: timefac_ndvi,sr,fpar,dead_lai,green_frac
   real, save :: sr_min=1.081,fpar_min=.001,fpar_max=.950,fpcon=-.3338082
   real, save :: ccc=-2.9657
   real, save :: bz=.91,hz=.0075,extinc_veg=.5

   real, dimension(nvtyp+nvtyp_teb), save :: dfpardsr

   !  Initialize dfpardsr array

   if (nvcall == 0) then
      nvcall = 1
      do nveg = 1,(nvtyp+nvtyp_teb)
         dfpardsr(nveg) = (fpar_max - fpar_min) / (sr_max(nveg) - sr_min)
      enddo
   endif

   !  Compute LAI, vegetation roughness, albedo, vegfrac from time-dependent NDVI

   nveg = nint(leaf_class)

   if (tai_max(nveg) < .1) then

      veg_lai = 0.
      veg_tai = 0.
      veg_rough = 0.
      veg_albedo = 0.
      veg_fracarea = 0.

   else

      if (iupdndvi == 0) then
         timefac_ndvi = 0.
      else
         
         ! HOW LARGE COULD THIS GET? SHOULD I MAKE IT A 64bit?  RGK
         timefac_ndvi = sngl((time - ndvitime1(ifm)) / (ndvitime2(ifm) - ndvitime1(ifm)))

      endif

   !  Time-interpolate ndvi to get current value veg_ndvic(i,j) for this patch

      veg_ndvic = veg_ndvip + (veg_ndvif - veg_ndvip) * timefac_ndvi

   !  Limit ndvi to prevent values > .99 to prevent division by zero.

      if (veg_ndvic > .99) veg_ndvic = .99

   ! Compute "simple ratio" and limit between sr_min and sr_max(nveg).

      sr = (1. + veg_ndvic) / (1. - veg_ndvic)

      if (sr < sr_min) then
         sr = sr_min
      elseif (sr > sr_max(nveg)) then
         sr = sr_max(nveg)
      endif

   ! Compute fpar

      fpar = fpar_min + (sr - sr_min) * dfpardsr(nveg)

   ! Compute green leaf area index (veg_lai), dead leaf area index (dead_lai),
   ! total area index (tai), and green fraction

      veg_lai = glai_max(nveg) * (veg_clump(nveg) * fpar / fpar_max  &
              + (1. - veg_clump(nveg)) * alog(1. - fpar) * fpcon)

      dead_lai = (glai_max(nveg) - veg_lai) * dead_frac(nveg)
      veg_tai = veg_lai + sai(nveg) + dead_lai
      green_frac = veg_lai / veg_tai

   ! Compute vegetation roughness height, albedo, and fractional area

      veg_rough = veg_height * (1. - bz * exp(-hz * veg_tai))
      veg_albedo = albv_green(nveg) * green_frac  &
                 + albv_brown(nveg) * (1. - green_frac)
      veg_fracarea = veg_frac(nveg) * (1. - exp(-extinc_veg * veg_tai))

   endif
return
end subroutine vegndvi
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This routine is called by the radiation parameterization and by LEAF.  It computes   !
! net surface albedo plus radiative exchange between the atmosphere, vegetation, and the   !
! snow/ground given previously computed downward longwave and shortwave fluxes from the    !
! atmosphere.  Also computed are functions of snowcover that are required for the above    !
! radiation calculations as well as other calculations in LEAF.                            !
!     The shortwave parameterizations are only valid if the cosine of the zenith angle is  !
! greater than .03 .  Water albedo from Atwater and Bell (1981) alg, als, and alv are the  !
! albedos of the ground, snow, and vegetation (als needs a better formula based on age of  !
! the surface snow).  absg and vctr32 are the actual fractions of shortwave incident on    !
! snow plus ground that get absorbed by the ground and each snow layer, respectively.      !
! They currently use the variable fractrans, which is the fraction of light transmitted    !
! through each layer based on mass per square meter.  algs is the resultant albedo from    !
! snow plus ground.                                                                        !
!------------------------------------------------------------------------------------------!
subroutine sfcrad(mzg,mzs,ip,soil_energy,soil_water,soil_text,sfcwater_energy              &
                 ,sfcwater_mass,sfcwater_depth,patch_area,leaf_class,veg_height            &
                 ,veg_fracarea,veg_albedo,sfcwater_nlev,veg_energy,veg_water,veg_hcap      &
                 ,can_prss,can_theta,can_rvap,rshort,rlong,albedt,rlongup,cosz             &
                 ,g_urban, etown, albtown, tstown)
   use mem_leaf
   use leaf_coms
   use rconstants
   use mem_scratch
   use node_mod     , only : mynum         ! ! intent(in)
   use therm_lib    , only : qwtk          & ! subroutine
                           , qtk           & ! subroutine
                           , ptqz2enthalpy & ! function
                           , idealdenssh   ! ! function
   use catt_start   , only : CATT          ! ! intent(in)
   use teb_spm_start, only : TEB_SPM       ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                , intent(in)    :: mzg,mzs,ip
   real   , dimension(mzg), intent(in)    :: soil_energy,soil_water,soil_text
   real   , dimension(mzs), intent(in)    :: sfcwater_energy,sfcwater_depth ,sfcwater_mass
   real                   , intent(in)    :: patch_area,leaf_class,veg_height,veg_fracarea
   real                   , intent(in)    :: veg_albedo,sfcwater_nlev
   real                   , intent(in)    :: veg_energy,veg_water,veg_hcap
   real                   , intent(in)    :: can_prss,can_theta,can_rvap
   real                   , intent(in)    :: rshort,rlong,cosz
   real                   , intent(in)    :: g_urban, etown, albtown, tstown
   real                   , intent(inout) :: albedt,rlongup
   !----- Local variables. ----------------------------------------------------------------!
   integer                                :: k,m,nsoil,nveg,ksn
   real                                   :: alb,vfc,fcpct,alg,rad,als,fractrans
   real                                   :: absg,algs,emv,emgs,gslong,vlong,alv
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     First we update the canopy air properties.                                        !
   !---------------------------------------------------------------------------------------!
   can_temp     = can_theta * (p00i * can_prss) ** rocp
   can_shv      = can_rvap / (can_rvap + 1.)
   can_enthalpy = ptqz2enthalpy(can_prss,can_temp,can_shv,can_depth)
   can_rhos     = idealdenssh(can_prss,can_temp,can_shv)

   if (ip == 1) then
      !----- Combuting the albedo and upward longwave for water patches. ------------------!
      if (cosz > .03) then
         alb = min(max(-.0139 + .0467 * tan(acos(cosz)),.03),.999)
         albedt = albedt + patch_area * alb
      end if
      call qtk(soil_energy(mzg),tempk(mzg),fracliq(mzg))
      rlongup = rlongup + patch_area * stefan * tempk(mzg) ** 4

   elseif (isfcl == 0) then
      !------ Not running a land surface model, use prescribed value of can_temp. ---------!
      albedt  = albedt  + patch_area * albedo
      rlongup = rlongup + patch_area * stefan * can_temp ** 4
   else
      !------ Running an actual land surface model... -------------------------------------!


      !------ Diagnose vegetation temperature and surface water liquid fraction. ----------!
      call qwtk(veg_energy,veg_water,veg_hcap,veg_temp,veg_fliq)

      !------ Diagnose soil temperature and liquid fraction. ------------------------------!
      do k = 1,mzg
         nsoil = nint(soil_text(k))
         call qwtk(soil_energy(k),soil_water(k)*wdns,slcpd(nsoil),tempk(k),fracliq(k))
      end do

      !------ Diagnose snow temperature and the influence of snow covering veg. -----------!
      nveg = nint(leaf_class)
      ksn  = nint(sfcwater_nlev)
      snowfac = 0.
      do k = 1,ksn
         snowfac = snowfac + sfcwater_depth(k)
         call qtk(sfcwater_energy(k),tempk(k+mzg),fracliq(k+mzg))
      end do
      snowfac = min(.99, snowfac / max(.001,veg_height))

      !------ Defining the exposed area. --------------------------------------------------!
      vf = veg_fracarea * (1. - snowfac)
      vfc = 1. - vf

      !------------------------------------------------------------------------------------!
      !     Shortwave radiation calculations.                                              !
      !------------------------------------------------------------------------------------!
      nsoil=nint(soil_text(mzg))
      fcpct = soil_water(mzg) / slmsts(nsoil)
      if (fcpct > .5) then
         alg = .14
      else
         alg = .31 - .34 * fcpct
      end if
      alv = veg_albedo

      rad = 1.
      if (ksn > 0) then
         !------ als = .14 (the wet soil value) for all-liquid. ---------------------------!
         als = .5 - .36 * fracliq(ksn+mzg)
         rad = 1. - als
      end if
      do k = ksn,1,-1
         fractrans = exp(-20. * sfcwater_depth(k))
         vctr32(k) = rad * (1. - fractrans)
         rad = rad * fractrans
      end do
      absg = (1. - alg) * rad
      algs = 1. - absg
      do k = ksn,1,-1
         algs = algs - vctr32(k)
         rshort_s(k) = rshort * vfc * vctr32(k)
      end do
      rshort_g = rshort * vfc * absg
      rshort_v = rshort * vf * (1. - alv + vfc * algs)
      alb      = vf * alv + vfc * vfc * algs

      !----- Adding urban contribution if running TEB. ------------------------------------!
      if (teb_spm==1) then
         if (nint(g_urban) == 0) then
            albedt = albedt + patch_area * alb
         else
            albedt = albedt + patch_area * albtown
         endif
      else
         albedt = albedt + patch_area * alb
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      !     Longwave radiation calculations.                                               !
      !------------------------------------------------------------------------------------!
      emv  = emisv(nveg)
      emgs = emisg(nsoil)
      if (ksn > 0) emgs = 1.0
      gslong = emgs * stefan * tempk(ksn+mzg) ** 4
      vlong  = emv * stefan * veg_temp ** 4

      rlonga_v  = rlong  * vf * (emv + vfc * (1. - emgs))
      rlonga_gs = rlong  * vfc * emgs
      rlongv_gs = vlong  * vf * emgs
      rlongv_a  = vlong  * vf * (2. - emgs - vf + emgs * vf)
      rlonggs_v = gslong * vf * emv
      rlonggs_a = gslong * vfc
      rlonga_a  = rlong  * (vf * (1. - emv) + vfc * vfc * (1. - emgs))

      !----- Adding urban contribution if running TEB. ------------------------------------!
      if (teb_spm==1) then
         if (nint(g_urban) == 0) then
            rlongup = rlongup + patch_area * (rlongv_a + rlonggs_a + rlonga_a)
         else
            rlongup = rlongup + patch_area * etown * stefan * tstown**4
         endif
      else
         rlongup = rlongup + patch_area * (rlongv_a + rlonggs_a + rlonga_a)
      endif

      !------------------------------------------------------------------------------------!
      !      In case rlong is not computed, zero out all longwave fluxes other than        !
      ! rlongup.  [On the first timestep, radiative fluxes may not be available until      !
      ! microphysics is called, and zeroing these fluxes avoids the imbalance of having    !
      ! upward longwave without downward longwave in LEAF-3.  Also, this allows LEAF-3 to  !
      ! run without radiation for all timesteps, if desired for testing.].                 !
      !------------------------------------------------------------------------------------!
      if (rlong < .1) then
         rlonga_v  = 0.
         rlonga_gs = 0.
         rlongv_gs = 0.
         rlongv_a  = 0.
         rlonggs_v = 0.
         rlonggs_a = 0.
         rlonga_a  = 0.
      end if
      !------------------------------------------------------------------------------------!
   end if

   return
end subroutine sfcrad
!==========================================================================================!
!==========================================================================================!
