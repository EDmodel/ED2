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
   real, dimension(mxp,myp)             :: l_vps2,l_zts2
   !---------------------------------------------------------------------------------------!
   
   if (nstbot == 0) return

   ng=ngrid


   !----- Calling LEAF-3 main driver ------------------------------------------------------!
   call leaf3(mzp,mxp,myp,nzg,nzs,npatch,ia,iz,ja,jz,leaf_g (ng),basic_g (ng),turb_g (ng)  &
             ,radiate_g(ng),grid_g (ng),cuparm_g(ng),micro_g(ng),l_ths2,l_rvs2,l_pis2      &
             ,l_dens2,l_ups2,l_vps2,l_zts2,teb_g(ng),tebc_g(ng))

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
                  , leaf_g(ng)%veg_albedo          , leaf_g(ng)%veg_fracarea               &
                  , leaf_g(ng)%veg_lai             , leaf_g(ng)%veg_tai                    &
                  , leaf_g(ng)%veg_rough           , leaf_g(ng)%veg_height                 &
                  , leaf_g(ng)%patch_area          , leaf_g(ng)%patch_rough                &
                  , leaf_g(ng)%patch_wetind        , leaf_g(ng)%leaf_class                 &
                  , leaf_g(ng)%soil_rough          , leaf_g(ng)%sfcwater_nlev              &
                  , leaf_g(ng)%stom_resist         , leaf_g(ng)%ground_rsat                &
                  , leaf_g(ng)%ground_rvap         , leaf_g(ng)%veg_water                  &
                  , leaf_g(ng)%veg_temp            , leaf_g(ng)%can_rvap                   &
                  , leaf_g(ng)%can_temp            , leaf_g(ng)%veg_ndvip                  &
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
                ,ths2,rvs2,pis2,dens2,ups2,vps2,zts2,teb,tebc)

   use mem_all
   use leaf_coms
   use rconstants
   use node_mod       , only : mynum       ! ! intent(in)
   use catt_start     , only : CATT        ! ! intent(in)
   use teb_spm_start  , only : TEB_SPM     ! ! intent(in)
   use mem_teb        , only : teb_vars    ! ! type
   use mem_teb_common , only : teb_common  ! ! type
   use therm_lib      , only : rslif       & ! function
                             , level       ! ! intent(in)
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
   real, dimension(m2,m3), intent(inout) :: ths2,rvs2,pis2,dens2,ups2,vps2,zts2
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
   ! longer than about 15 seconds.  This allows values of hcapcan = 2.e4, wcapcan = 2.e1,  !
   ! and hcapveg = 3.e4 as are now defined below.                                          !
   !---------------------------------------------------------------------------------------!
   niter_leaf  = max(1,nint(dtlt/40.+.4))
   niter_can   = max(1,nint(dtll/15.+.4))
   dtll_factor = 1. / float(niter_leaf)
   dtll        = dtlt * dtll_factor
   dtlc_factor = 1. / float(niter_can)
   dtlc        = dtll * dtlc_factor

   !----- Defining some variables to speed up calculations --------------------------------!
   dtllohcc = dtll / hcapcan
   dtllowcc = dtll / wcapcan
   dtlcohcc = dtlc / hcapcan
   dtlcowcc = dtlc / wcapcan
   dtlcohcv = dtlc / hcapveg

   !----- Copy surface atmospheric variables into 2d arrays for input to LEAF. ------------!

   if (if_adap == 1) then
      call sfc_fields_adap(m1,m2,m3,ia,iz,ja,jz,jdim      , grid%flpu     , grid%flpv      &
                          , grid%flpw     , grid%topma    , grid%aru      , grid%arv       &
                          , basic%theta   , basic%rv      , basic%up      , basic%vp       &
                          , basic%dn0     , basic%pp      , basic%pi0                      &
                          , zt,zm,dzt,ths2,rvs2,ups2,vps2,pis2,dens2,zts2)
   else
      call sfc_fields(m1,m2,m3,ia,iz,ja,jz,jdim           , basic%theta   , basic%rv       &
                          , basic%up      , basic%vp      , basic%dn0     ,  basic%pp      &
                          , basic%pi0     , grid%rtgt                                      &
                          , zt,ths2,rvs2,ups2,vps2,pis2,dens2,zts2)
   end if

   !----- Big domain loop -----------------------------------------------------------------!
   jloop1: do j = ja,jz
      iloop1: do i = ia,iz

         !----- Copy surface variables to single-column values ----------------------------!
         ups  = ups2(i,j)
         vps  = vps2(i,j)
         ths  = ths2(i,j)
         rvs  = rvs2(i,j)
         zts  = zts2(i,j)
         pis  = pis2(i,j)
         dens = dens2(i,j)

         prss     = pis ** cpor * p00
         vels     = sqrt(ups ** 2 + vps ** 2)
         gzotheta = g * zts / ths

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

         !----- For no soil model (patch 2) fill "canopy" temperature and moisture. -------!
         if (isfcl == 0) then
            leaf%can_temp  (i,j,2) = (ths - dthcon) * pis
            leaf%can_rvap  (i,j,2) = rvs - drtcon
            leaf%patch_area(i,j,1) = 1. - pctlcon
            leaf%patch_area(i,j,2) = pctlcon
         end if

         !----- Begin patch loop. ---------------------------------------------------------!
         ploop1: do ip = 1,np

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
                       , leaf%can_temp        (  i,j,ip) , leaf%veg_temp        (  i,j,ip) &
                       , leaf%leaf_class      (  i,j,ip) , leaf%veg_height      (  i,j,ip) &
                       , leaf%veg_fracarea    (  i,j,ip) , leaf%veg_albedo      (  i,j,ip) &
                       , leaf%sfcwater_nlev   (  i,j,ip) , radiate%rshort       (  i,j   ) &
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

                  leaf%ground_rsat(i,j,ip) = rslif(prss,tempk(mzg))
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
               !     Calculate turbulent fluxes between atmosphere and canopy (or          !
               ! "canopy").                                                                !
               !---------------------------------------------------------------------------!
               if (leaf%patch_area(i,j,ip) >= .009 .or. isfcl == 0) then
                  thetacan = leaf%can_temp(i,j,ip) / pis
               end if

               call stars( leaf%ustar          (i,j,ip) , leaf%tstar              (i,j,ip) &
                         , leaf%rstar          (i,j,ip) , ths                              &
                         , rvs                          , thetacan                         &
                         , leaf%can_rvap       (i,j,ip) , zts                              &
                         , leaf%patch_rough    (i,j,ip) , leaf%patch_area         (i,j,ip) &
                         , vels                         , vels_pat                         &
                         , vonk                         , dtllohcc                         &
                         , dens                         , dtll                             &
                         , leaf%R_aer          (i,j,ip) )

               if (teb_spm==1) g_urban = leaf%g_urban(i,j,ip)

               call sfclmcv( leaf%ustar        (i,j,ip) , leaf%tstar              (i,j,ip) &
                           , leaf%rstar        (i,j,ip) , vels                             &
                           , vels_pat                   , ups                              &
                           , vps                        , gzotheta                         &
                           , leaf%patch_area   (i,j,ip) , turb%sflux_u            (i,j   ) &
                           , turb%sflux_v      (i,j   ) , turb%sflux_w            (i,j   ) &
                           , turb%sflux_t      (i,j   ) , turb%sflux_r            (i,j   ) &
                           , g_urban                    )

               !---------------------------------------------------------------------------!
               !     For water patches, update temperature and moisture of "canopy" from   !
               ! divergence of fluxes with water surface and atmosphere.  rdi = ustar/5    !
               ! is the viscous sublayer conductivity from Garratt (1992).                 !
               !---------------------------------------------------------------------------!
               if (ip == 1) then

                  rdi = .2 * leaf%ustar(i,j,1)

                  leaf%can_temp(i,j,1) = leaf%can_temp(i,j,1)                              &
                                       + dtllohcc * dens * cp * rdi                        &
                                       * ((tempk(mzg) - leaf%can_temp(i,j,1))              &
                                       + leaf%ustar(i,j,1) * leaf%tstar(i,j,1) * pis)

                  leaf%can_rvap(i,j,1) = leaf%can_rvap(i,j,1)                              &
                                       + dtllowcc * dens * rdi                             &
                                       * ((leaf%ground_rsat(i,j,1) - leaf%can_rvap(i,j,1)) &
                                       + leaf%ustar(i,j,1) * leaf%rstar(i,j,1))

                  !----- KML - drydep -----------------------------------------------------!
                  if (catt == 1) leaf%R_aer(i,j,ip) = rdi
               end if

               if (teb_spm == 1) then

                  !----- TEB: Urban canopy parameterization starts here -------------------!
                  if (ip >= 2) then

                     if (nint(leaf%g_urban(i,j,ip))/=0.) then

                        !edmilson
                        psup1  = ((basic%pi0(1,i,j)+basic%pp(1,i,j))*cpi)**cpor*p00
                        psup2  = ((basic%pi0(2,i,j)+basic%pp(2,i,j))*cpi)**cpor*p00
                        depe   = psup2-psup1
                        alt2   = zt(1)*grid%rtgt(i,j)
                        deze   = zts-alt2
                        dpdz   = depe/deze
                        exn1st = (psup2/p00)**rocp
                        
                        airt= basic%theta(2,i,j)*exn1st 

                        ! TEB - defining pointers
                        g_urban   = leaf%g_urban(i,j,ip)

                        call leaf3_teb_interface(istp,dtlt,dtll                            &
                           , radiate%cosz(i,j)          , zts                              &
                           , radiate%rlong(i,j)         , radiate%rshort(i,j)              &
                           , psup2                      , airt                             &
                           , ups                        , vps                              &
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
                           , dens                       )
                        
                        turb%sflux_u(i,j) = turb%sflux_u(i,j)                              &
                                          + leaf%patch_area(i,j,ip) * zsfu_town
                        turb%sflux_v(i,j) = turb%sflux_v(i,j)                              &
                                          + leaf%patch_area(i,j,ip) * zsfv_town
                        turb%sflux_t(i,j) = turb%sflux_t(i,j)                              &
                                          + leaf%patch_area(i,j,ip) * zh_town/(cp*dens)
                        turb%sflux_r(i,j) = turb%sflux_r(i,j)                              &
                                          + leaf%patch_area(i,j,ip) * zle_town/(alvl*dens)
                     end if
                  end if
               end if

               !---------------------------------------------------------------------------!
               !     For soil model patches, update temperature and moisture of soil,      !
               ! vegetation, and canopy                                                    !
               !---------------------------------------------------------------------------!
               if (isfcl >= 1 .and. ip >= 2) then
                  if (leaf%patch_area(i,j,ip) >= .009) then

                     call leaftw(mzg,mzs,np                                                &
                        , leaf%soil_water     (:,i,j,ip) , leaf%soil_energy    (:,i,j,ip)  &
                        , leaf%soil_text      (:,i,j,ip) , leaf%sfcwater_mass  (:,i,j,ip)  &
                        , leaf%sfcwater_depth (:,i,j,ip) , leaf%ustar            (i,j,ip)  &
                        , leaf%tstar            (i,j,ip) , leaf%rstar            (i,j,ip)  &
                        , leaf%veg_albedo       (i,j,ip) , leaf%veg_fracarea     (i,j,ip)  &
                        , leaf%veg_lai          (i,j,ip) , leaf%veg_tai          (i,j,ip)  &
                        , leaf%veg_rough        (i,j,ip) , leaf%veg_height       (i,j,ip)  &
                        , leaf%patch_area       (i,j,ip) , leaf%patch_rough      (i,j,ip)  &
                        , leaf%patch_wetind     (i,j,ip) , leaf%leaf_class       (i,j,ip)  &
                        , leaf%soil_rough       (i,j,ip) , leaf%sfcwater_nlev    (i,j,ip)  &
                        , leaf%stom_resist      (i,j,ip) , leaf%ground_rsat      (i,j,ip)  &
                        , leaf%ground_rvap      (i,j,ip) , leaf%veg_water        (i,j,ip)  &
                        , leaf%veg_temp         (i,j,ip) , leaf%can_rvap         (i,j,ip)  &
                        , leaf%can_temp         (i,j,ip) , leaf%veg_ndvip        (i,j,ip)  &
                        , leaf%veg_ndvic        (i,j,ip) , leaf%veg_ndvif        (i,j,ip)  &
                        , radiate%rshort        (i,j)    , radiate%cosz          (i,j)     &
                        , ip,i,j)

                 end if
               end if

               !---------------------------------------------------------------------------!
               !    Converting the surface water internal energy back to an intensive      !
               ! variable.                                                                 !
               !---------------------------------------------------------------------------!
               ksn = nint(leaf%sfcwater_nlev(i,j,ip))
               leaf%sfcwater_energy(:,i,j,ip) = 0.
               do k=1,ksn
                  leaf%sfcwater_energy(k,i,j,ip) = sfcw_energy_ext(k)                      &
                                                 / leaf%sfcwater_mass(k,i,j,ip)
               end do
               !---------------------------------------------------------------------------!

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
subroutine leaftw(mzg,mzs,np,soil_water, soil_energy, soil_text,sfcwater_mass              &
                 ,sfcwater_depth,ustar, tstar,rstar,veg_albedo,veg_fracarea,veg_lai        &
                 ,veg_tai,veg_rough,veg_height,patch_area,patch_rough,patch_wetind         &
                 ,leaf_class,soil_rough,sfcwater_nlev,stom_resist,ground_rsat,ground_rvap  &
                 ,veg_water,veg_temp,can_rvap,can_temp,veg_ndvip,veg_ndvic,veg_ndvif       &
                 ,rshort,cosz,ip,i,j)

   use leaf_coms
   use mem_leaf
   use rconstants
   use mem_scratch
   use therm_lib   , only : qwtk & ! subroutine
                          , qtk  ! ! subroutine
   use catt_start  , only : CATT ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                , intent(in)    :: mzg,mzs,np
   real   , dimension(mzg), intent(inout) :: soil_water,soil_energy,soil_text
   real   , dimension(mzs), intent(inout) :: sfcwater_mass,sfcwater_depth

   real                 :: ustar        ,tstar         ,rstar        &
                          ,veg_albedo   ,veg_fracarea  ,veg_lai      &
                          ,veg_tai                                   &
                          ,veg_rough    ,veg_height    ,patch_area   &
                          ,patch_rough  ,patch_wetind  ,leaf_class   &
                          ,soil_rough   ,sfcwater_nlev ,stom_resist  &
                          ,ground_rsat  ,ground_rvap   ,veg_water    &
                          ,veg_temp     ,can_rvap      ,can_temp     &
                          ,veg_ndvip    ,veg_ndvic     ,veg_ndvif    &
                          ,rshort       ,cosz

   integer :: ip,k,nveg,ksn,nsoil,ksnnew,newlayers,nlayers,kold,ktrans,nsl,kzs
   integer, save :: ncall=0

   integer :: i,j
   real :: sfcw_energy_int

   real :: stretch,thik,pf,snden,vegfracc,qwfree,wfree          &
      ,depthgain,totsnow,qw,w,qwt,wt,soilhcap,fac,wfreeb,depthloss,soilcap    &
      ,sndenmax,sndenmin,snowmin,hold,wtnew,wtold,wdiff,watermid,availwat,wg  &
      ,wloss,soilcond,waterfrac

   real, save, dimension(20) :: thicknet
   real, save, dimension(20,20) :: thick

   do k = 1,mzg
      dslz   (k) = slz(k+1) - slz(k)
      dslzi  (k) = 1. / dslz(k)
      dslzidt(k) = dslzi(k) * dtll
      slzt   (k) = .5 * (slz(k) + slz(k+1))
   enddo

   do k = 2,mzg
      dslzt   (k) = slzt(k) - slzt(k-1)
      dslzti  (k) = 1. / dslzt(k)
      dslztidt(k) = dslzti(k) * dtll
   enddo

   ! Initialize snow thickness scaling array

   if (ncall /= 40) then
      ncall = 40
      stretch = 2.0
      do kzs = 1,mzs
         thik = 1.0
         thicknet(kzs) = 0.0
         do k = 1,(kzs+1)/2
            thick(k,kzs) = thik
            thick(kzs+1-k,kzs) = thik
            thicknet(kzs) = thicknet(kzs) + 2. * thik
            thik = thik * stretch
         enddo
         if ((kzs+1)/2 /= kzs/2) thicknet(kzs) = thicknet(kzs) - thik/stretch
         do k = 1,kzs
            thick(k,kzs) = thick(k,kzs) / thicknet(kzs)
         enddo
      enddo
   endif

   nveg = nint(leaf_class)
   ksn = nint(sfcwater_nlev)

   ! Evaluate any exchanges of heat and moisture to or from vegetation, apply
   ! moisture and heat changes to vegetation, and evaluate the resistance
   ! parameter rd between canopy air and the top soil or snow surface.

   call canopy(mzg,mzs,ksn,nveg  &
      ,soil_energy,soil_water,soil_text,sfcwater_mass  &
      ,ustar,tstar,rstar,soil_rough,veg_rough,veg_height  &
      ,veg_lai,veg_tai,veg_water,veg_temp,leaf_class,veg_fracarea  &
      ,stom_resist,can_temp,can_rvap,ground_rsat,ground_rvap,rshort  &
      ,i,j,ip)

   ! Compute soil and effective snow heat resistance times layer depth (rfactor).

   do k = 1,mzg
      nsoil = nint(soil_text(k))
      waterfrac = soil_water(k) / slmsts(nsoil)
      soilcond = soilcond0(nsoil)  &
         + waterfrac * (soilcond1(nsoil) + waterfrac * soilcond2(nsoil))
      rfactor(k) = dslz(k) / soilcond
   enddo

   do k = 1,ksn
      snden = sfcwater_mass(k) / sfcwater_depth(k)
      rfactor(k+mzg) = sfcwater_depth(k)  &
         / (1.093e-3 * exp(.028 * tempk(k+mzg)) * (.030 + snden  &
         * (.303e-3 + snden * (-.177e-6 + snden * 2.25e-9))))
   enddo

   ! Find soil and snow internal sensible heat fluxes [W/m2]

   hfluxgsc(1) = 0.

   do k = 2,mzg+ksn
      hfluxgsc(k) = - (tempk(k) - tempk(k-1)) / ((rfactor(k) + rfactor(k-1)) * .5)
   enddo

   ! Heat flux at soil or snow top from longwave, sensible, and
   ! upward latent heat fluxes [W/m^2]

   hfluxgsc(ksn+1+mzg) = hflxgc + wflxgc * alvl - rlonga_gs  &
      - rlongv_gs + rlonggs_v + rlonggs_a

   ! Update soil Q values [J/m3] and snow Q values [J/kg] from sensible heat,
   ! upward water vapor (latent heat), longwave, and shortwave fluxes.
   ! This excludes effects of dew/frost formation, precipitation, shedding,
   ! and percolation.  Update top soil or snow moisture from evaporation only.

   do k = 1,mzg
      soil_energy(k) = soil_energy(k) + dslzidt(k) * (hfluxgsc(k) - hfluxgsc(k+1))
   enddo

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

   ! New moisture, qw, and depth from dew/frost formation, precipitation,
   ! shedding, and percolation.  ksnnew is the layer that receives the new
   ! condensate that comes directly from the air above.  If there is no
   ! pre-existing snowcover, this is a temporary "snow" layer.

   if (pcpgl + wshed + dewgnd > min_sfcwater_mass) then
      ksnnew = max(ksn,1)
      vegfracc = 1. - veg_fracarea
      qwfree    = qpcpgl * vegfracc + (dewgnd * alvl + qwshed )* veg_fracarea
      wfree     =  pcpgl * vegfracc + (dewgnd + wshed) * veg_fracarea
      depthgain = dpcpgl * vegfracc + (dewgnd + wshed) * veg_fracarea * wdnsi
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

   wloss = min(transp * dslzidt(ktrans) * wdnsi, wg, soil_water(ktrans) - soilcp(nsl))
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
               ,sfcwater_nlev,ground_rsat,ground_rvap,can_temp,can_rvap,prss)

   return
end subroutine leaftw
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine canopy(mzg,mzs,ksn,nveg  &
   ,soil_energy,soil_water,soil_text,sfcwater_mass  &
   ,ustar,tstar,rstar,soil_rough,veg_rough,veg_height &
   ,veg_lai,veg_tai,veg_water,veg_temp,leaf_class,veg_fracarea  &
   ,stom_resist,can_temp,can_rvap,ground_rsat,ground_rvap,rshort  &
   ,i,j,ip)

   use leaf_coms
   use rconstants

   use therm_lib, only: eslif,rslif,qwtk

   !CATT
   use catt_start, only: CATT

   implicit none

   integer :: mzg,mzs,ksn,nveg

   real, dimension(mzg) :: soil_energy,soil_water,soil_text
   real, dimension(mzs) :: sfcwater_mass
   real :: ustar,tstar,rstar,soil_rough,veg_rough,veg_height             &
          ,veg_lai,veg_tai,veg_water,veg_temp,leaf_class,veg_fracarea    &
          ,stom_resist,can_temp,can_rvap,ground_rsat,ground_rvap,rshort

   integer :: k,kk,nsoil,iter_can

   integer :: i,j,ip

   real :: aux,brad,bswp,bthi,btlo,bvpd,c2,c3,c4,dsm,es,fac,factv         &
      ,fracliqv,fthi,ftlo,frad,fswp,fvpd,qwtot,rasgnd,rasveg,rleaf,rsatveg     &
      ,sigmaw ,slai,stai,slpotv,srad,sswp,sthi,stlo,svpd,swp,tveg,tvegk         &
      ,vpd,wtemp,wtroot,x,zognd,zoveg,zdisp,zveg,wflx,dewgndflx,ustar0         &
      ,transp_test,rc,rc_inf
   ! From RAMS 6.0
   real :: wshed0, transp0

   !     Note: c1=261.5*sqrt((1.-exp(-2.*exar))/(2.*exar))
   !     from Lee's dissertation, Eq. 3.36.  The factor of 261.5 is
   !     100 * ln((h-d)/zo) / vonk   where d = .63 * h and zo = .13 * h.
   !     The factor of 100 is 1/L in Eq. 3.37.  Thus, c1 * ustar is the
   !     total expression inside the radical in Eq. 3.37.
   !bob      parameter(exar=3.5,covr=2.16,c1=98.8)
   real, parameter :: exar=2.5,covr=2.16,c1=116.6

   data   brad,  srad/    196.,   0.047/,  &
          btlo,  stlo/   281.5,    0.26/,  &
          bthi,  sthi/   310.1,  -0.124/,  &
          bvpd,  svpd/  4850.0, -0.0051/,  &
          bswp,  sswp/ -1.07e6, 7.42e-6/
   !     save

   ! Compute ground-canopy resistance rd.  Assume zognd not affected by snow.
   ! Assume (zoveg,zdisp) decrease linearly with snow depth, attaining
   ! the values (zognd,0) when veg covered.

   zognd = soil_rough
   zoveg = veg_rough * (1.-snowfac) + zognd * snowfac
   zdisp = veg_height * (1.-snowfac)
   zveg = zdisp / 0.63
   !bob      rasgnd = log(zts / zognd) * log((zdisp + zoveg) / zognd)
   !bob     +      / (vonk * vonk * vels)

   ! The following value of ground-canopy resistance for the
   ! nonvegetated (bare soil or water) surface is from John Garratt.
   ! It is 5/ustar and replaces the one from old leaf.

   rasgnd = 5. / ustar

   if (veg_tai >= .1 .and. snowfac < .9) then

   ! If vegetation is sufficiently abundant and not covered by snow, compute
   ! heat and moisture fluxes from vegetation to canopy, and flux resistance
   ! from soil or snow to canopy.

      factv = log(zts / zoveg) / (vonk * vonk * vels)
      aux = exp(exar * (1. - (zdisp + zoveg) / zveg))
      rasveg = factv * zveg / (exar * (zveg - zdisp)) * (exp(exar) - aux)
      c2 = max(0.,min(1., 1.1 * veg_tai / covr))
      rd = rasgnd * (1. - c2) + rasveg * c2
      ! From RAMS 6.0
      wshed = 0.
      qwshed = 0.
      transp = 0.

   else

   ! If the TAI is very small or if snow mostly covers the vegetation, bypass
   ! vegetation computations.  Set heat and moisture flux resistance rd between
   ! the "canopy" and snow or soil surface to its bare soil value.  Set shed
   ! precipitation heat and moisture to unintercepted values.

      wshed = pcpgl
      qwshed = qpcpgl
      transp = 0.
      rd = rasgnd

   endif

   ! Compute sensible heat and moisture fluxes between top soil or snow surface
   ! and canopy air.  wflxgc [kg/m2/s] is the upward vapor flux from soil or snow
   ! evaporation and dewgnd is the mass of dew that forms on the snow/soil
   ! surface this timestep; both are defined as always positive or zero.

   hflxgc = cp * dens * (tempk(mzg+ksn) - can_temp) / rd
   wflx = dens * (ground_rsat - can_rvap) / rd
   dewgndflx = max(0.,-wflx)
   dewgnd = dewgndflx * dtll
   if (ksn == 0) then
      wflxgc = max(0.,dens * (ground_rvap - can_rvap) / rd)
   else
      wflxgc = max(0.,min(wflx,sfcwater_mass(ksn)/dtll))
   endif

   if (veg_tai >= .1 .and. snowfac < .9) then

   ! If vegetation is sufficiently abundant and not covered by snow, compute
   ! heat and moisture fluxes from vegetation to canopy, and flux resistance
   ! from soil or snow to canopy.

   ! Compute veg-canopy resistance rb.  Note that rb and rc are both defined
   ! WITHOUT the LAI factor; the LAI factor is included later in computing
   ! fluxes that involve rb and/or rc.

      ustar0 = max(.1,ustar)
      c4 = .01 * sqrt(ustar0 * c1)

      rb  = (1. + .5 * veg_tai) / (.01 * sqrt(ustar0 * c1))

   ! Soil water potential factor for stomatal control

      swp = -200.0
      nveg = nint(leaf_class)

      do k = kroot(nveg),mzg
         nsoil = nint(soil_text(k))
         slpotv = slpots(nsoil) * (slmsts(nsoil) / soil_water(k)) ** slbs(nsoil)
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !  Bypass this computation, and replace with wtroot = 1 for root zone
   !      wtroot = float(int(root(k) + .95))
   !      wtroot = 1.
   !      swp = swp + wtroot * (max(swp,slpotv) - swp)
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (slpotv > swp) swp = slpotv
      enddo
      swp = swp * g * wdns

   ! Begin canopy time-split iterations

      do iter_can = 1,niter_can

   ! Calculate the saturation vapor pressure at TVEG

         tveg = veg_temp
         es = eslif(tveg)
         rsatveg = rslif(prss,tveg)
         
   ! Compute mixing ratio at leaf surface using previous rc

         rc = stom_resist
         rleaf = (rb * rsatveg + rc * can_rvap) / (rb + rc)
         vpd = max((es - rleaf * prss / (ep + rleaf)),0.)

   ! Evaluate 5 environmental factors and new rc

         ftlo = 1. + exp(-stlo * (tveg - btlo))
         fthi = 1. + exp(-sthi * (tveg - bthi))
         frad = 1. + exp(-srad * (rshort - brad))
         fswp = 1. + exp(-sswp * (swp - bswp))
         fvpd = 1. + exp(-svpd * (vpd - bvpd))

   ! 15-minute response time for stomatal conductance (must have dtlc <= 900.)

         rc_inf = ftlo * fthi * frad * fvpd * fswp * rcmin(nveg)
         rc = 1. / (1. / rc + .0011 * dtlc * (1. / rc_inf - 1. / rc))

   ! Limit maximum transpiration to be <= 400 w m-2 by increasing rc if necessary.

         transp_test = alvl * dens * veg_lai * (rsatveg - can_rvap) / (rb + rc)
         if (transp_test > 400.) then
            rc = (rb + rc) * transp_test * .0025 - rb
         endif

         stom_resist = rc

   ! Flux of heat and moisture from vegetation to canopy

         stai = veg_tai * (1. - snowfac)
         slai = veg_lai * (1. - snowfac)

         hflxvc = 2. * stai * cp * dens * (tveg - can_temp) / rb
         c3 = dens * (rsatveg - can_rvap)
         if (c3 >= 0.) then
            sigmaw = min(1.,(veg_water / (.2 * stai)) ** .66667)
            wflxvc = min(c3 * 2. * stai * sigmaw / rb,veg_water/dtlc)
            ! From RAMS 6.0
            !transp = max(0.,c3 * slai * (1.-sigmaw) / (rb + rc))
            transp0 = max(0.,c3 * slai * (1.-sigmaw) / (rb + rc))
         else
            wflxvc = c3 * 2. * stai / rb
            ! From RAMS 6.0
            !transp = 0.
            transp0 = 0.
         endif

   ! Update vegetation moisture from precipitation, vapor flux with canopy,
   ! and shedding of excess moisture.

         wtemp = veg_water - wflxvc * dtlc + pcpgc * vf
         ! From RAMS 6.0
         !wshed = max(0.,wtemp - 0.2 * stai)
         !veg_water = max(0.,wtemp - wshed)
         wshed0 = max(0.,wtemp - 0.2 * stai)
         wshed = wshed + wshed0
         transp = transp + transp0
         veg_water = max(0.,wtemp - wshed0)

   ! Update vegetation temperature from radiative fluxes and sensible and
   ! latent heat transfer with canopy air.  Heat capacity of the vegetation
   ! is defined as hcapveg; it requires further testing and needs to be
   ! increased for dense vegetation.  A value of hcapveg of 3.e4 corresponds
   ! to the heat capacity of a layer of water 7.5 mm thick.

         ! From RAMS 6.0
         tvegk = veg_temp  + dtlcohcv  &
            * (rshort_v + rlonga_v + rlonggs_v - rlongv_gs  &
            -  rlongv_a - hflxvc - (wflxvc + transp0) * alvl)

   ! Exchange heat between vegetation and precipitation in implicit scheme

         qwtot = qpcpgc * vf + hcapveg * tvegk
         call qwtk(qwtot,pcpgc * vf,hcapveg,tvegk,fracliqv)
         veg_temp = tvegk
         ! From RAMS 6.0
         qwshed = qwshed + wshed0 * ( fracliqv * cliq * (tvegk - tsupercool) &
                                    + (1.-fracliqv) * cice * tvegk )

   ! Update temperature and moisture of canopy.  hcapcan [J/m2/K] and
   ! wcapcan [kg_air/m2] are
   ! the heat and moisture capacities of the canopy, while c2 and c3 are deltat
   ! times their inverses, respectively.  They require further testing.

         can_temp = can_temp + dtlcohcc * (hflxgc + hflxvc  &
                  + dens * cp * ustar * tstar * pis)

         ! From RAMS 6.0
         !can_rvap = can_rvap + dtlcowcc * (wflxgc - dewgndflx + wflxvc  &
         !         + transp + dens * ustar * rstar)
         can_rvap = can_rvap + dtlcowcc * (wflxgc - dewgndflx + wflxvc  &
              + transp0 + dens * ustar * rstar)

      enddo

   else

      can_temp = can_temp + dtllohcc * (hflxgc + dens * cp * ustar * tstar * pis)
      can_rvap = can_rvap + dtllowcc * (wflxgc - dewgndflx + dens * ustar * rstar)

   endif

   return
end subroutine canopy

!*****************************************************************************

!subroutine stars(t,ng,mynum,i0,j0,i,j,ip  &
subroutine stars(ustar,tstar,rstar,ths,rvs,thetacan,can_rvap,zts,patch_rough  &
!kml drydep
   ,patch_area,vels,vels_pat,vonk,dtllohcc,dens,dtll                        &
   ,R_aer)
!kml drydep

   !CATT
   use rconstants, only: g
   use catt_start, only: CATT

   implicit none

   real :: ustar,tstar,rstar,ths,rvs,thetacan,can_rvap,zts,patch_rough  &
          ,patch_area,vels,vels_pat,vonk,dtllohcc,dens,dtll,t

   !kml drydep
   real :: R_aer
   !kml drydep

   !integer i,j,ip,ng,mynum,i0,j0

   real :: b,csm,csh,d,a2,c1,ri,fm,fh,c2,cm,ch,c3
   !real :: tstaro,d_theta,d_thetan,theta_new
   real :: ustaro,d_vel,d_veln,vel_new,delz
   integer :: ifix,ifixu

   real, parameter :: ubmin = .25   & ! should use ubmin=1.0 for convec case
                     ,ustmin = .1     !            ubmin=0.1 for stable case

   ! Routine to compute Louis (1981) surface layer parameterization.

   vels_pat = max(vels,ubmin)
   b = 5.
   csm = 7.5
   csh = 5.
   d = 5.

   ! a2 is the drag coefficient in neutral conditions, here same for h/m
   ! ri is the bulk richardson numer, eq. 3.45 in Garratt

   a2 = (vonk / log(zts / patch_rough)) ** 2
   c1 = a2 * vels_pat
   ri = g * zts * (ths - thetacan)  &
      / (.5 * (ths + thetacan) * vels_pat * vels_pat)

   if (ths - thetacan > 0.) then   ! STABLE CASE

      fm = 1. / (1. + (2 * b * ri / sqrt(1 + d * ri)))
      fh = 1. / (1. + (3 * b * ri * sqrt(1 + d * ri)))

   else                            ! UNSTABLE CASE

      c2 = b * a2 * sqrt(zts / patch_rough * (abs(ri)))
      cm = csm * c2
      ch = csh * c2
      fm = (1. - 2 * b * ri / (1. + 2 * cm))
      fh = (1. - 3 * b * ri / (1. + 3 * ch))

   endif

   ustar = max(ustmin,sqrt(c1 * vels_pat * fm))
   c3 = c1 * fh / ustar
   rstar = c3 * (rvs - can_rvap)
   tstar = c3 * (ths - thetacan)

   if (CATT==1) then
      !kml drydep
      R_aer = 1./(a2*vels_pat*fh)
      !kml drydep
   endif

   ! Limit ustar so that the flux cannot take more than 1/2 velocity in a timestep

   ifixu=0
   ustaro=ustar
   delz = 2.*zts
   d_vel =  - ustar * ustar *dtll / delz
   vel_new = vels_pat + d_vel
   if(vel_new < .5 * vels_pat) then
         ifixu=1
         d_veln = .5 * vels_pat
         ustar=sqrt(d_veln*delz/dtll)
   endif

   ! Limit tstar and rstar so that the direction of the gradients cannot change
   !   sign due to the fluxes
      !   can_thetan = thetacan + dtlcohcc * dens * cp * ustar * tstar * pis
   !ifix=0
   !tstaro=tstar
   !delz = 2.*zts
   !d_theta = dtllohcc * dens * cp * ustar * tstar
   !theta_new = thetacan + d_theta
   !if(ths - thetacan > 0.) then
   !   if(theta_new > ths) then
   !      ifix=1
   !      d_thetan = ths - thetacan
   !      tstar=d_thetan/(dtllohcc * dens * cp * ustar)
   !   endif
   !else
   !   if(theta_new < ths) then
   !      ifix=2
   !      d_thetan = ths - thetacan
   !      tstar=d_thetan/(dtllohcc * dens * cp * ustar)
   !   endif
   !endif
   !if( (ng == 1 .and. i0+i == 53 .and. j0+j == 32) .or. ifix /= 0 .or. ifixu /= 0) then
   !if( (ng == 1 .and. i0+i == 53 .and. j0+j == 32) ) then
   !write(11,'(a7,f8.0,6i3,2f8.0,9f10.3)')  &
   !    'leaf--f',t,ng,i0+i,j0+j,ip,ifixu,ifix,-1004.*1.15*ustar*tstar,-2.5e6*1.15*ustar*rstar  &
   !   ,ustar,ustaro,d_vel,d_veln,vels_pat,vel_new,patch_area,patch_rough
      !,ustar,tstar,tstaro,ths,thetacan,theta_new,d_theta,d_thetan,vels_pat,patch_area,patch_rough
   !endif

   !if(ng == 1 .and. i0+i == 53 .and. j0+j == 32)  &
   !if(abs(-1004.*1.15*ustar*tstar) > 1000. .or. abs(-2.5e6*1.15*ustar*rstar) > 1000.)  &
   !write(11,'(a7,f8.0,4i3,2f8.0,2f8.3,f10.5,f9.4,2f7.1,2f9.4,2f7.2,2f8.4,f8.3)')  &
   !    'leaf--3',t,ng,i0+i,j0+j,ip,-1004.*1.15*ustar*tstar,-2.5e6*1.15*ustar*rstar  &
   !   ,ustar,tstar,rstar,ri,ths,thetacan,rvs,can_rvap,vels_pat,patch_rough,c1,fh,patch_area

   return
end subroutine stars

!*****************************************************************************

subroutine sfclmcv(ustar,tstar,rstar,vels,vels_pat,ups,vps,gzotheta,patch_area&
     ,sflux_u,sflux_v,sflux_w,sflux_t,sflux_r                                 &
     ! For TEB
     ,G_URBAN)
     !
!  +---------------------------------------------------------------+
!  \  This routine computes the turbulent fluxes of momentum,      \
!  \  heat and moisture from the surface layer using the           \
!  \  Manton-Cotton algebraic surface layer equations.             \
!  +---------------------------------------------------------------+

   use rconstants

   ! TEB_SPM
   use teb_spm_start, only: TEB_SPM !INTENT(IN)

   implicit none

   ! Arguments:
   real, intent(in)    :: ustar, tstar, rstar, vels,vels_pat, ups, vps, gzotheta,&
        patch_area
   real, intent(inout) :: sflux_u, sflux_v, sflux_w, sflux_t, sflux_r
   ! For TEB
   real:: G_URBAN

   ! Local Variables:
   real :: zoverl
   real :: wtol,cosine1,sine1,vtscr,cx,psin



   data wtol/1e-20/

   cosine1 = ups / vels_pat
   sine1   = vps / vels_pat

   vtscr = ustar * patch_area

   ! For TEB
   ! For TEB is necessary to jump the lines above
   IF (TEB_SPM /=1) then
      sflux_u = sflux_u - ustar * vtscr * cosine1
      sflux_v = sflux_v - ustar * vtscr * sine1
      sflux_t = sflux_t - tstar * vtscr
      sflux_r = sflux_r - rstar * vtscr
   else
      if (NINT(G_URBAN) == 0) then
         sflux_u = sflux_u - ustar * vtscr * cosine1
         sflux_v = sflux_v - ustar * vtscr * sine1
         sflux_t = sflux_t - tstar * vtscr
         sflux_r = sflux_r - rstar * vtscr
      end if
   end if


   zoverl = gzotheta * vonk * tstar / (ustar * ustar)

   if (zoverl < 0.)then
      cx = zoverl * sqrt(sqrt(1. - 15. * zoverl))
   else
      cx = zoverl / (1.0 + 4.7 * zoverl)
   endif

   psin = sqrt((1.-2.86 * cx) / (1. + cx * (-5.39 + cx * 6.998 )))
   sflux_w = sflux_w + (0.27 * max(6.25 * (1. - cx) * psin,wtol)  &
      - 1.18 * cx * psin) * ustar * vtscr

   return
end subroutine sfclmcv

!*****************************************************************************

subroutine grndvap(soil_energy,soil_water,soil_text,sfcw_energy_int  &
   ,sfcwater_nlev,ground_rsat,ground_rvap,can_temp,can_rvap,prsg)

  use leaf_coms
  use rconstants
  use therm_lib, only: rslif,qwtk,qtk

  implicit none

  real :: soil_energy,soil_water,soil_text,sfcw_energy_int
  real :: sfcwater_nlev,ground_rsat,ground_rvap,can_temp,can_rvap,prsg

  integer :: ksn,nsoil

  real :: tempkk,fracliqq,slpotvn,alpha,beta


  ksn = nint(sfcwater_nlev)

  ! ground_rsat(i,j) is the saturation mixing ratio of the top soil/snow surface
  ! and is used for dew formation and snow evaporation.

  if (ksn > 0 .and. sfcw_energy_int > 0.) then
     call qtk(sfcw_energy_int,tempkk,fracliqq)
     ground_rsat = rslif(prsg,tempkk)
  else

  ! Without snowcover, ground_rvap is the effective saturation mixing
  ! ratio of soil and is used for soil evaporation.  First, compute the
  ! "alpha" term or soil "relative humidity" and the "beta" term.

     nsoil = nint(soil_text)

     call qwtk(soil_energy,soil_water*wdns,slcpd(nsoil),tempkk,fracliqq)
     ground_rsat = rslif(prsg,tempkk)

     slpotvn = slpots(nsoil) * (slmsts(nsoil) / soil_water) ** slbs(nsoil)
     alpha = exp(gorm * slpotvn / tempkk)
     beta = .25 * (1. - cos (min(1.,soil_water / sfldcap(nsoil)) * pi1)) ** 2
     ground_rvap = ground_rsat * alpha * beta + (1. - beta) * can_rvap

  end if

  return
end subroutine grndvap

!*****************************************************************************

subroutine leaf_bcond(m2,m3,mzg,mzs,npat,jdim  &

   ,soil_water       ,sfcwater_mass  ,soil_energy     &
   ,sfcwater_energy  ,soil_text      ,sfcwater_depth  &
   ,ustar            ,tstar          ,rstar           &
   ,veg_albedo       ,veg_fracarea   ,veg_lai         &
   ,veg_tai                                           &
   ,veg_rough        ,veg_height     ,patch_area      &
   ,patch_rough      ,patch_wetind   ,leaf_class      &
   ,soil_rough       ,sfcwater_nlev  ,stom_resist     &
   ,ground_rsat      ,ground_rvap    ,veg_water       &
   ,veg_temp         ,can_rvap       ,can_temp        &
   ,veg_ndvip        ,veg_ndvic      ,veg_ndvif)

   implicit none

   integer :: m2,m3,mzg,mzs,npat,jdim

   real, dimension(mzg,m2,m3,npat) :: soil_water,soil_energy,soil_text
   real, dimension(mzs,m2,m3,npat) :: sfcwater_mass,sfcwater_energy            &
                                     ,sfcwater_depth
   real, dimension(m2,m3,npat)     :: ustar        ,tstar         ,rstar        &
                                     ,veg_albedo   ,veg_fracarea  ,veg_lai      &
                                     ,veg_tai                                   &
                                     ,veg_rough    ,veg_height    ,patch_area   &
                                     ,patch_rough  ,patch_wetind  ,leaf_class   &
                                     ,soil_rough   ,sfcwater_nlev ,stom_resist  &
                                     ,ground_rsat  ,ground_rvap   ,veg_water    &
                                     ,veg_temp     ,can_rvap      ,can_temp     &
                                     ,veg_ndvip    ,veg_ndvic     ,veg_ndvif

   integer :: i,j,k,ipat

   do ipat = 1,npat
      do j = 1,m3

         ustar          (1,j,ipat) = ustar            (2,j,ipat)
         tstar          (1,j,ipat) = tstar            (2,j,ipat)
         rstar          (1,j,ipat) = rstar            (2,j,ipat)
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
         veg_temp       (1,j,ipat) = veg_temp         (2,j,ipat)
         can_rvap       (1,j,ipat) = can_rvap         (2,j,ipat)
         can_temp       (1,j,ipat) = can_temp         (2,j,ipat)
         veg_ndvip      (1,j,ipat) = veg_ndvip        (2,j,ipat)
         veg_ndvic      (1,j,ipat) = veg_ndvic        (2,j,ipat)
         veg_ndvif      (1,j,ipat) = veg_ndvif        (2,j,ipat)

         ustar         (m2,j,ipat) = ustar         (m2-1,j,ipat)
         tstar         (m2,j,ipat) = tstar         (m2-1,j,ipat)
         rstar         (m2,j,ipat) = rstar         (m2-1,j,ipat)
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
         veg_temp      (m2,j,ipat) = veg_temp      (m2-1,j,ipat)
         can_rvap      (m2,j,ipat) = can_rvap      (m2-1,j,ipat)
         can_temp      (m2,j,ipat) = can_temp      (m2-1,j,ipat)
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
         enddo

         do k = 1,mzs
            sfcwater_mass    (k,1,j,ipat) = sfcwater_mass      (k,2,j,ipat)
            sfcwater_energy  (k,1,j,ipat) = sfcwater_energy    (k,2,j,ipat)
            sfcwater_depth   (k,1,j,ipat) = sfcwater_depth     (k,2,j,ipat)

            sfcwater_mass   (k,m2,j,ipat) = sfcwater_mass   (k,m2-1,j,ipat)
            sfcwater_energy (k,m2,j,ipat) = sfcwater_energy (k,m2-1,j,ipat)
            sfcwater_depth  (k,m2,j,ipat) = sfcwater_depth  (k,m2-1,j,ipat)
         enddo

      enddo

      if (jdim == 1) then

         do i = 1,m2
            ustar          (i,1,ipat) = ustar            (i,2,ipat)
            tstar          (i,1,ipat) = tstar            (i,2,ipat)
            rstar          (i,1,ipat) = rstar            (i,2,ipat)
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
            veg_temp       (i,1,ipat) = veg_temp         (i,2,ipat)
            can_rvap       (i,1,ipat) = can_rvap         (i,2,ipat)
            can_temp       (i,1,ipat) = can_temp         (i,2,ipat)
            veg_ndvip      (i,1,ipat) = veg_ndvip        (i,2,ipat)
            veg_ndvic      (i,1,ipat) = veg_ndvic        (i,2,ipat)
            veg_ndvif      (i,1,ipat) = veg_ndvif        (i,2,ipat)

            ustar         (i,m3,ipat) = ustar         (i,m3-1,ipat)
            tstar         (i,m3,ipat) = tstar         (i,m3-1,ipat)
            rstar         (i,m3,ipat) = rstar         (i,m3-1,ipat)
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
            veg_temp      (i,m3,ipat) = veg_temp      (i,m3-1,ipat)
            can_rvap      (i,m3,ipat) = can_rvap      (i,m3-1,ipat)
            can_temp      (i,m3,ipat) = can_temp      (i,m3-1,ipat)
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
            enddo

            do k = 1,mzs
               sfcwater_mass    (k,i,1,ipat) = sfcwater_mass      (k,i,2,ipat)
               sfcwater_energy  (k,i,1,ipat) = sfcwater_energy    (k,i,2,ipat)
               sfcwater_depth   (k,i,1,ipat) = sfcwater_depth     (k,i,2,ipat)

               sfcwater_mass   (k,i,m3,ipat) = sfcwater_mass   (k,i,m3-1,ipat)
               sfcwater_energy (k,i,m3,ipat) = sfcwater_energy (k,i,m3-1,ipat)
               sfcwater_depth  (k,i,m3,ipat) = sfcwater_depth  (k,i,m3-1,ipat)
            enddo

         enddo

      endif

   enddo
   return
end subroutine leaf_bcond

!****************************************************************************

subroutine sfc_fields(m1,m2,m3,ia,iz,ja,jz,jd  &
   ,theta,rv,up,vp,dn0,pp,pi0,rtgt,zt,ths2,rvs2,ups2,vps2,pis2,dens2,zts2)

   use leaf_coms
   use rconstants

   implicit none

   integer :: m1,m2,m3,ia,iz,ja,jz,jd
   real, dimension(m1,m2,m3) :: theta,rv,up,vp,dn0,pp,pi0
   real, dimension(m2,m3) :: rtgt,ths2,rvs2,ups2,vps2,pis2,dens2,zts2
   real, dimension(m1) :: zt

   integer :: i,j
   real :: hcpi

   ! Compute surface atmospheric conditions

   hcpi = .5 * cpi

   do j = ja,jz
      do i = ia,iz
         ths2(i,j) = theta(2,i,j)
         rvs2(i,j) = rv(2,i,j)
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
   ,topma,aru,arv,theta,rv,up,vp,dn0,pp,pi0,zt,zm,dzt       &
   ,ths2,rvs2,ups2,vps2,pis2,dens2,zts2)

   use leaf_coms
   use rconstants

   implicit none

   integer :: m1,m2,m3,ia,iz,ja,jz,jd
   real, dimension(m2,m3) :: flpu,flpv,flpw
   real, dimension(m1,m2,m3) :: aru,arv,theta,rv,up,vp,dn0,pp,pi0
   real, dimension(m2,m3) :: topma,ths2,rvs2,ups2,vps2,pis2,dens2,zts2
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

         ths2(i,j) =  wtw * theta(k2,i,j) + (1. - wtw)  * theta(k3,i,j)

         rvs2(i,j) =  wtw * rv(k2,i,j)    + (1. - wtw)  * rv(k3,i,j)

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
      qpcpgl = pcpgl  * cliq * (ths * pis - tsupercool)
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

! Remaining issues:
!
! 1. Relationship between clumping, V, vegfrac
! 2. Impact of V on radiation
! 3. Put in 900 second response time for rc?
! 4. Build lookup tables, especially for things with exponentials?

!****************************************************************************

subroutine sfcrad(mzg,mzs,ip  &
   ,soil_energy,soil_water,soil_text,sfcwater_energy,sfcwater_mass    &
   ,sfcwater_depth,patch_area,can_temp,veg_temp,leaf_class,veg_height &
   ,veg_fracarea,veg_albedo,sfcwater_nlev,rshort,rlong,albedt         &
   ,rlongup,cosz,G_URBAN, ETOWN, ALBTOWN, TSTOWN )
   use node_mod       , only : mynum       ! ! intent(in)
   use mem_leaf
   use leaf_coms
   use rconstants
   use mem_scratch
   use therm_lib, only : qwtk,qtk
   !CATT
   use catt_start, only: CATT ! INTENT(IN)

   !TEB_SPM
   use teb_spm_start, only: TEB_SPM ! INTENT(IN)

   implicit none

   integer, intent(in) :: mzg,mzs,ip
   real, dimension(mzg) :: soil_energy,soil_water,soil_text !(mzg)
   real, dimension(mzs) :: sfcwater_energy,sfcwater_depth ,sfcwater_mass   !(mzs)
   real :: patch_area,can_temp,veg_temp,leaf_class,veg_height,veg_fracarea  &
          ,veg_albedo,sfcwater_nlev,rshort,rlong,albedt,rlongup,cosz

   ! for TEB
   real :: G_URBAN, ETOWN, ALBTOWN, TSTOWN
   !

   integer :: k,m,nsoil,nveg,ksn
   real :: alb,vfc,fcpct,alg,rad,als,fractrans,absg,algs,emv,emgs,gslong,vlong,alv

   ! This routine is called by the radiation parameterization and by leaf.
   ! It computes net surface albedo plus radiative exchange between the
   ! atmosphere, vegetation, and the snow/ground given previously computed
   ! downward longwave and shortwave fluxes from the atmosphere.
   ! Also computed are functions of snowcover that are required for the above
   ! radiation calculations as well as other calculations in leaf.

   ! The shortwave parameterizations are only valid if the cosine of the zenith
   ! angle is greater than .03 .  Water albedo from Atwater and Bell (1981)

   ! alg, als, and alv are the albedos of the ground, snow, and vegetation
   ! (als needs a better formula based on age of the surface snow).

   ! absg and vctr32 are the actual fractions of shortwave incident on snow
   ! plus ground that get absorbed by the ground and each snow layer,
   ! respectively.  They currently use the variable fractrans, which is the
   ! fraction of light transmitted through each layer based on mass per square
   ! meter.  algs is the resultant albedo from snow plus ground.

   if (ip == 1) then
      if (cosz > .03) then
         alb = min(max(-.0139 + .0467 * tan(acos(cosz)),.03),.999)
         albedt = albedt + patch_area * alb
      endif

      call qtk(soil_energy(mzg),tempk(mzg),fracliq(mzg))
      rlongup = rlongup + patch_area * stefan * tempk(mzg) ** 4

   elseif (isfcl == 0) then

      albedt = albedt + patch_area * albedo
      rlongup = rlongup + patch_area * stefan * can_temp ** 4

   else

   ! Diagnose soil temperature and liquid fraction

      do k = 1,mzg
         nsoil = nint(soil_text(k))
         call qwtk(soil_energy(k),soil_water(k)*wdns  &
            ,slcpd(nsoil),tempk(k),fracliq(k))
      enddo

   ! Diagnose snow temperature and the influence of snow covering veg.

      nveg = nint(leaf_class)
      ksn = nint(sfcwater_nlev)
      snowfac = 0.
      do k = 1,ksn
         snowfac = snowfac + sfcwater_depth(k)
         call qtk(sfcwater_energy(k),tempk(k+mzg),fracliq(k+mzg))
      enddo
      snowfac = min(.99, snowfac / max(.001,veg_height))

      vf = veg_fracarea * (1. - snowfac)
      vfc = 1. - vf

   ! Shortwave radiation calculations

      !srf-25-02-2005
      nsoil=nint(soil_text(mzg))
      !srf

      fcpct = soil_water(mzg) / slmsts(nsoil)
      if (fcpct > .5) then
         alg = .14
      else
         alg = .31 - .34 * fcpct
      endif
      alv = veg_albedo

      rad = 1.
      if (ksn > 0) then

   ! als = .14 (the wet soil value) for all-liquid

         als = .5 - .36 * fracliq(ksn+mzg)
         rad = 1. - als
      endif
      do k = ksn,1,-1
         fractrans = exp(-20. * sfcwater_depth(k))
         vctr32(k) = rad * (1. - fractrans)
         rad = rad * fractrans
      enddo
      absg = (1. - alg) * rad
      algs = 1. - absg
      do k = ksn,1,-1
         algs = algs - vctr32(k)
         rshort_s(k) = rshort * vfc * vctr32(k)
      enddo
      rshort_g = rshort * vfc * absg
      rshort_v = rshort * vf * (1. - alv + vfc * algs)
   !  rshort_a = rshort * (vf * alv + vfc * vfc * algs)

      alb = vf * alv + vfc * vfc * algs

      ! TEB
      if (TEB_SPM==1) then
         if (NINT(G_URBAN) == 0) then
            albedt = albedt + patch_area * alb
         else
            albedt = albedt + patch_area * ALBTOWN
         endif
      else
         albedt = albedt + patch_area * alb
      endif

   ! Longwave radiation calculations

      emv = emisv(nveg)
      emgs = emisg(nsoil)
      if (ksn > 0) emgs = 1.0
      gslong = emgs * stefan * tempk(ksn+mzg) ** 4
      vlong = emv * stefan * veg_temp ** 4

      rlonga_v  = rlong * vf * (emv + vfc * (1. - emgs))
      rlonga_gs = rlong * vfc * emgs
      rlongv_gs = vlong * vf * emgs
      rlongv_a  = vlong * vf * (2. - emgs - vf + emgs * vf)
      rlonggs_v = gslong * vf * emv
      rlonggs_a = gslong * vfc
      rlonga_a  = rlong * (vf * (1. - emv) + vfc * vfc * (1. - emgs))

      ! TEB
      if (TEB_SPM==1) then
         IF (NINT(G_URBAN) == 0) THEN
            rlongup = rlongup + patch_area * (rlongv_a + rlonggs_a + rlonga_a)
         else
            rlongup = rlongup + patch_area * ETOWN * STEFAN * TSTOWN**4
         ENDIF
      else
         rlongup = rlongup + patch_area * (rlongv_a + rlonggs_a + rlonga_a)
      endif

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! In case rlong is not computed, zero out all longwave fluxes other
   ! than rlongup.  [On the first timestep, radiative fluxes may not
   ! be available until microphysics is called, and zeroing these fluxes
   ! avoids the imbalance of having upward longwave without downward longwave
   ! in LEAF-2.  Also, this allows LEAF-2 to run without radiation for all
   ! timesteps, if desired for testing.]

      if (rlong < .1) then
         rlonga_v  = 0.
         rlonga_gs = 0.
         rlongv_gs = 0.
         rlongv_a  = 0.
         rlonggs_v = 0.
         rlonggs_a = 0.
         rlonga_a  = 0.
      endif
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !------------------------------------------------------------------------------------!
      !   800 W/m² corresponds to an emission at roughly 71°C, which is too hot. In case   !
      ! this happens, we'd rather stop the run and check what is going on. It will write   !
      ! the output in a fort.?? file.                                                      !
      !------------------------------------------------------------------------------------!
      if (rlongup > 800.) then
         m=mynum+20
         write (unit=m,fmt='(80a)')         ('-',k=1,80)
         write (unit=m,fmt='(a)')           '    NON-SENSE RLONGUP:'
         write (unit=m,fmt='(80a)')         ('-',k=1,80)
         write (unit=m,fmt='(a,1x,i6)')     ' + PATCH:      ',ip
         write (unit=m,fmt='(a,1x,i6)')     ' + SOIL_TEXT:  ',nsoil
         write (unit=m,fmt='(a,1x,i6)')     ' + KSN:        ',ksn
         write (unit=m,fmt='(a,1x,es12.5)') ' + RLONGUP:    ',rlongup
         write (unit=m,fmt='(a,1x,es12.5)') ' + PATCH_AREA: ',patch_area
         write (unit=m,fmt='(a,1x,es12.5)') ' + ALBEDT:     ',albedt
         write (unit=m,fmt='(a,1x,es12.5)') ' + RLONGV_A:   ',rlongv_a
         write (unit=m,fmt='(a,1x,es12.5)') ' + RLONGGS_A:  ',rlonggs_a
         write (unit=m,fmt='(a,1x,es12.5)') ' + RLONGA_A:   ',rlonga_a
         write (unit=m,fmt='(a,1x,es12.5)') ' + VEG_TEMP:   ',veg_temp
         write (unit=m,fmt='(a,1x,es12.5)') ' + CAN_TEMP:   ',can_temp
         write (unit=m,fmt='(80a)')         ('-',k=1,80)
         write (unit=m,fmt='(a3,1x,4(a12,1x))') 'LEV',' SOIL_ENERGY','  SOIL_WATER'        &
                                                     ,'       TEMPK','     FRACLIQ'
         write (unit=m,fmt='(80a)')         ('-',k=1,80)
         do k=1,mzg
            write (unit=m,fmt='(i3,1x,4(es12.5,1x))') k,soil_energy(k),soil_water(k)       &
                                                       ,tempk(k),fracliq(k)
         end do
         write (unit=m,fmt='(80a)')         ('-',k=1,80)
         if (ksn > 0) then
            write (unit=m,fmt='(a3,1x,5(a12,1x))') 'LEV',' SFCW_ENERGY','  SFCW_DEPTH'     &
                                         ,'   SFCW_MASS','       TEMPK','     FRACLIQ'
            write (unit=m,fmt='(80a)')         ('-',k=1,80)
            do k=1,ksn
               write (unit=m,fmt='(i3,1x,5(es12.5,1x))')  k,sfcwater_energy(k)             &
                                                           ,sfcwater_depth(k)              &
                                                           ,sfcwater_mass(k)               &
                                                           ,tempk(mzg+k),fracliq(mzg+k)
            end do
            write (unit=m,fmt='(80a)')         ('-',k=1,80)
         end if
         call abort_run('Non-sense rlongup','sfcrad','leaf3.f90')
      end if

   endif
   return
end subroutine sfcrad

