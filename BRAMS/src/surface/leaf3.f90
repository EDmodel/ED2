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
subroutine leaf3_timestep()
   use mem_grid       , only : nstbot            & ! intent(in)
                             , ngrid             & ! intent(in)
                             , grid_g            & ! intent(in)
                             , time              & ! intent(in)
                             , dtlt              & ! intent(in)
                             , dzt               & ! intent(in)
                             , zt                & ! intent(in)
                             , zm                & ! intent(in)
                             , nzg               & ! intent(in)
                             , nzs               & ! intent(in)
                             , istp              & ! intent(in)
                             , npatch            & ! intent(in)
                             , jdim              & ! intent(in)
                             , itimea            & ! intent(in)
                             , if_adap           ! ! intent(in)
   use io_params      , only : iupdsst           & ! intent(in)
                             , iupdndvi          & ! intent(in)
                             , ssttime1          & ! intent(in)
                             , ssttime2          & ! intent(in)
                             , ndvitime1         & ! intent(in)
                             , ndvitime2         ! ! intent(in)
   use teb_spm_start  , only : teb_spm           ! ! intent(in)
   use mem_teb        , only : teb_g             & ! data type
                             , teb_vars          ! ! type
   use mem_teb_common , only : tebc_g            & ! data type
                             , teb_common        ! ! type
   use mem_scratch    , only : scratch           ! ! intent(inout)
   use node_mod       , only : mzp               & ! intent(in)
                             , mxp               & ! intent(in)
                             , myp               & ! intent(in)
                             , ia                & ! intent(in)
                             , iz                & ! intent(in)
                             , ja                & ! intent(in)
                             , jz                & ! intent(in)
                             , ibcon             & ! intent(in)
                             , mynum             ! ! intent(in)
   use mem_leaf       , only : dtleaf            & ! intent(in)
                             , isfcl             & ! intent(in)
                             , leaf_g            ! ! intent(inout)
   use leaf_coms      , only : timefac_sst       & ! intent(in)
                             , timefac_ndvi      & ! intent(in)
                             , niter_leaf        & ! intent(in)
                             , dtll_factor       & ! intent(in)
                             , dtll              & ! intent(in)
                             , atm_theta         & ! intent(out)
                             , atm_enthalpy      & ! intent(out)
                             , atm_shv           & ! intent(out)
                             , atm_rvap          & ! intent(out)
                             , atm_co2           & ! intent(out)
                             , can_enthalpy      & ! intent(out)
                             , can_shv           & ! intent(out)
                             , can_rhos          & ! intent(out)
                             , geoht             & ! intent(out)
                             , atm_up            & ! intent(out)
                             , atm_vp            & ! intent(out)
                             , atm_vels          & ! intent(out)
                             , pcpgl             & ! intent(out)
                             , estar             & ! intent(out)
                             , qstar             & ! intent(out)
                             , min_patch_area    & ! intent(out)
                             , g_urban           & ! intent(out)
                             , emis_town         & ! intent(out)
                             , alb_town          & ! intent(out)
                             , ts_town           & ! intent(out)
                             , ggbare            & ! intent(out)
                             , flush_leaf_coms   ! ! sub-routine
   use mem_basic      , only : basic_g           & ! intent(in)
                             , co2_on            & ! intent(in)
                             , co2con            ! ! intent(in)
   use mem_turb       , only : turb_g            ! ! intent(inout)
   use mem_cuparm     , only : cuparm_g          & ! intent(in)
                             , nnqparm           & ! intent(in)
                             , nclouds           ! ! intent(in)
   use mem_micro      , only : micro_g           ! ! intent(in)
   use mem_radiate    , only : radiate_g         & ! intent(inout)
                             , iswrtyp           & ! intent(in)
                             , ilwrtyp           ! ! intent(in)
   use therm_lib      , only : bulk_on           & ! intent(in)
                             , press2exner       & ! function
                             , exner2press       & ! function
                             , extheta2temp      ! ! function
   use rconstants     , only : alvl3             & ! intent(in)
                             , cpdry             ! ! intent(in)
   implicit none
   !----- Local variables. ----------------------------------------------------------------!
   integer  :: i
   integer  :: j
   integer  :: ip
   integer  :: iter_leaf
   real     :: psup1
   real     :: psup2
   real     :: depe
   real     :: alt2
   real     :: deze
   real     :: dpdz
   real     :: exn1st
   real     :: airt
   real     :: zh_town
   real     :: zle_town
   real     :: zsfu_town
   real     :: zsfv_town
   !---------------------------------------------------------------------------------------!
   
   !------Nothing to do here if the bottom is not at the ground. --------------------------!
   if (nstbot == 0) return
   !---------------------------------------------------------------------------------------!



   !----- This will reset some variables, so we make sure they are properly initialised. --!
   call flush_leaf_coms('INITIAL')
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Define LEAF-3 time-split timesteps here.  This ensures that LEAF-3 will never use  !
   ! a timestep longer than about 30 seconds, but the actual number depends on the user's  !
   ! choice and the actual time step.                                                      !
   !---------------------------------------------------------------------------------------!
   niter_leaf  = max(1,nint(dtlt/dtleaf + .4))
   dtll_factor = 1. / float(niter_leaf)
   dtll        = dtlt * dtll_factor
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Here we copy a few variables to scratch arrays, as they may not be always exist. !
   ! In case they don't, we fill the scratch arrays with the default values.  The follow-  !
   ! ing scratch arrays will contain the following fields.                                 !
   !                                                                                       !
   ! vt3do => CO2 mixing ratio                                                             !
   ! vt3dp => precipitation rate from cumulus parametrisation.                             !
   ! vt2dq => precipitation rate from bulk microphysics                                    !
   ! vt2dr => internal energy associated with precipitation rate from bulk microphysics    !
   ! vt2ds => depth associated with precipitation rate from bulk microphysics              !
   !---------------------------------------------------------------------------------------!
   !----- Check whether we have CO2, and copy to an scratch array. ------------------------!
   if (co2_on) then
      call atob(mzp*mxp*myp,basic_g(ngrid)%co2p,scratch%vt3do)
   else
      call ae0(mzp*mxp*myp,scratch%vt3do,co2con(1))
   end if
   !----- Check whether cumulus parametrisation was used, and copy to scratch array. ------!
   if (nnqparm(ngrid) /= 0) then
      call atob(mxp*myp*nclouds,cuparm_g(ngrid)%conprr,scratch%vt3dp)
   else
      call azero(mxp*myp*nclouds,scratch%vt3dp)
   end if
   !----- Check whether bulk microphysics was used, and copy values to scratch array. -----!
   if (bulk_on) then
      call atob(mxp*myp,micro_g(ngrid)%pcpg ,scratch%vt2dq)
      call atob(mxp*myp,micro_g(ngrid)%qpcpg,scratch%vt2dr)
      call atob(mxp*myp,micro_g(ngrid)%dpcpg,scratch%vt2ds)
   else
      call azero(mxp*myp,scratch%vt2dq)
      call azero(mxp*myp,scratch%vt2dr)
      call azero(mxp*myp,scratch%vt2ds)
   end if 
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Copy surface atmospheric variables into 2-D arrays for input to LEAF.  The 2-D    !
   ! arrays are save as the following:                                                     !!
   !                                                                                       !
   ! vt2da => ice-liquid potential temperature                                             !
   ! vt2db => potential temperature                                                        !
   ! vt2dc => water vapour mixing ratio                                                    !
   ! vt2dd => total water mixing ratio (ice + liquid + vapour)                             !
   ! vt2de => CO2 mixing ratio                                                             !
   ! vt2df => zonal wind speed                                                             !
   ! vt2dg => meridional wind speed                                                        !
   ! vt2dh => Exner function                                                               !
   ! vt2di => Air density                                                                  !
   ! vt2dj => Reference height                                                             !
   ! vt2dk => Precipitation rate                                                           !
   ! vt2dl => Internal energy of precipitation rate                                        !
   ! vt2dm => Depth associated with the precipitation rate                                 !
   !---------------------------------------------------------------------------------------!
   select case (if_adap)
   case (0)
      call sfc_fields( mzp,mxp,myp,ia,iz,ja,jz,jdim                                        &
                     , basic_g(ngrid)%thp   , basic_g(ngrid)%theta , basic_g(ngrid)%rv     &
                     , basic_g(ngrid)%rtp   , scratch%vt3do        , basic_g(ngrid)%up     &
                     , basic_g(ngrid)%vp    , basic_g(ngrid)%dn0   , basic_g(ngrid)%pp     &
                     , basic_g(ngrid)%pi0   , grid_g(ngrid)%rtgt   , zt                    &
                     , zm                   , scratch%vt2da        , scratch%vt2db         &
                     , scratch%vt2dc        , scratch%vt2dd        , scratch%vt2de         &
                     , scratch%vt2df        , scratch%vt2dg        , scratch%vt2dh         &
                     , scratch%vt2di        , scratch%vt2dj                                ) 
   case (1)
      call sfc_fields_adap(mzp,mxp,myp,ia,iz,ja,jz,jdim                                    &
                     , grid_g(ngrid)%flpu   , grid_g(ngrid)%flpv   , grid_g(ngrid)%flpw    &
                     , grid_g(ngrid)%topma  , grid_g(ngrid)%aru    , grid_g(ngrid)%arv     &
                     , basic_g(ngrid)%thp   , basic_g(ngrid)%theta , basic_g(ngrid)%rv     &
                     , basic_g(ngrid)%rtp   , scratch%vt3do        , basic_g(ngrid)%up     &
                     , basic_g(ngrid)%vp    , basic_g(ngrid)%dn0   , basic_g(ngrid)%pp     &
                     , basic_g(ngrid)%pi0   , zt                   , zm                    &
                     , dzt                  , scratch%vt2da        , scratch%vt2db         &
                     , scratch%vt2dc        , scratch%vt2dd        , scratch%vt2de         &
                     , scratch%vt2df        , scratch%vt2dg        , scratch%vt2dh         &
                     , scratch%vt2di        , scratch%vt2dj                                )
   end select
   !---------------------------------------------------------------------------------------!



   !----- Fill surface precipitation arrays for input to LEAF-3 ---------------------------!
   call sfc_pcp(mxp,myp,nclouds,ia,iz,ja,jz,dtll,dtll_factor,scratch%vt2db,scratch%vt2dh   &
               ,scratch%vt3dp,scratch%vt2dq,scratch%vt2dr,scratch%vt2ds,scratch%vt2dk      &
               ,scratch%vt2dl,scratch%vt2dm)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Reset fluxes, albedo, and upwelling long-wave radiation.                          !
   !---------------------------------------------------------------------------------------!
   call azero(mxp*myp       ,turb_g(ngrid)%sflux_u    )
   call azero(mxp*myp       ,turb_g(ngrid)%sflux_v    )
   call azero(mxp*myp       ,turb_g(ngrid)%sflux_w    )
   call azero(mxp*myp       ,turb_g(ngrid)%sflux_t    )
   call azero(mxp*myp       ,turb_g(ngrid)%sflux_r    )
   call azero(mxp*myp       ,turb_g(ngrid)%sflux_c    )
   call azero(mxp*myp*npatch,leaf_g(ngrid)%sensible_gc)
   call azero(mxp*myp*npatch,leaf_g(ngrid)%sensible_vc)
   call azero(mxp*myp*npatch,leaf_g(ngrid)%evap_gc    )
   call azero(mxp*myp*npatch,leaf_g(ngrid)%evap_vc    )
   call azero(mxp*myp*npatch,leaf_g(ngrid)%transp     )
   call azero(mxp*myp*npatch,leaf_g(ngrid)%gpp        )
   call azero(mxp*myp*npatch,leaf_g(ngrid)%plresp     )
   call azero(mxp*myp*npatch,leaf_g(ngrid)%resphet    )
   call azero(mxp*myp*npatch,leaf_g(ngrid)%rshort_gnd )
   call azero(mxp*myp*npatch,leaf_g(ngrid)%rlong_gnd  )
   if (iswrtyp > 0 .or. ilwrtyp > 0) then
      call azero(mxp*myp,radiate_g(ngrid)%albedt)
      call azero(mxp*myp,radiate_g(ngrid)%rlongup)
   end if
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Big loop over the horizontal domain.                                              !
   !---------------------------------------------------------------------------------------!
   latloop: do j=ja,jz
      lonloop: do i=ia,iz


         !---------------------------------------------------------------------------------!
         !      This will reset some variables, so we make sure they are properly          !
         ! initialised for every grid point.                                               !
         !---------------------------------------------------------------------------------!
         call flush_leaf_coms('GRID_POINT')
         !---------------------------------------------------------------------------------!


         !----- Copy the surface variables to single column values. -----------------------!
         call leaf3_atmo1d(mxp,myp,i,j,scratch%vt2da,scratch%vt2db,scratch%vt2dc           &
                          ,scratch%vt2dd,scratch%vt2de,scratch%vt2df,scratch%vt2dg         &
                          ,scratch%vt2dh,scratch%vt2di, scratch%vt2dj,scratch%vt2dk        &
                          ,scratch%vt2dl,scratch%vt2dm)
         !---------------------------------------------------------------------------------!





         !----- For no soil model (patch 2) fill "canopy" temperature and moisture. -------!
         if (isfcl == 0) then
            call leaf0(mxp,myp,npatch,i,j,leaf_g(ngrid)%can_theta,leaf_g(ngrid)%can_rvap   &
                      ,leaf_g(ngrid)%can_co2,leaf_g(ngrid)%can_prss                        &
                      ,leaf_g(ngrid)%can_theiv,leaf_g(ngrid)%patch_area)
         end if
         !---------------------------------------------------------------------------------!



         !---------------------------------------------------------------------------------!
         !      Begin the first big patch loop.                                            !
         !---------------------------------------------------------------------------------!
         patloop1: do ip = 1, npatch


            !------------------------------------------------------------------------------!
            !    If this patch is not to be solved, skip it.                               !
            !------------------------------------------------------------------------------!
            if (ip /= 1 .and. leaf_g(ngrid)%patch_area(i,j,ip) < min_patch_area) then
               cycle patloop1
            end if

            !------------------------------------------------------------------------------!
            !      This will reset some variables, so we make sure they are properly       !
            ! initialised for every grid point.                                            !
            !------------------------------------------------------------------------------!
            call flush_leaf_coms('PATCH')
            !------------------------------------------------------------------------------!
            


            !----- Update time-dependent SST, vegetation LAI and fractional coverage ------!
            if (ip == 1) then
               call leaf3_ocean_diag(ngrid,nzg                                             &
                                    ,leaf_g(ngrid)%seatp (i,j)                             &
                                    ,leaf_g(ngrid)%seatf (i,j)                             &
                                    ,leaf_g(ngrid)%soil_energy(:,i,j,1))
            elseif (isfcl >= 1) then
               call vegndvi( ngrid                                                         &
                           , leaf_g(ngrid)%patch_area  (i,j,ip)                            &
                           , leaf_g(ngrid)%leaf_class  (i,j,ip)                            &
                           , leaf_g(ngrid)%veg_fracarea(i,j,ip)                            &
                           , leaf_g(ngrid)%veg_lai     (i,j,ip)                            &
                           , leaf_g(ngrid)%veg_tai     (i,j,ip)                            &
                           , leaf_g(ngrid)%veg_rough   (i,j,ip)                            &
                           , leaf_g(ngrid)%veg_height  (i,j,ip)                            &
                           , leaf_g(ngrid)%veg_displace(i,j,ip)                            &
                           , leaf_g(ngrid)%veg_albedo  (i,j,ip)                            &
                           , leaf_g(ngrid)%veg_ndvip   (i,j,ip)                            &
                           , leaf_g(ngrid)%veg_ndvic   (i,j,ip)                            &
                           , leaf_g(ngrid)%veg_ndvif   (i,j,ip)                            &
                           , leaf_g(ngrid)%psibar_10d  (i,j,ip)                            )
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !     Initialise various canopy air space and temporary surface water vari-    !
            ! ables.                                                                       !
            !------------------------------------------------------------------------------!
            call leaf3_soilsfcw_diag(ip,nzg,nzs,leaf_g(ngrid)%soil_energy      (:,i,j,ip)  &
                                               ,leaf_g(ngrid)%soil_water       (:,i,j,ip)  &
                                               ,leaf_g(ngrid)%soil_text        (:,i,j,ip)  &
                                               ,leaf_g(ngrid)%sfcwater_nlev    (  i,j,ip)  &
                                               ,leaf_g(ngrid)%sfcwater_energy  (:,i,j,ip)  &
                                               ,leaf_g(ngrid)%sfcwater_mass    (:,i,j,ip)  &
                                               ,leaf_g(ngrid)%sfcwater_depth   (:,i,j,ip)  &
                                               ,.true.)

            call leaf3_solve_veg(ip,nzs,leaf_g(ngrid)%leaf_class               (  i,j,ip)  &
                                       ,leaf_g(ngrid)%veg_height               (  i,j,ip)  &
                                       ,leaf_g(ngrid)%patch_area               (  i,j,ip)  &
                                       ,leaf_g(ngrid)%veg_fracarea             (  i,j,ip)  &
                                       ,leaf_g(ngrid)%veg_tai                  (  i,j,ip)  &
                                       ,leaf_g(ngrid)%sfcwater_nlev            (  i,j,ip)  &
                                       ,leaf_g(ngrid)%sfcwater_depth           (:,i,j,ip)  &
                                       ,.true.)

            call leaf3_can_diag(ip ,leaf_g(ngrid)%can_theta        (i,j,ip)                &
                                   ,leaf_g(ngrid)%can_theiv        (i,j,ip)                &
                                   ,leaf_g(ngrid)%can_rvap         (i,j,ip)                &
                                   ,leaf_g(ngrid)%leaf_class       (i,j,ip)                &
                                   ,leaf_g(ngrid)%can_prss         (i,j,ip)                &
                                   ,.true.                                                 )

            call leaf3_veg_diag(leaf_g(ngrid)%veg_energy       (i,j,ip)                    &
                               ,leaf_g(ngrid)%veg_water        (i,j,ip)                    &
                               ,leaf_g(ngrid)%veg_hcap         (i,j,ip)                    )

            if (ip /= 1) then
               call leaf3_grndvap(leaf_g(ngrid)%soil_energy       (nzg,i,j,ip)             &
                                 ,leaf_g(ngrid)%soil_water        (nzg,i,j,ip)             &
                                 ,leaf_g(ngrid)%soil_text         (nzg,i,j,ip)             &
                                 ,leaf_g(ngrid)%sfcwater_energy   (  1,i,j,ip)             &
                                 ,leaf_g(ngrid)%sfcwater_nlev     (    i,j,ip)             &
                                 ,leaf_g(ngrid)%can_rvap          (    i,j,ip)             &
                                 ,leaf_g(ngrid)%can_prss          (    i,j,ip)             &
                                 ,leaf_g(ngrid)%ground_rsat       (    i,j,ip)             &
                                 ,leaf_g(ngrid)%ground_rvap       (    i,j,ip)             &
                                 ,leaf_g(ngrid)%ground_temp       (    i,j,ip)             &
                                 ,leaf_g(ngrid)%ground_fliq       (    i,j,ip)             )
            end if
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Begin the time step.                                                    !
            !------------------------------------------------------------------------------!
            tloop: do iter_leaf = 1,niter_leaf

               !---- If TEB is on, copy the values to single variables. -------------------!
               if (teb_spm == 1) then
                  g_urban   = leaf_g(ngrid)%g_urban(i,j,ip)
                  emis_town = tebc_g(ngrid)%emis_town(i,j)
                  alb_town  = tebc_g(ngrid)%alb_town(i,j)
                  ts_town   = tebc_g(ngrid)%ts_town(i,j)
               end if

               !---------------------------------------------------------------------------!
               !    Calculate radiative fluxes between atmosphere, vegetation, and ground/ !
               ! /snow based on already-computed downward shortwave and longwave fluxes    !
               ! from the atmosphere.  Fill tempk array with soil and snow temperature (C) !
               ! and fracliq array with liquid fraction of water content in soil and snow. !
               ! Other snowcover properties are also computed here.                        !
               !---------------------------------------------------------------------------!
               if (iswrtyp > 0 .or. ilwrtyp > 0) then
                  call leaf3_sfcrad(nzg,nzs,ip                                             &
                                   ,leaf_g(ngrid)%soil_water       (:,i,j,ip)              &
                                   ,leaf_g(ngrid)%soil_color       (  i,j,ip)              &
                                   ,leaf_g(ngrid)%soil_text        (:,i,j,ip)              &
                                   ,leaf_g(ngrid)%sfcwater_depth   (:,i,j,ip)              &
                                   ,leaf_g(ngrid)%patch_area       (  i,j,ip)              &
                                   ,leaf_g(ngrid)%veg_fracarea     (  i,j,ip)              &
                                   ,leaf_g(ngrid)%leaf_class       (  i,j,ip)              &
                                   ,leaf_g(ngrid)%veg_albedo       (  i,j,ip)              &
                                   ,leaf_g(ngrid)%sfcwater_nlev    (  i,j,ip)              &
                                   ,radiate_g(ngrid)%rshort        (  i,j   )              &
                                   ,radiate_g(ngrid)%rlong         (  i,j   )              &
                                   ,radiate_g(ngrid)%cosz          (  i,j   )              &
                                   ,radiate_g(ngrid)%albedt        (  i,j   )              &
                                   ,radiate_g(ngrid)%rlongup       (  i,j   )              &
                                   ,leaf_g(ngrid)%rshort_gnd       (  i,j,ip)              &
                                   ,leaf_g(ngrid)%rlong_gnd        (  i,j,ip)              )
               end if
               !---------------------------------------------------------------------------!


               !----- Find the roughness length. ------------------------------------------!
               call leaf3_roughness(ip,leaf_g(ngrid)%veg_fracarea (i,j,ip)                 &
                                      ,leaf_g(ngrid)%patch_area   (i,j,ip)                 &
                                      ,leaf_g(ngrid)%ustar        (i,j,ip)                 &
                                      ,grid_g(ngrid)%topzo        (i,j)                    &
                                      ,leaf_g(ngrid)%veg_rough    (i,j,ip)                 &
                                      ,leaf_g(ngrid)%soil_rough   (i,j,ip)                 &
                                      ,leaf_g(ngrid)%patch_rough  (i,j,ip)                 )
               !---------------------------------------------------------------------------!



               !----- Compute the characteristic scales. ----------------------------------!
               call leaf3_stars(atm_theta                                                  &
                               ,atm_enthalpy                                               &
                               ,atm_shv                                                    &
                               ,atm_rvap                                                   &
                               ,atm_co2                                                    &
                               ,leaf_g(ngrid)%can_theta   (i,j,ip)                         &
                               ,can_enthalpy                                               &
                               ,can_shv                                                    &
                               ,leaf_g(ngrid)%can_rvap    (i,j,ip)                         &
                               ,leaf_g(ngrid)%can_co2     (i,j,ip)                         &
                               ,geoht                                                      &
                               ,leaf_g(ngrid)%veg_displace(i,j,ip)                         &
                               ,atm_vels                                                   &
                               ,dtll                                                       &
                               ,leaf_g(ngrid)%patch_rough (i,j,ip)                         &
                               ,leaf_g(ngrid)%ustar       (i,j,ip)                         &
                               ,leaf_g(ngrid)%tstar       (i,j,ip)                         &
                               ,estar                                                      &
                               ,qstar                                                      &
                               ,leaf_g(ngrid)%rstar       (i,j,ip)                         &
                               ,leaf_g(ngrid)%cstar       (i,j,ip)                         &
                               ,leaf_g(ngrid)%zeta        (i,j,ip)                         &
                               ,leaf_g(ngrid)%ribulk      (i,j,ip)                         &
                               ,leaf_g(ngrid)%R_aer       (i,j,ip)                         )
               !---------------------------------------------------------------------------!




               !----- Find the turbulent fluxes of momentum, heat and moisture.  ----------!
               call leaf3_sfclmcv(leaf_g(ngrid)%ustar         (i,j,ip)                     &
                                 ,leaf_g(ngrid)%tstar         (i,j,ip)                     &
                                 ,leaf_g(ngrid)%rstar         (i,j,ip)                     &
                                 ,leaf_g(ngrid)%cstar         (i,j,ip)                     &
                                 ,leaf_g(ngrid)%zeta          (i,j,ip)                     &
                                 ,atm_vels                                                 &
                                 ,atm_up                                                   &
                                 ,atm_vp                                                   &
                                 ,leaf_g(ngrid)%patch_area    (i,j,ip)                     &
                                 ,turb_g(ngrid)%sflux_u       (i,j   )                     &
                                 ,turb_g(ngrid)%sflux_v       (i,j   )                     &
                                 ,turb_g(ngrid)%sflux_w       (i,j   )                     &
                                 ,turb_g(ngrid)%sflux_t       (i,j   )                     &
                                 ,turb_g(ngrid)%sflux_r       (i,j   )                     &
                                 ,turb_g(ngrid)%sflux_c       (i,j   )                     )
               !---------------------------------------------------------------------------!


               !----- Solve the canopy and soil, depending on which patch we are. ---------!
               select case (ip)
               case (1)
                  !----- Call the solver for ocean (water) patches. -----------------------!
                  call leaf3_ocean(nzg,leaf_g(ngrid)%ustar                       (i,j,ip)  &
                                      ,leaf_g(ngrid)%tstar                       (i,j,ip)  &
                                      ,leaf_g(ngrid)%rstar                       (i,j,ip)  &
                                      ,leaf_g(ngrid)%cstar                       (i,j,ip)  &
                                      ,leaf_g(ngrid)%patch_rough                 (i,j,ip)  &
                                      ,leaf_g(ngrid)%can_prss                    (i,j,ip)  &
                                      ,leaf_g(ngrid)%can_rvap                    (i,j,ip)  &
                                      ,leaf_g(ngrid)%can_co2                     (i,j,ip)  &
                                      ,leaf_g(ngrid)%sensible_gc                 (i,j,ip)  &
                                      ,leaf_g(ngrid)%evap_gc                     (i,j,ip)  &
                                      ,leaf_g(ngrid)%plresp                      (i,j,ip)  &
                                      ,leaf_g(ngrid)%ground_temp                 (i,j,ip)  &
                                      ,leaf_g(ngrid)%ground_rsat                 (i,j,ip)  &
                                      ,leaf_g(ngrid)%ground_rvap                 (i,j,ip)  &
                                      ,leaf_g(ngrid)%ground_fliq                 (i,j,ip)  )

                  !----- Update some diagnostic variables. --------------------------------!
                  call flush_leaf_coms('TIMESTEP')
                  call leaf3_ocean_diag(ngrid,nzg                                          &
                                             ,leaf_g(ngrid)%seatp                   (i,j)  &
                                             ,leaf_g(ngrid)%seatf                   (i,j)  &
                                             ,leaf_g(ngrid)%soil_energy         (:,i,j,1)  )

                  call leaf3_can_diag(ip ,leaf_g(ngrid)%can_theta                (i,j,ip)  &
                                         ,leaf_g(ngrid)%can_theiv                (i,j,ip)  &
                                         ,leaf_g(ngrid)%can_rvap                 (i,j,ip)  &
                                         ,leaf_g(ngrid)%leaf_class               (i,j,ip)  &
                                         ,leaf_g(ngrid)%can_prss                 (i,j,ip)  &
                                         ,.false.                                          )
                  call leaf3_veg_diag(leaf_g(ngrid)%veg_energy                   (i,j,ip)  &
                                     ,leaf_g(ngrid)%veg_water                    (i,j,ip)  &
                                     ,leaf_g(ngrid)%veg_hcap                     (i,j,ip)  )
               case default

                  if (teb_spm == 1) then

                     !----- TEB: Urban canopy parameterization starts here ----------------!
                     if (nint(leaf_g(ngrid)%g_urban(i,j,ip)) /= 0.) then

                        !---- Initialise a few variables. ---------------------------------!
                        psup1  = exner2press( basic_g(ngrid)%pi0(1,i,j)                    &
                                            + basic_g(ngrid)%pp (1,i,j))
                        psup2  = exner2press( basic_g(ngrid)%pi0(2,i,j)                    &
                                            + basic_g(ngrid)%pp (2,i,j))
                        depe   = psup2-psup1
                        alt2   = zt(1)*grid_g(ngrid)%rtgt(i,j)
                        deze   = geoht-alt2
                        dpdz   = depe/deze
                        exn1st = press2exner(psup2)
                        
                        airt= extheta2temp(exn1st,basic_g(ngrid)%theta(2,i,j))

                        g_urban   = leaf_g(ngrid)%g_urban(i,j,ip)

                        call leaf3_teb_interface(istp,dtlt,dtll                            &
                                                ,radiate_g(ngrid)%cosz       (i,j)         &
                                                ,geoht                                     &
                                                ,radiate_g(ngrid)%rlong      (i,j)         &
                                                ,radiate_g(ngrid)%rshort     (i,j)         &
                                                ,psup2                                     &
                                                ,airt                                      &
                                                ,atm_up                                    &
                                                ,atm_vp                                    &
                                                ,basic_g(ngrid)%rv         (2,i,j)         &
                                                ,pcpgl/dtlt                                &
                                                ,teb_g(ngrid)%fuso           (i,j)         &
                                                ,teb_g(ngrid)%t_canyon       (i,j)         &
                                                ,teb_g(ngrid)%r_canyon       (i,j)         &
                                                ,teb_g(ngrid)%ts_roof        (i,j)         &
                                                ,teb_g(ngrid)%ts_road        (i,j)         &
                                                ,teb_g(ngrid)%ts_wall        (i,j)         &
                                                ,teb_g(ngrid)%ti_road        (i,j)         &
                                                ,teb_g(ngrid)%ti_bld         (i,j)         &
                                                ,teb_g(ngrid)%ws_roof        (i,j)         &
                                                ,teb_g(ngrid)%ws_road        (i,j)         &
                                                ,teb_g(ngrid)%t_roof     (2:4,i,j)         &
                                                ,teb_g(ngrid)%t_road     (2:4,i,j)         &
                                                ,teb_g(ngrid)%t_wall     (2:4,i,j)         &
                                                ,zh_town                                   &
                                                ,zle_town                                  &
                                                ,tebc_g(ngrid)%emis_town     (i,j)         &
                                                ,zsfu_town                                 &
                                                ,zsfv_town                                 &
                                                ,tebc_g(ngrid)%ts_town       (i,j)         &
                                                ,tebc_g(ngrid)%alb_town      (i,j)         &
                                                ,nint(g_urban)                             &
                                                ,teb_g(ngrid)%h_traffic      (i,j)         &
                                                ,teb_g(ngrid)%h_industry     (i,j)         &
                                                ,teb_g(ngrid)%le_traffic     (i,j)         &
                                                ,teb_g(ngrid)%le_industry    (i,j)         &
                                                ,teb_g(ngrid)%t2m_town       (i,j)         &
                                                ,teb_g(ngrid)%r2m_town       (i,j)         &
                                                ,time                                      &
                                                ,itimea                                    &
                                                ,dpdz                                      &
                                                ,can_rhos                                  )
                        
                        turb_g(ngrid)%sflux_u(i,j) = turb_g(ngrid)%sflux_u(i,j)            &
                                          + leaf_g(ngrid)%patch_area(i,j,ip) * zsfu_town
                        turb_g(ngrid)%sflux_v(i,j) = turb_g(ngrid)%sflux_v(i,j)            &
                                          + leaf_g(ngrid)%patch_area(i,j,ip) * zsfv_town
                        turb_g(ngrid)%sflux_t(i,j) = turb_g(ngrid)%sflux_t(i,j)            &
                                          + leaf_g(ngrid)%patch_area(i,j,ip) * zh_town     &
                                          / (cpdry * can_rhos)
                        turb_g(ngrid)%sflux_r(i,j) = turb_g(ngrid)%sflux_r(i,j)            &
                                          + leaf_g(ngrid)%patch_area(i,j,ip) * zle_town    &
                                          / (alvl3 * can_rhos)
                     end if
                  end if

                  if (isfcl >= 1) then
                     !---------------------------------------------------------------------!
                     !     For soil model patches, update temperature and moisture of      !
                     ! soil, vegetation, and canopy                                        !
                     !---------------------------------------------------------------------!
                     if (leaf_g(ngrid)%patch_area(i,j,ip) >= min_patch_area) then
                        call leaf3_prognostic(ngrid,i,j,ip)
                     end if
                     !---------------------------------------------------------------------!
                  end if
                  !------------------------------------------------------------------------!




                  !------------------------------------------------------------------------!
                  !     Update all the diagnostic variables.                               !
                  !------------------------------------------------------------------------!
                  call flush_leaf_coms('TIMESTEP')
                  call leaf3_soilsfcw_diag(ip,nzg,nzs                                       &
                                          ,leaf_g(ngrid)%soil_energy        (  :,i,j,ip)   &
                                          ,leaf_g(ngrid)%soil_water         (  :,i,j,ip)   &
                                          ,leaf_g(ngrid)%soil_text          (  :,i,j,ip)   &
                                          ,leaf_g(ngrid)%sfcwater_nlev      (    i,j,ip)   &
                                          ,leaf_g(ngrid)%sfcwater_energy    (  :,i,j,ip)   &
                                          ,leaf_g(ngrid)%sfcwater_mass      (  :,i,j,ip)   &
                                          ,leaf_g(ngrid)%sfcwater_depth     (  :,i,j,ip)   &
                                          ,.false.)
                  call leaf3_solve_veg(ip,nzs,leaf_g(ngrid)%leaf_class          (i,j,ip)   &
                                             ,leaf_g(ngrid)%veg_height      (    i,j,ip)   &
                                             ,leaf_g(ngrid)%patch_area          (i,j,ip)   &
                                             ,leaf_g(ngrid)%veg_fracarea        (i,j,ip)   &
                                             ,leaf_g(ngrid)%veg_tai         (    i,j,ip)   &
                                             ,leaf_g(ngrid)%sfcwater_nlev   (    i,j,ip)   &
                                             ,leaf_g(ngrid)%sfcwater_depth  (  :,i,j,ip)   &
                                             ,.false.)

                  call leaf3_can_diag(ip ,leaf_g(ngrid)%can_theta           (    i,j,ip)   &
                                         ,leaf_g(ngrid)%can_theiv           (    i,j,ip)   &
                                         ,leaf_g(ngrid)%can_rvap            (    i,j,ip)   &
                                         ,leaf_g(ngrid)%leaf_class          (    i,j,ip)   &
                                         ,leaf_g(ngrid)%can_prss            (    i,j,ip)   &
                                         ,.false.                                          )

                  call leaf3_veg_diag(leaf_g(ngrid)%veg_energy               (    i,j,ip)   &
                                     ,leaf_g(ngrid)%veg_water               (    i,j,ip)   &
                                     ,leaf_g(ngrid)%veg_hcap                (    i,j,ip)   )

                  call leaf3_grndvap(leaf_g(ngrid)%soil_energy               (nzg,i,j,ip)   &
                                    ,leaf_g(ngrid)%soil_water               (nzg,i,j,ip)   &
                                    ,leaf_g(ngrid)%soil_text                (nzg,i,j,ip)   &
                                    ,leaf_g(ngrid)%sfcwater_energy          (  1,i,j,ip)   &
                                    ,leaf_g(ngrid)%sfcwater_nlev            (    i,j,ip)   &
                                    ,leaf_g(ngrid)%can_rvap                 (    i,j,ip)   &
                                    ,leaf_g(ngrid)%can_prss                 (    i,j,ip)   &
                                    ,leaf_g(ngrid)%ground_rsat              (    i,j,ip)   &
                                    ,leaf_g(ngrid)%ground_rvap              (    i,j,ip)   &
                                    ,leaf_g(ngrid)%ground_temp              (    i,j,ip)   &
                                    ,leaf_g(ngrid)%ground_fliq              (    i,j,ip)   )

               end select
               !---------------------------------------------------------------------------!
            end do tloop
            !------------------------------------------------------------------------------!
         end do patloop1
         !---------------------------------------------------------------------------------!

      end do lonloop
   end do latloop
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !     Normalise accumulated fluxes and albedo seen by atmosphere over one full BRAMS    !
   ! timestep (dtlt).                                                                      !
   !---------------------------------------------------------------------------------------!
   call normal_accfluxes(mxp,myp,npatch,ia,iz,ja,jz,scratch%vt2di,leaf_g(ngrid)%patch_area &
                        ,turb_g(ngrid)%sflux_u,turb_g(ngrid)%sflux_v,turb_g(ngrid)%sflux_w &
                        ,turb_g(ngrid)%sflux_t,turb_g(ngrid)%sflux_r,turb_g(ngrid)%sflux_c &
                        ,radiate_g(ngrid)%albedt,radiate_g(ngrid)%rlongup )
   !---------------------------------------------------------------------------------------!






   !---- Call TOPMODEL if the user wants so -----------------------------------------------!
   if (isfcl == 2) then
      call hydro(mxp,myp,nzg,nzs,npatch            , leaf_g(ngrid)%soil_water              &
                , leaf_g(ngrid)%soil_energy        , leaf_g(ngrid)%soil_text               &
                , leaf_g(ngrid)%sfcwater_mass      , leaf_g(ngrid)%sfcwater_energy         &
                , leaf_g(ngrid)%patch_area         , leaf_g(ngrid)%patch_wetind            )
   end if
   !---------------------------------------------------------------------------------------!


   !----- Update the 10-day running average for phenology. --------------------------------!
   call update_psibar(mxp,myp,nzg,npatch,ia,iz,ja,jz,dtlt                                  &
                     , leaf_g(ngrid)%soil_energy        , leaf_g(ngrid)%soil_water         &
                     , leaf_g(ngrid)%soil_text          , leaf_g(ngrid)%leaf_class         &
                     , leaf_g(ngrid)%psibar_10d         )
   !---------------------------------------------------------------------------------------!


   !----- Apply lateral boundary conditions to leaf3 arrays -------------------------------!
   call leaf3_bcond(mxp,myp,nzg,nzs,npatch,ia,iz,ja,jz,jdim,ibcon                          &
                      , leaf_g(ngrid)%soil_water         , leaf_g(ngrid)%sfcwater_mass     &
                      , leaf_g(ngrid)%soil_energy        , leaf_g(ngrid)%sfcwater_energy   &
                      , leaf_g(ngrid)%soil_color         , leaf_g(ngrid)%soil_text         &
                      , leaf_g(ngrid)%psibar_10d         , leaf_g(ngrid)%sfcwater_depth    &
                      , leaf_g(ngrid)%ustar              , leaf_g(ngrid)%tstar             &
                      , leaf_g(ngrid)%rstar              , leaf_g(ngrid)%cstar             &
                      , leaf_g(ngrid)%zeta               , leaf_g(ngrid)%ribulk            &
                      , leaf_g(ngrid)%veg_albedo         , leaf_g(ngrid)%veg_fracarea      &
                      , leaf_g(ngrid)%veg_lai            , leaf_g(ngrid)%veg_tai           &
                      , leaf_g(ngrid)%veg_rough          , leaf_g(ngrid)%veg_height        &
                      , leaf_g(ngrid)%veg_displace       , leaf_g(ngrid)%patch_area        &
                      , leaf_g(ngrid)%patch_rough        , leaf_g(ngrid)%patch_wetind      &   
                      , leaf_g(ngrid)%leaf_class         , leaf_g(ngrid)%soil_rough        &
                      , leaf_g(ngrid)%sfcwater_nlev      , leaf_g(ngrid)%stom_condct       &
                      , leaf_g(ngrid)%ground_rsat        , leaf_g(ngrid)%ground_rvap       &
                      , leaf_g(ngrid)%ground_temp        , leaf_g(ngrid)%ground_fliq       &
                      , leaf_g(ngrid)%veg_water          , leaf_g(ngrid)%veg_hcap          &
                      , leaf_g(ngrid)%veg_energy         , leaf_g(ngrid)%can_prss          &
                      , leaf_g(ngrid)%can_theiv          , leaf_g(ngrid)%can_theta         &
                      , leaf_g(ngrid)%can_rvap           , leaf_g(ngrid)%can_co2           &
                      , leaf_g(ngrid)%sensible_gc        , leaf_g(ngrid)%sensible_vc       &
                      , leaf_g(ngrid)%evap_gc            , leaf_g(ngrid)%evap_vc           &
                      , leaf_g(ngrid)%transp             , leaf_g(ngrid)%gpp               &
                      , leaf_g(ngrid)%plresp             , leaf_g(ngrid)%resphet           &
                      , leaf_g(ngrid)%veg_ndvip          , leaf_g(ngrid)%veg_ndvic         &
                      , leaf_g(ngrid)%veg_ndvif          , turb_g(ngrid)%sflux_u           &
                      , turb_g(ngrid)%sflux_v            , turb_g(ngrid)%sflux_w           &
                      , turb_g(ngrid)%sflux_t            , turb_g(ngrid)%sflux_r           &
                      , turb_g(ngrid)%sflux_c            , radiate_g(ngrid)%albedt         &
                      , radiate_g(ngrid)%rlongup         , leaf_g(ngrid)%rshort_gnd        &
                      , leaf_g(ngrid)%rlong_gnd          )
   !---------------------------------------------------------------------------------------!
   return
end subroutine leaf3_timestep
!==========================================================================================!
!==========================================================================================!
