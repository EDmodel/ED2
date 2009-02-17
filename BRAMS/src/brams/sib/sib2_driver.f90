!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine sfclyr_sib(mzp,mxp,myp,ia,iz,ja,jz,ibcon)

  use mem_all
  use therm_lib, only: level

  implicit none

  integer :: mzp,mxp,myp,ia,iz,ja,jz,ibcon

  integer :: ng
  !integer, save :: ncall=0

  if (nstbot == 0) return

  !print*,'ncall=',ncall
  !if(ncall == 0) then
  !   ncall=1
  !   open(11,file='leaf.lis', status='unknown')
  !   rewind 11

  !print*,'calling alloc=',ncall
  !   call alloc_leafcol(nzg,nzs)

  !endif


  ng=ngrid

  call sib_driver(mzp,mxp,myp,nzg,nzs,npatch,ia,iz,ja,jz           &
       ,leaf_g (ng), basic_g (ng), turb_g (ng), radiate_g(ng)      &
       ,grid_g (ng), cuparm_g(ng), micro_g(ng)                     &
       ,scratch%vt2da  ,scratch%vt2db  ,scratch%vt2dc              &
       ,scratch%vt2dd  ,scratch%vt2de  ,scratch%vt2df              &
       ,scratch%vt3da                                              &
       ,ng)

  ! Apply lateral boundary conditions to leaf3/SiB arrays

  call leaf_bcond(mxp,myp,nzg,nzs,npatch,jdim                        &
       ,leaf_g(ng)%soil_water           ,leaf_g(ng)%sfcwater_mass    &
       ,leaf_g(ng)%soil_energy          ,leaf_g(ng)%sfcwater_energy  &
       ,leaf_g(ng)%soil_text            ,leaf_g(ng)%sfcwater_depth   &
       ,leaf_g(ng)%ustar                ,leaf_g(ng)%tstar            &
       ,leaf_g(ng)%rstar                ,leaf_g(ng)%veg_albedo       &
       ,leaf_g(ng)%veg_fracarea         ,leaf_g(ng)%veg_lai          &
       ,leaf_g(ng)%veg_tai                                           &
       ,leaf_g(ng)%veg_rough            ,leaf_g(ng)%veg_height       &
       ,leaf_g(ng)%patch_area           ,leaf_g(ng)%patch_rough      &
       ,leaf_g(ng)%patch_wetind         ,leaf_g(ng)%leaf_class       &
       ,leaf_g(ng)%soil_rough           ,leaf_g(ng)%sfcwater_nlev    &
       ,leaf_g(ng)%stom_resist          ,leaf_g(ng)%ground_rsat      &
       ,leaf_g(ng)%ground_rvap          ,leaf_g(ng)%veg_water        &
       ,leaf_g(ng)%veg_temp             ,leaf_g(ng)%can_rvap         &
       ,leaf_g(ng)%can_temp             ,leaf_g(ng)%veg_ndvip        &
       ,leaf_g(ng)%veg_ndvic            ,leaf_g(ng)%veg_ndvif        )

  return
end subroutine sfclyr_sib

!*****************************************************************************

subroutine sib_driver(m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz   &
     ,leaf,basic,turb,radiate,grid,cuparm,micro    &
     ,ths2,rvs2,pis2,dens2,ups2,vps2,zts2, ng      )

  use mem_all
  use leaf_coms
  use rconstants

  use mem_sib_co2   ! For SiB

  use mem_sib, only : sib_brams_g
  ! TEB_SPM
  use teb_spm_start, only: TEB_SPM !INTENT(IN)
  use therm_lib, only: rslif,qwtk,qtk,level

  implicit none

  integer :: m1,m2,m3,mzg,mzs,np,ia,iz,ja,jz

  integer :: ng

  type (leaf_vars)    leaf
  type (basic_vars)   basic
  type (turb_vars)    turb
  type (radiate_vars) radiate
  type (grid_vars)    grid
  type (cuparm_vars)  cuparm
  type (micro_vars)   micro

  real, dimension(m2,m3) :: ths2,rvs2,pis2,dens2,ups2,vps2,zts2

  integer :: i,j,ip,iter_leaf


  !itb...stuff for SiB...
  integer, external :: julday
  integer :: doy

  real :: dvelu,dvelv,velnew,sflux_uv,cosine1,sine1

  !srf - SIB2 - temporary diagnostic/output SIB2 variables
  real :: cupr,lspr

!!$  real :: co2flx,pco2c,assimn,FSS,FWS

  real :: co2flx,pco2c,pco2m_sib

  ! Local Variables
  ! Needed by TEB_SPM
  real :: L_EMIS_TOWN, L_ALB_TOWN, L_TS_TOWN, L_G_URBAN

  integer :: ksn,nsoil,k
  integer, save :: i_init_stars_sib
  data i_init_stars_sib/0/
  

  ! TEB_SPM
  

  !srf - SIB2

  ! Time interpolation factor for updating SST

  if (iupdsst == 0) then
     timefac_sst = 0.
  else
    timefac_sst = real((time-ssttime1(ngrid))/(ssttime2(ngrid)-ssttime1(ngrid)))
  endif

  ! Define SIB timestep as the same as dtlong (dtll)
  niter_leaf = 1
  niter_can  = 1

  dtll_factor = 1. / float(niter_leaf)
  dtll        = dtlt * dtll_factor
  dtlc_factor = 1. / float(niter_can)
  dtlc        = dtll * dtlc_factor

  dtllohcc = dtll / hcapcan
  dtllowcc = dtll / wcapcan
  dtlcohcc = dtlc / hcapcan
  dtlcowcc = dtlc / wcapcan
  dtlcohcv = dtlc / hcapveg


  ! Copy surface atmospheric variables into 2d arrays for input to leaf

  if (if_adap == 1) then
     call sfc_fields_adap(m1,m2,m3,ia,iz,ja,jz,jdim               &
          ,grid%flpu          ,grid%flpv        ,grid%flpw        &
          ,grid%topma         ,grid%aru         ,grid%arv         &
          ,basic%theta        ,basic%rv         ,basic%up         &
          ,basic%vp           ,basic%dn0        ,basic%pp         &
          ,basic%pi0          ,zt,zm,dzt                          &
          ,ths2,rvs2,ups2,vps2,pis2,dens2,zts2                    )
  else
     call sfc_fields(m1,m2,m3,ia,iz,ja,jz,jdim             &
          ,basic%theta        ,basic%rv         ,basic%up  &
          ,basic%vp           ,basic%dn0        ,basic%pp  &
          ,basic%pi0          ,grid%rtgt        ,zt        &
          ,ths2,rvs2,ups2,vps2,pis2,dens2,zts2             )
  endif

  do j = ja,jz
     do i = ia,iz

        ! Copy surface variables to single-column values

        ups = ups2(i,j)
        vps = vps2(i,j)
        ths = ths2(i,j)
        rvs = rvs2(i,j)
        zts = zts2(i,j)
        pis = pis2(i,j)
        dens = dens2(i,j)

        prss = pis ** cpor * p00
        vels = sqrt(ups ** 2 + vps ** 2)
        gzotheta = g * zts / ths

        ! Update water internal energy from time-dependent SST

        leaf%soil_energy(mzg,i,j,1) = alli  &
             + cliq * (leaf%seatp(i,j) + (leaf%seatf(i,j) - leaf%seatp(i,j)) &
             * timefac_sst - tsupercool)

        ! Fill surface precipitation arrays for input to SiB

        call sfc_pcp_sib(nnqparm(ngrid),level,i,j,cuparm,micro,cupr,lspr)

        !Zero out albedo,upward surface longwave,and momentum,heat,and moisture
        !flux arrays before summing over patches

        if (ilwrtyp > 0 .or. iswrtyp > 0) then
           radiate%albedt(i,j) = 0.
           radiate%rlongup(i,j) = 0.
        endif

        turb%sflux_u(i,j) = 0.
        turb%sflux_v(i,j) = 0.
        turb%sflux_w(i,j) = 0.
        turb%sflux_t(i,j) = 0.
        turb%sflux_r(i,j) = 0.

        ! Begin patch loop

        do ip = 1,np

           ! Update time-dependent vegetation LAI and fractional coverage

           if (ip >= 2 .and. leaf%patch_area(i,j,ip) >= .009) then

              if (ip >= 2) call vegndvi(ngrid                           &
                   ,leaf%patch_area  (i,j,ip) ,leaf%leaf_class(i,j,ip)  &
                   ,leaf%veg_fracarea(i,j,ip) ,leaf%veg_lai   (i,j,ip)  &
                   ,leaf%veg_tai     (i,j,ip) ,leaf%veg_rough (i,j,ip)  &
                   ,leaf%veg_height  (i,j,ip) ,leaf%veg_albedo(i,j,ip)  &
                   ,leaf%veg_ndvip   (i,j,ip) ,leaf%veg_ndvic (i,j,ip)  &
                   ,leaf%veg_ndvif   (i,j,ip)                           )

           endif

           ! Begin leaf small timestep here.

           do iter_leaf = 1,niter_leaf

              ! Calculate radiative fluxes between atmosphere, vegetation, and
              ! ground/snow based on already-computed downward shortwave and
              ! longwave fluxes from the atmosphere.  Fill tempk array with
              ! soil and snow temperature (C) and fracliq array with liquid
              ! fraction of water content in soil and snow.
              ! Other snowcover properties are also computed here.

              if (iswrtyp > 0 .or. ilwrtyp > 0) then

                 if (ip == 1 .or. leaf%patch_area(i,j,ip) >= .009) then
                    ! TEB_SPM - This is just to avoid problems with pointers in 
                    !           subroutines. 
                    if (TEB_SPM==1) then
                       L_G_URBAN   = leaf%G_URBAN(i,j,ip)
                    else
                       L_G_URBAN   = 0. 
                    end if
                    L_EMIS_TOWN = 0.
                    L_ALB_TOWN  = 0.
                    L_TS_TOWN   = 0.

                    call sfcrad(mzg,mzs,ip                 &
                         ,leaf%soil_energy     (:,i,j,ip)  &
                         ,leaf%soil_water      (:,i,j,ip)  &
                         ,leaf%soil_text       (:,i,j,ip)  &
                         ,leaf%sfcwater_energy (:,i,j,ip)  &
                         ,leaf%sfcwater_mass   (:,i,j,ip)  &
                         ,leaf%sfcwater_depth  (:,i,j,ip)  &
                         ,leaf%patch_area      (i,j,ip)    &
                         ,leaf%can_temp        (i,j,ip)    &
                         ,leaf%veg_temp        (i,j,ip)    &
                         ,leaf%leaf_class      (i,j,ip)    &
                         ,leaf%veg_height      (i,j,ip)    &
                         ,leaf%veg_fracarea    (i,j,ip)    &
                         ,leaf%veg_albedo      (i,j,ip)    &   
                         ,leaf%sfcwater_nlev   (i,j,ip)    &
                         ,radiate%rshort       (i,j)       &
                         ,radiate%rlong        (i,j)       &
                         ,radiate%albedt       (i,j)       &
                         ,radiate%rlongup      (i,j)       &
                         ,radiate%cosz         (i,j)       &
                         ! TEB_SPM
                         ,L_G_URBAN                                              &
                         ,L_EMIS_TOWN               ,L_ALB_TOWN                  &
                         ,L_TS_TOWN                                              &
                         !
                         )

                 endif

              endif

              ! For water surface (patch 1), compute surface saturation mixing
              ! ratio and roughness length based on previous ustar.
              ! For soil patches, compute roughness length based on vegetation
              ! and snow.

              if (ip == 1) then

                 leaf%ground_rsat(i,j,ip) = rslif(prss,tempk(mzg))   
                 leaf%patch_rough(i,j,ip)  &
                      = max(z0fac_water * leaf%ustar(i,j,ip) ** 2,.0001)

              else

                 if (leaf%patch_area(i,j,ip) >= .009) then
                    leaf%patch_rough(i,j,ip)   &
                         = max(grid%topzo(i,j),leaf%soil_rough(i,j,ip)  &
                         ,leaf%veg_rough(i,j,ip)) * (1. - snowfac)   &
                         + snowrough * snowfac

                 endif

              endif

              ! Calculate turbulent fluxes between atmosphere and canopy
              ! (or "canopy")

              if (leaf%patch_area(i,j,ip) >= .009) then
                 thetacan = leaf%can_temp(i,j,ip) / pis
              endif


              !srf   Use stars routine inside SIB2 code for ip > 1
              if (ip > 1 ) then
                 !
                 !srf - but, at first time initalize stars parameter avoiding
                 ! NAN for the turbulent fluxes. 
                 ! After that, use stars from SIB parameterization
                 if(i_init_stars_sib == 0 ) then
                    if(i==iz .and. j==jz .and. ip==np .and. ngrid==ngrids)   &
                         i_init_stars_sib = 1
                    call stars(leaf%ustar(i,j,ip),leaf%tstar(i,j,ip)         &
                         ,leaf%rstar(i,j,ip),ths,rvs,thetacan                &
                         ,leaf%can_rvap(i,j,ip)                              &
                         ,zts,leaf%patch_rough(i,j,ip)                       &
                         ,leaf%patch_area(i,j,ip)                            &
                         ,vels,vels_pat,vonk,dtllohcc,dens,dtll              &
	                 ,leaf%R_aer(i,j,ip))
                 endif
              else

                 call stars(leaf%ustar(i,j,ip),leaf%tstar(i,j,ip)            &
                      ,leaf%rstar(i,j,ip),ths,rvs,thetacan                   &
                      ,leaf%can_rvap(i,j,ip)                                 &
                      ,zts,leaf%patch_rough(i,j,ip),leaf%patch_area(i,j,ip)  &
                      ,vels,vels_pat,vonk,dtllohcc,dens,dtll                 &
	              ,leaf%R_aer(i,j,ip))
                 if (TEB_SPM == 1) then
                   L_G_URBAN = leaf%G_URBAN(i,j,ip)
                 end if
                 call sfclmcv(leaf%ustar(i,j,ip),leaf%tstar(i,j,ip)       &
                      ,leaf%rstar(i,j,ip),vels,vels_pat,ups,vps,gzotheta  &
                      ,leaf%patch_area(i,j,ip),turb%sflux_u(i,j)          &
                      ,turb%sflux_v(i,j),turb%sflux_w(i,j)                &
                      ,turb%sflux_t(i,j),turb%sflux_r(i,j)                &
                      ,L_G_URBAN)
              endif

              ! For water patches, update temperature and moisture of "canopy"
              ! from divergence of fluxes with water surface and atmosphere.
              ! rdi = ustar/5 is the viscous sublayer conductivity from
              ! Garratt (1992).

              if (ip == 1) then

                 rdi = .2 * leaf%ustar(i,j,1)

                 leaf%can_temp(i,j,1) = leaf%can_temp(i,j,1)          &
                      + dtllohcc * dens * cp                          &
                      * ((tempk(mzg) - leaf%can_temp(i,j,1)) * rdi    &
                      + leaf%ustar(i,j,1) * leaf%tstar(i,j,1) * pis)

                 leaf%can_rvap(i,j,1) = leaf%can_rvap(i,j,1) + dtllowcc*dens  &
                      *((leaf%ground_rsat(i,j,1) - leaf%can_rvap(i,j,1))*rdi  &
                      + leaf%ustar(i,j,1) * leaf%rstar(i,j,1))

              endif

              ! For soil model patches,update temperature and moisture of soil,
              ! vegetation, and canopy

              !srf -modification for include sib2
              !srf ------- SIB2 (isfcl =3)
              if (ip >= 2) then

                 if (leaf%patch_area(i,j,ip) >= .009) then
                    if (TEB_SPM == 1) then
                      L_G_URBAN = leaf%G_URBAN(i,j,ip)
                    end if

                    call sfclmcv(leaf%ustar(i,j,ip),leaf%tstar(i,j,ip)      &
                         ,leaf%rstar(i,j,ip),vels,vels_pat,ups,vps,gzotheta &
                         ,leaf%patch_area(i,j,ip),turb%sflux_u(i,j)         &
                         ,turb%sflux_v(i,j),turb%sflux_w(i,j)               &
                         ,turb%sflux_t(i,j),turb%sflux_r(i,j)               &
                         ,L_G_URBAN                                   &
                         !
                         )

                    !! Diagnose soil temperature and liquid fraction
                    !
                    do k = 1,mzg
                       nsoil = nint(leaf%soil_text(k,i,j,ip))
                       call qwtk(leaf%soil_energy(k,i,j,ip),  &
                            leaf%soil_water(k,i,j,ip)*wdns,   &
                            slcpd(nsoil), tempk(k), fracliq(k))
                    enddo
                    ! Diagnose snow temperature.

                    ksn = nint(leaf%sfcwater_nlev(i,j,ip))
                    do k = 1,ksn
                       call qtk(leaf%sfcwater_energy(k,i,j,ip),  &
                            tempk(k+mzg), fracliq(k+mzg))
                    enddo
                    if(ksn == 0) then
                       do k = 1,mzs
                          tempk(k+mzg) = tempk(mzg) !if there is not snow
                       enddo
                    endif

                    doy=julday(imontha,idatea,iyeara)

                    pco2m_sib = scalar_g(1, ng)%sclp(2, i, j) * (29./44.) *1.e6
!MLO - Commented out since the sib subroutine is full of bugs (I tried to add it back
!      with no luck :-/ ). Since we are developing ED, this is not a priority for us...
!                   call sib(mzg,mzs,dtll,ths,rvs,0.01*prss,dens,ths*pis     & 
!                        ,cupr,lspr,radiate%rshort(i,j),vels                 &
!                        ,radiate%rlong(i,j)                                 &
!                        ,zts,radiate%cosz(i,j),leaf%leaf_class(i,j,ip)      &
!                        ,leaf%veg_ndvic(i,j,ip),leaf%veg_ndvif(i,j,ip)      &
!                        ! respfactor_in ::: no futuro acrescente
!                        ,leaf%soil_text(mzg,i,j,ip)                         &
!                        !Soil type (k-dep)
!                        ,leaf%can_temp(i,j,ip),leaf%can_rvap(i,j,ip)        &
!                        ! canopy temp & water
!                        ,leaf%veg_temp(i,j,ip)                              &
!                        ! veg temp
!                        ! ,tempk(1+mzg),tempk                               &
!                        !soil temp (surf and deep layers)
!                        ,tempk                                              &
!                        !soil + snow temp
!                        ,leaf%soil_water(1,i,j,ip)                          &
!                        !soil moisture
!                        , sib_brams_g(ng)%snow1(i,j)                        &
!                        , sib_brams_g(ng)%snow2(i,j)                        &
!                        , sib_brams_g(ng)%capac1(i,j)                       &
!                        , sib_brams_g(ng)%capac2(i,j)                       &
!                        , sib_brams_g(ng)%rst(i,j)                          &
!                        !SiB prognostic variables
!                        !,FSS,FWS,co2flx                                     &
!                        !,FSS,FWS, sib_g(ng, 1)%SRC_CO2(i, j)                &
!                        , sib_g(ng, 1)%SRC_CO2(i, j)                        &
!                        ! Change "co2flx" by "sib_g(ng, 1)%SRC_CO2(i, j)"
!                        !sib FSS and FWS
!                        ,leaf%ustar(i,j,ip),leaf%rstar(i,j,ip)              &
!                        ,leaf%tstar(i,j,ip)                                 &
!                        , sib_brams_g(ng)%pco2ap(i,j)                       &
!!$                         ,pco2c,assimn,i,j                              &
!!$                         , scalar_g(1, ng)%sclp(2, i, j)                &
!!$                         , sib_brams_g(ng)%pco2m(i,j))
!                        , pco2c,i,j                          &
!                        , pco2m_sib                          &
!                        , scalar_g(1, ng)%sclp(2, i, j)       &
!                        , sib_brams_g(ng)%pco2m(i,j)         &
!                        , grid_g(ng)%glat(i,j),doy           &!latitude, time 
!                        , sib_brams_g(ng)%fss(i,j)           &!new diagnostics
!                        , sib_brams_g(ng)%fws(i,j)           &!new diagnostics
!                        , sib_brams_g(ng)%assimn(i,j)        &!new diagnostics
!                        , sib_brams_g(ng)%respg(i,j)         &!new diagnostics
!                        , sib_brams_g(ng)%rstfac4(i,j)       &!new diagnostics
!                        , sib_brams_g(ng)%rstfac1(i,j)       &!new diagnostics
!                        , sib_brams_g(ng)%rstfac2(i,j)       &!new diagnostics
!                        , sib_brams_g(ng)%rstfac3(i,j)       &!new diagnostics
!                        , sib_brams_g(ng)%ect(i,j)           &!new diagnostics
!                        , sib_brams_g(ng)%eci(i,j)           &!new diagnostics
!                        , sib_brams_g(ng)%egi(i,j)           &!new diagnostics
!                        , sib_brams_g(ng)%egs(i,j)           &!new diagnostics
!                        , sib_brams_g(ng)%hc(i,j)            &!new diagnostics
!                        , sib_brams_g(ng)%hg(i,j)            &
!                        , sib_brams_g(ng)%w1(i,j)            &
!                        , sib_brams_g(ng)%w2(i,j)            &
!                        , sib_brams_g(ng)%w3(i,j)            &
!                        , sib_brams_g(ng)%ww1(i,j)           &
!                        , sib_brams_g(ng)%ww2(i,j)           &
!                        , sib_brams_g(ng)%ww3(i,j)           &
!                        , sib_brams_g(ng)%exo(i,j)           &
!                        , sib_brams_g(ng)%ta(i,j)            &
!                        , sib_brams_g(ng)%tc(i,j)            &
!                        , sib_brams_g(ng)%tg(i,j)            &
!                        , sib_brams_g(ng)%td1(i,j)           &
!                        , sib_brams_g(ng)%td2(i,j)           &
!                        , sib_brams_g(ng)%td3(i,j)           &
!                        , sib_brams_g(ng)%td4(i,j)           &
!                        , sib_brams_g(ng)%td5(i,j)           &
!                        , sib_brams_g(ng)%td6(i,j)           &
!                        , sib_brams_g(ng)%ra(i,j)            &
!                        , sib_brams_g(ng)%rb(i,j)            &
!                        , sib_brams_g(ng)%rc(i,j)            &
!                        , sib_brams_g(ng)%rd(i,j)            &
!                        , sib_brams_g(ng)%roff(i,j)          &
!                        , sib_brams_g(ng)%zlt(i,j)           &
!                        , sib_brams_g(ng)%green(i,j)         &
!                        , sib_brams_g(ng)%apar(i,j)          &
!                        , sib_brams_g(ng)%nee(i,j)           &
!                        , sib_brams_g(ng)%cu(i,j)            &
!                        , sib_brams_g(ng)%ct(i,j)            &
!                        , sib_brams_g(ng)%ventmf(i,j)        &
!                        , sib_brams_g(ng)%pco2c(i,j)         &
!                        , sib_brams_g(ng)%pco2i(i,j)         &
!                        , sib_brams_g(ng)%pco2s(i,j)         &
!                        , sib_brams_g(ng)%ea(i,j)            &
!                        , sib_brams_g(ng)%sha(i,j)           &
!                        , sib_brams_g(ng)%em(i,j)            &
!                        , sib_brams_g(ng)%rha(i,j)           &
!                        , sib_brams_g(ng)%radvbc(i,j)        &
!                        , sib_brams_g(ng)%radvdc(i,j)        &
!                        , sib_brams_g(ng)%radnbc(i,j)        &
!                        , sib_brams_g(ng)%radndc(i,j)        &
!                        , sib_brams_g(ng)%dlwbot(i,j)        &
!                        , sib_brams_g(ng)%cp(i,j)            &
!                        , sib_brams_g(ng)%rho(i,j)           &
!                        , sib_brams_g(ng)%psy(i,j)           &
!                        , sib_brams_g(ng)%cupr(i,j)          &
!                        , sib_brams_g(ng)%lspr(i,j))          !new diagnostics


                    !calculate back soil energy 
                    do k=1,mzg
                       nsoil = nint(leaf%soil_text(k,i,j,ip))
                       leaf%soil_energy(k,i,j,ip) = tempk(k) * slcpd(nsoil)                &
                                                  + leaf%soil_water(k,i,j,ip) * cliqvlme   &
                                                  * (tempk(k) - tsupercool)
                    enddo

                 endif
              endif

           enddo
        enddo

     enddo
  enddo

  ! Normalize accumulated fluxes and albedo seen by atmosphere over model
  ! timestep dtlt.

  do j = ja,jz
     do i = ia,iz
        turb%sflux_u(i,j) = turb%sflux_u(i,j) * dtll_factor * dens2(i,j)
        turb%sflux_v(i,j) = turb%sflux_v(i,j) * dtll_factor * dens2(i,j)
        turb%sflux_w(i,j) = turb%sflux_w(i,j) * dtll_factor * dens2(i,j)
        turb%sflux_t(i,j) = turb%sflux_t(i,j) * dtll_factor * dens2(i,j)
        turb%sflux_r(i,j) = turb%sflux_r(i,j) * dtll_factor * dens2(i,j)
     enddo
  enddo

  if (ilwrtyp > 0 .or. iswrtyp > 0) then
     do j = ja,jz
        do i = ia,iz
           radiate%albedt (i,j) = radiate%albedt (i,j) * dtll_factor
           radiate%rlongup(i,j) = radiate%rlongup(i,j) * dtll_factor
        enddo
     enddo
  endif

  return
end subroutine sib_driver

!***************************************************************************


subroutine sfc_pcp_sib(nqparm,level,i,j,cuparm,micro,cupr,lspr)

  use mem_basic
  use mem_micro
  use mem_cuparm
  use leaf_coms

  implicit none

  integer :: nqparm,level,i,j,icld
  real :: cupr,lspr
  type (cuparm_vars)  cuparm
  type (micro_vars)   micro

  cupr = 0.
  lspr = 0.
  if (nqparm > 0) then
     do icld = 1,nclouds
        cupr = cupr + cuparm%conprr(i,j,icld) ! conv   precip at mm/s
     end do
  end if
  if (level >= 3) lspr = micro%   pcpg(i,j) ! explic precip at mm/s

!  if(cupr.gt. 0.) print*,'CUPR=',i,j,cupr


  return
end subroutine sfc_pcp_sib


