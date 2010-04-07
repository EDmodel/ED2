!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine toptnest(ngra,ngrb)

  use mem_mksfc
  use mem_grid
  use io_params

  implicit none

  integer :: ngra,ngrb 

  integer :: ifm,icm,ipat,i,j,k,indfm,ivtime,nc1,mynum

  do ifm = ngra,ngrb
     icm = nxtnest(ifm)
     ! Initialize TOPOGRAPHY in toptinit.

     call toptinit(nnxp(ifm),nnyp(ifm),ifm  &
          ,sfcfile_p(ifm)%topt,sfcfile_p(ifm)%topzo)

     if (icm .ge. 1 .and. itoptflg(ifm) .eq. 0) then

        ! Interpolate TOPO from coarser grid:
        call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
             ,scr1,sfcfile_p(icm)%topt)
        call eintp(scr1,scr2,1,maxnxp,maxnyp  &
             ,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
        call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
             ,scr2,sfcfile_p(ifm)%topt)

        ! Interpolate TOPO ZO from coarser grid:
        call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
             ,scr1,sfcfile_p(icm)%topzo)
        call eintp(scr1,scr2,1,maxnxp,maxnyp  &
             ,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
        call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
             ,scr2,sfcfile_p(ifm)%topzo)

     elseif (itoptflg(ifm) .eq. 1) then

        ! Interpolate TOPO from standard dataset:
        call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%topt  &
             ,itoptfn(ifm),itoptfn(ifm),vt2da,vt2db,ifm,'TOP')
        ! Interpolate TOPO ZO from standard dataset:
        call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%topzo  &
             ,itoptfn(ifm),itoptfn(ifm),vt2da,vt2db,ifm,'ZOT')

     elseif (itoptflg(ifm) .eq. 3) then

        ! Interpolate TOPO from dted dataset:
        call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%topt  &
             ,itoptfn(ifm),itoptfn(ifm),vt2da,vt2db,ifm,'TOD')

        ! Interpolate TOPO ZO from dted dataset:
        call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%topzo  &
             ,itoptfn(ifm),itoptfn(ifm),vt2da,vt2db,ifm,'ZOD')

     endif

     ! If desired, override current values of TOPOGRAPHY in ruser.f subroutine.

     call toptinit_user(nnxp(ifm),nnyp(ifm),ifm  &
          ,sfcfile_p(ifm)%topt ,sfcfile_p(ifm)%topzo)

  enddo

  if (ngra .eq. ngrb) return

  ! In case topography data have been independently reassigned on any grid,
  ! average fine mesh topography sequentially to the coarser grids.

  do ifm = ngrb,ngra,-1
     if (nxtnest(ifm) .gt. ngridsh .and. ifm .ge. 2) then
        icm = nxtnest(ifm)

        call fdback(sfcfile_p(icm)%topt,sfcfile_p(ifm)%topt  &
             ,vt2da,scr2,1,nnxp(icm),nnyp(icm)  &
             ,1,nnxp(ifm),nnyp(ifm),ifm,'terr',vt2db)

     endif
  enddo

  ! In case terrain heights have been independently reassigned on
  ! any grid, interpolate coarse grid terrain heights to a temporary
  ! fine mesh array.  Fill the fine mesh boundary terrain heights
  ! from the temporary array.

  do ifm = ngra,ngrb
     icm = nxtnest(ifm)
     if (icm .ge. 1) then
        call fillscr(1,nxpmax,nypmax,1,nnxp(icm),nnyp(icm),1,1  &
             ,scr1,sfcfile_p(icm)%topt)

        call eintp(scr1,scr2,1,nxpmax,nypmax  &
             ,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
        call fillvar(1,nxpmax,nypmax,1,nnxp(ifm),nnyp(ifm),1,1  &
             ,scr2,scr1)

        nc1 = jdim * (nstraty(ifm) + 1)
        call ae2(nnxp(ifm),nnyp(ifm),2+nstratx(ifm)  &
             ,nnxp(ifm)-1-nstratx(ifm),1+nc1,nnyp(ifm)-nc1  &
             ,scr1,sfcfile_p(ifm)%topt)
        call ae1(nnxp(ifm)*nnyp(ifm),sfcfile_p(ifm)%topt,scr1)

     endif
  enddo
  return
end subroutine toptnest

!*************************************************************************

subroutine geonest_file(ifm)

  use mem_mksfc
  use mem_grid
  use io_params
  use mem_leaf

  implicit none

  integer :: ifm 

  integer :: icm,ipat,i,j,k,indfm,ivtime,nc1,mynum

  icm = nxtnest(ifm)

  ! Initialize PATCH AREA, LANDUSE CLASS, and SOIL TEXTURAL CLASS
  ! in subroutine sfcinit.

  call sfcinit_file(nnxp(ifm),nnyp(ifm),nzg,npatch,ifm      &
       ,sfcfile_p(ifm)%patch_area   &
       ,sfcfile_p(ifm)%leaf_class   &
       ,sfcfile_p(ifm)%soil_text)


  print*, ' '
  print*,'====================================================='
  print*,'Starting landuse data input on grid ',ifm
  print*,'====================================================='

  if (icm .ge. 1 .and. ivegtflg(ifm) .eq. 0) then

     ! Assign PATCH AREAS and PATCH CLASSES from coarser grid:

     do ipat = 1,npatch
        do j = 1,nnyp(ifm)
           do i = 1,nnxp(ifm)

              sfcfile_p(ifm)%patch_area(i,j,ipat) =  &
                   sfcfile_p(icm)%patch_area(ipm(i,ifm),jpm(j,ifm),ipat)
              sfcfile_p(ifm)%leaf_class(i,j,ipat) =  &
                   sfcfile_p(icm)%leaf_class(ipm(i,ifm),jpm(j,ifm),ipat)

           enddo
        enddo
     enddo


  elseif (ivegtflg(ifm) .eq. 1) then

     ! Assign PATCH AREAS and PATCH CLASSES from standard dataset:

     call landuse_opqr(nnxp(ifm),nnyp(ifm),nzg,npatch,nvegpat  &
          ,ivegtflg(ifm),ivegtfn(ifm),isoilflg(ifm),isoilfn(ifm) &
          ,ndviflg(ifm),ndvifn(ifm),vndvifil(1,ifm)  &
          ,'veg',platn(ifm),plonn(ifm)        &
          ,sfcfile_p(ifm)%soil_text  &
          ,sfcfile_p(ifm)%patch_area   &
          ,sfcfile_p(ifm)%leaf_class   &
          ,sfcfile_p(ifm)%veg_ndvif)


  endif

  if (icm .ge. 1 .and. isoilflg(ifm) .eq. 0) then

     ! Assign SOIL TEXTURE CLASS from coarser grid

     do ipat = 2,npatch
        do k = 1,nzg
           do j = 1,nnyp(ifm)
              do i = 1,nnxp(ifm)
                 sfcfile_p(ifm)%soil_text(k,i,j,ipat) =  &
                      sfcfile_p(icm)%soil_text(k,ipm(i,ifm),jpm(j,ifm),ipat)
              enddo
           enddo
        enddo
     enddo

  elseif (isoilflg(ifm) .eq. 1) then

     ! Assign SOIL TEXTURE CLASS from standard dataset:

     call landuse_opqr(nnxp(ifm),nnyp(ifm),nzg,npatch,nvegpat  &
          ,ivegtflg(ifm),ivegtfn(ifm),isoilflg(ifm),isoilfn(ifm) &
          ,ndviflg(ifm),ndvifn(ifm),vndvifil(1,ifm)  &
          ,'soil',platn(ifm),plonn(ifm)        &
          ,sfcfile_p(ifm)%soil_text  &
          ,sfcfile_p(ifm)%patch_area   &
          ,sfcfile_p(ifm)%leaf_class   &
          ,sfcfile_p(ifm)%veg_ndvif)

  endif

  ! If desired, override current values of PATCH AREA, PATCH CLASS, 
  ! LEAF-2 VEGETATION CLASS, SOIL TEXTURAL CLASS, and/or
  ! NDVI in ruser.f subroutines.

  call sfcinit_file_user(nnxp(ifm),nnyp(ifm),nzg,npatch,ifm &
       ,sfcfile_p(ifm)%patch_area  &
       ,sfcfile_p(ifm)%leaf_class  &
       ,sfcfile_p(ifm)%soil_text  )

  ! As a final initialization step, eliminate any land patch area that is less 
  ! than 1% of the total grid cell area.  Set its area to zero, and compensate
  ! by enlarging areas of remaining patches.

  call patch_minsize(nnxp(ifm),nnyp(ifm),npatch  &
       ,sfcfile_p(ifm)%patch_area)


  return
end subroutine geonest_file

!*************************************************************************

subroutine geonest_nofile(ngra,ngrb)

  use mem_leaf
  use mem_basic
  use mem_scratch
  use mem_grid
  use io_params

  ! For Soil Moisture Init.
  use mem_soil_moisture, only : SOIL_MOIST  ! INTENT(IN)

  implicit none

  integer :: ngra,ngrb

  integer :: isiz,ifm,icm,ipat,i,j,k,indfm,ivtime,nc1,mynum,ic,jc,can_shv,veg_fliq

  ! Initialization/interpolation of leaf-2 variables for which standard RAMS
  ! datasets never exist.

  isiz = maxnxp * maxnyp

  do ifm = ngra,ngrb
     icm = nxtnest(ifm)

     if (co2_on) then
        call atob(nnzp(ifm)*nnxp(ifm)*nnyp(ifm),basic_g(ifm)%co2p,scratch%vt3do)
     else
        call ae0(nnzp(ifm)*nnxp(ifm)*nnyp(ifm),scratch%vt3do,co2con(1))
     end if

     ! First, fill NOFILE LEAF-2 variables with default values in SFCINIT.

     call sfcinit_nofile(nnzp(ifm),nnxp(ifm),nnyp(ifm),nzg,nzs,npatch,ifm                  &
          , basic_g(ifm)%theta                    , basic_g(ifm)%pi0                       &
          , basic_g(ifm)%pp                       , basic_g(ifm)%rv                        &
          , scratch%vt3do                         , leaf_g(ifm)%seatp                      &
          , leaf_g(ifm)%seatf                     , leaf_g(ifm)%soil_water                 &
          , leaf_g(ifm)%soil_energy               , leaf_g(ifm)%soil_text                  &
          , leaf_g(ifm)%sfcwater_mass             , leaf_g(ifm)%sfcwater_energy            &
          , leaf_g(ifm)%sfcwater_depth            , leaf_g(ifm)%ustar                      &
          , leaf_g(ifm)%tstar                     , leaf_g(ifm)%rstar                      &
          , leaf_g(ifm)%cstar                     , leaf_g(ifm)%veg_fracarea               &
          , leaf_g(ifm)%veg_agb                   , leaf_g(ifm)%veg_lai                    &
          , leaf_g(ifm)%veg_tai                   , leaf_g(ifm)%veg_rough                  &
          , leaf_g(ifm)%veg_height                , leaf_g(ifm)%veg_albedo                 &
          , leaf_g(ifm)%patch_area                , leaf_g(ifm)%patch_rough                &
          , leaf_g(ifm)%patch_wetind              , leaf_g(ifm)%leaf_class                 &
          , leaf_g(ifm)%soil_rough                , leaf_g(ifm)%sfcwater_nlev              &
          , leaf_g(ifm)%stom_resist               , leaf_g(ifm)%ground_rsat                &
          , leaf_g(ifm)%ground_rvap               , leaf_g(ifm)%ground_temp                &
          , leaf_g(ifm)%ground_fliq               , leaf_g(ifm)%veg_water                  &
          , leaf_g(ifm)%veg_hcap                  , leaf_g(ifm)%veg_energy                 &
          , leaf_g(ifm)%can_prss                  , leaf_g(ifm)%can_theta                  &
          , leaf_g(ifm)%can_rvap                  , leaf_g(ifm)%can_co2                    &
          , leaf_g(ifm)%sensible                  , leaf_g(ifm)%evap                       &
          , leaf_g(ifm)%transp                    , leaf_g(ifm)%gpp                        &
          , leaf_g(ifm)%plresp                    , leaf_g(ifm)%resphet                    &
          , leaf_g(ifm)%veg_ndvip                 , leaf_g(ifm)%veg_ndvic                  &
          , leaf_g(ifm)%veg_ndvif                 , leaf_g(ifm)%snow_mass                  &
          , leaf_g(ifm)%snow_depth                , scratch%vt2da                          &
          , scratch%vt2db                         , scratch%vt2dc                          &
          , scratch%vt2dd                         , scratch%vt2de                          &
          , grid_g(ifm)%glat                      , grid_g(ifm)%glon                       &
          , grid_g(ifm)%topzo                     , grid_g(ifm)%flpw                       &
          , grid_g(ifm)%rtgt                      )

     ! Assignment section for NOFILE leaf-2 variables

     if (icm > 0) then
        if( nofilflg(ifm) == 0) then
           ! Assign values from coarse grid cells and patches
           do ipat = 1,npatch
              do j = 1,nnyp(ifm)
                 do i = 1,nnxp(ifm)
                    ic = ipm(i,ifm)
                    jc = jpm(j,ifm) 

                    do k = 1,nzg
                       leaf_g(ifm)%soil_water           (k,i,j,ipat) = &
                            leaf_g(icm)%soil_water      (k,ic,jc,ipat)
                       leaf_g(ifm)%soil_energy          (k,i,j,ipat) = &
                            leaf_g(icm)%soil_energy     (k,ic,jc,ipat)
                    enddo

                    do k = 1,nzs
                       leaf_g(ifm)%sfcwater_mass        (k,i,j,ipat) = &
                            leaf_g(icm)%sfcwater_mass   (k,ic,jc,ipat)
                       leaf_g(ifm)%sfcwater_energy      (k,i,j,ipat) = &
                            leaf_g(icm)%sfcwater_energy (k,ic,jc,ipat)  
                       leaf_g(ifm)%sfcwater_depth       (k,i,j,ipat) = &
                            leaf_g(icm)%sfcwater_depth  (k,ic,jc,ipat)
                    enddo

                    leaf_g(ifm)%veg_fracarea         (i,j,ipat) = &
                         leaf_g(icm)%veg_fracarea    (ic,jc,ipat)
                    leaf_g(ifm)%veg_agb              (i,j,ipat) = &
                         leaf_g(icm)%veg_agb         (ic,jc,ipat)
                    leaf_g(ifm)%veg_lai              (i,j,ipat) = &
                         leaf_g(icm)%veg_lai         (ic,jc,ipat)
                    leaf_g(ifm)%veg_tai              (i,j,ipat) = &
                         leaf_g(icm)%veg_tai         (ic,jc,ipat)
                    leaf_g(ifm)%veg_rough            (i,j,ipat) = &
                         leaf_g(icm)%veg_rough       (ic,jc,ipat)
                    leaf_g(ifm)%veg_height           (i,j,ipat) = &
                         leaf_g(icm)%veg_height      (ic,jc,ipat)
                    leaf_g(ifm)%veg_albedo           (i,j,ipat) = &
                         leaf_g(icm)%veg_albedo      (ic,jc,ipat)
                    leaf_g(ifm)%patch_rough          (i,j,ipat) = &
                         leaf_g(icm)%patch_rough     (ic,jc,ipat)
                    leaf_g(ifm)%patch_wetind         (i,j,ipat) = &
                         leaf_g(icm)%patch_wetind    (ic,jc,ipat)
                    leaf_g(ifm)%soil_rough           (i,j,ipat) = &
                         leaf_g(icm)%soil_rough      (ic,jc,ipat)
                    leaf_g(ifm)%sfcwater_nlev        (i,j,ipat) = &
                         leaf_g(icm)%sfcwater_nlev   (ic,jc,ipat)
                    leaf_g(ifm)%stom_resist          (i,j,ipat) = &
                         leaf_g(icm)%stom_resist     (ic,jc,ipat) 
                    leaf_g(ifm)%ground_rsat          (i,j,ipat) = &
                         leaf_g(icm)%ground_rsat     (ic,jc,ipat)
                    leaf_g(ifm)%ground_rvap          (i,j,ipat) = &
                         leaf_g(icm)%ground_rvap     (ic,jc,ipat)
                    leaf_g(ifm)%ground_temp          (i,j,ipat) = &
                         leaf_g(icm)%ground_temp     (ic,jc,ipat)
                    leaf_g(ifm)%ground_fliq          (i,j,ipat) = &
                         leaf_g(icm)%ground_fliq     (ic,jc,ipat)
                    leaf_g(ifm)%veg_water            (i,j,ipat) = &
                         leaf_g(icm)%veg_water       (ic,jc,ipat)
                    leaf_g(ifm)%veg_energy           (i,j,ipat) = &
                         leaf_g(icm)%veg_energy      (ic,jc,ipat) 
                    leaf_g(ifm)%veg_hcap             (i,j,ipat) = &
                         leaf_g(icm)%veg_hcap        (ic,jc,ipat) 
                    leaf_g(ifm)%can_rvap             (i,j,ipat) = &
                         leaf_g(icm)%can_rvap        (ic,jc,ipat)
                    leaf_g(ifm)%can_co2              (i,j,ipat) = &
                         leaf_g(icm)%can_co2         (ic,jc,ipat)
                    leaf_g(ifm)%can_theta            (i,j,ipat) = &
                         leaf_g(icm)%can_theta       (ic,jc,ipat) 
                    leaf_g(ifm)%can_prss             (i,j,ipat) = &
                         leaf_g(icm)%can_prss        (ic,jc,ipat) 
                    leaf_g(ifm)%veg_ndvic            (i,j,ipat) = &
                         leaf_g(icm)%veg_ndvic       (ic,jc,ipat)
                    leaf_g(ifm)%sensible             (i,j,ipat) = &
                         leaf_g(icm)%sensible        (ic,jc,ipat)
                    leaf_g(ifm)%evap                 (i,j,ipat) = &
                         leaf_g(icm)%evap            (ic,jc,ipat)
                    leaf_g(ifm)%transp               (i,j,ipat) = &
                         leaf_g(icm)%transp          (ic,jc,ipat)
                    leaf_g(ifm)%gpp                  (i,j,ipat) = &
                         leaf_g(icm)%gpp             (ic,jc,ipat)
                    leaf_g(ifm)%plresp               (i,j,ipat) = &
                         leaf_g(icm)%plresp          (ic,jc,ipat)
                    leaf_g(ifm)%resphet              (i,j,ipat) = &
                         leaf_g(icm)%resphet         (ic,jc,ipat)

                 end do
              end do
           end do

        elseif(nofilflg(ifm) == 1) then

           ! Interpolate from coarse grid. 
           ! We can interpolate water patch directly.
           !   For land patches, do this by first averaging all
           !   coarse grid land patches, interpolate, then assign back to all
           !   fine grid land patches.
           call patch_interp_driver(icm,ifm)

        end if

     end if

     ! Heterogeneous Soil Moisture Initialization
     if ((SOIL_MOIST == 'i').or.(SOIL_MOIST == 'I').or.  &
         (SOIL_MOIST == 'a').or.(SOIL_MOIST == 'A')) then
        call soil_moisture_init(nnzp(ifm),nnxp(ifm),nnyp(ifm),nzg,npatch,ifm               &
                             ,leaf_g(ifm)%can_theta         ,leaf_g(ifm)%can_prss          &
                             ,grid_g(ifm)%glat              ,grid_g(ifm)%glon              &
                             ,leaf_g(ifm)%soil_water        ,leaf_g(ifm)%soil_energy       &
                             ,leaf_g(ifm)%soil_text         )
     end if

     ! Override any of the above variable assignments by user-specified changes
     ! to subroutine sfcinit_nofile_user.

     call sfcinit_nofile_user(nnzp(ifm),nnxp(ifm),nnyp(ifm),nzg,nzs,npatch,ifm             &
               , basic_g(ifm)%theta                  , basic_g(ifm)%pi0                    &
               , basic_g(ifm)%pp                     , basic_g(ifm)%rv                     &
               , scratch%vt3do                       , leaf_g(ifm)%soil_water              &
               , leaf_g(ifm)%soil_energy             , leaf_g(ifm)%soil_text               &
               , leaf_g(ifm)%sfcwater_mass           , leaf_g(ifm)%sfcwater_energy         &
               , leaf_g(ifm)%sfcwater_depth          , leaf_g(ifm)%ustar                   &
               , leaf_g(ifm)%tstar                   , leaf_g(ifm)%rstar                   &
               , leaf_g(ifm)%cstar                   , leaf_g(ifm)%veg_fracarea            &
               , leaf_g(ifm)%veg_agb                 , leaf_g(ifm)%veg_lai                 &
               , leaf_g(ifm)%veg_tai                 , leaf_g(ifm)%veg_rough               &
               , leaf_g(ifm)%veg_height              , leaf_g(ifm)%veg_albedo              &
               , leaf_g(ifm)%patch_area              , leaf_g(ifm)%patch_rough             &
               , leaf_g(ifm)%patch_wetind            , leaf_g(ifm)%leaf_class              &
               , leaf_g(ifm)%soil_rough              , leaf_g(ifm)%sfcwater_nlev           &
               , leaf_g(ifm)%stom_resist             , leaf_g(ifm)%ground_rsat             &
               , leaf_g(ifm)%ground_rvap             , leaf_g(ifm)%ground_temp             &
               , leaf_g(ifm)%ground_fliq             , leaf_g(ifm)%veg_water               &
               , leaf_g(ifm)%veg_hcap                , leaf_g(ifm)%veg_energy              &
               , leaf_g(ifm)%can_prss                , leaf_g(ifm)%can_theta               &
               , leaf_g(ifm)%can_rvap                , leaf_g(ifm)%can_co2                 &
               , leaf_g(ifm)%sensible                , leaf_g(ifm)%evap                    &
               , leaf_g(ifm)%transp                  , leaf_g(ifm)%gpp                     &
               , leaf_g(ifm)%plresp                  , leaf_g(ifm)%resphet                 &
               , leaf_g(ifm)%veg_ndvip               , leaf_g(ifm)%veg_ndvic               &
               , leaf_g(ifm)%veg_ndvif               , leaf_g(ifm)%snow_mass               &
               , leaf_g(ifm)%snow_depth              , scratch%vt2da                       &
               , scratch%vt2db                       , scratch%vt2dc                       &
               , scratch%vt2dd                       , scratch%vt2de                       &
               , grid_g(ifm)%glat                    , grid_g(ifm)%glon                    &
               , grid_g(ifm)%topzo                   , grid_g(ifm)%flpw                    &
               , grid_g(ifm)%rtgt                    )
  end do

  return
end subroutine geonest_nofile




!******************************************************************************

subroutine patch_interp(icm,ifm,nc1,nc2,nc3,nc4,nf1,nf2,nf3,nf4 &
     ,ac,af,pareac,pareaf,avgc,avgf,slabc,slabf)

  use mem_scratch

  implicit none

  integer :: icm,ifm,nc1,nc2,nc3,nc4,nf1,nf2,nf3,nf4
  real :: ac(nc1,nc2,nc3,nc4), af(nf1,nf2,nf3,nf4)
  real :: pareac(nc2,nc3,nc4), pareaf(nf2,nf3,nf4)
  real :: avgc(nc1,nc2,nc3),   avgf(nf1,nf2,nf3)
  real :: slabc(nc2,nc3), slabf(nf2,nf3)

  integer :: k,i,j

  ! Average coarse grid field over all land patches

  call patch_land_average(nc1,nc2,nc3,nc4,pareac,ac,avgc)

  ! Interpolate patch-averaged to fine grid

  do k=1,nc1    ! nc1 and nf1 are the same

     do j=1,nc3
        do i=1,nc2
           slabc(i,j)=avgc(k,i,j)
        enddo
     enddo

     call fmint2d(icm,ifm,'t',slabc,slabf)

     do j=1,nf3
        do i=1,nf2
           avgf(k,i,j)=slabf(i,j)
        enddo
     enddo

  enddo

  ! Fill fine grid field back into all land patches

  call patch_land_unaverage(nf1,nf2,nf3,nf4,avgf,af)

  return
end subroutine patch_interp


!     *****************************************************************

subroutine fmint5(var1,var2,dn0xc,dn0xf,vt2da,ifm,icm,vpnt,idwt)

  use mem_scratch
  use mem_grid

  implicit none
  integer :: ifm,icm,idwt

  real, dimension(*) :: var1,var2,vt2da,dn0xc,dn0xf
  character(len=*) :: vpnt

  if (icm .eq. 0) return

  call fillscr(maxnzp,maxnxp,maxnyp,nnzp(icm),nnxp(icm),nnyp(icm)  &
       ,1,nnzp(icm),scratch%scr1,var1)


  if (idwt .eq. 1) then
     call dnswt2(maxnzp,maxnxp,maxnyp,nnzp(icm),nnxp(icm),nnyp(icm)  &
          ,scratch%scr1,dn0xc,vpnt,1)
  endif

  call eintp(scratch%scr1,scratch%scr2,maxnzp,maxnxp,maxnyp  &
       ,nnzp(ifm),nnxp(ifm),nnyp(ifm),ifm,3,vpnt,0,0)

  call fillvar(maxnzp,maxnxp,maxnyp,nnzp(ifm),nnxp(ifm),nnyp(ifm)  &
       ,1,nnzp(ifm),scratch%scr2,var2)

  if (idwt .eq. 1) call dnswt2(nnzp(ifm),nnxp(ifm),nnyp(ifm),nnzp(ifm)  &
       ,nnxp(ifm),nnyp(ifm),var2,dn0xf,vpnt,2)

  call rtgintrp(nnzp(ifm),nnxp(ifm),nnyp(ifm),var2,vt2da  &
       ,grid_g(ifm)%topt,ifm,vpnt)

  return
end subroutine fmint5


!******************************************************************************

subroutine patch_minsize(n2,n3,npat,patch_area)

  implicit none
  integer :: n2,n3,npat,i,j,ipat,jpat

  real :: orig_size
  real, dimension(n2,n3,npat) :: patch_area

  do j = 1,n3
     do i = 1,n2
        do ipat = 2,npat
           if (patch_area(i,j,ipat) > 0. .and. patch_area(i,j,ipat) < .01) then
              orig_size = patch_area(i,j,ipat)
              patch_area(i,j,ipat) = 0.
              do jpat = 1,npat
                 if (jpat /= ipat) then
                    patch_area(i,j,jpat) = patch_area(i,j,jpat) / (1. - orig_size)
                 endif
              end do
           end if
        end do
     end do
  end do

  return
end subroutine patch_minsize

! TEB
!#############################################################################
subroutine fusonest(ngra,ngrb)
!#############################################################################

  use mem_mksfc
  use mem_grid
  use io_params

  implicit none

  integer :: ngra,ngrb 

  integer :: ifm,icm,ipat,i,j,k,indfm,ivtime,nc1,mynum



  do ifm = ngra,ngrb
     icm = nxtnest(ifm)
     ! Initialize FUSO in fusoinit.
     
!     write(*,*)'glon =',grid_g(ifm)%glon (1,1)
!     call fusoinit(nnxp(ifm),nnyp(ifm),ifm  &
!          ,sfcfile_p(ifm)%fuso(1,1),grid_g(ifm)%glon (1,1))

     if (icm .ge. 1 .and. ifusflg(ifm) .eq. 0) then

        ! Interpolate FUSO from coarser grid:
        call fillscr(1,maxnxp,maxnyp,1,nnxp(icm),nnyp(icm),1,1  &
             ,scr1,sfcfile_p(icm)%fuso)
        call eintp(scr1,scr2,1,maxnxp,maxnyp  &
             ,1,nnxp(ifm),nnyp(ifm),ifm,2,'t',0,0)
        call fillvar(1,maxnxp,maxnyp,1,nnxp(ifm),nnyp(ifm),1,1  &
             ,scr2,sfcfile_p(ifm)%fuso)

     elseif (ifusflg(ifm) .eq. 1) then

        ! Interpolate FUSO from standard dataset:
        call geodat(nnxp(ifm),nnyp(ifm),sfcfile_p(ifm)%fuso  &
             ,ifusfn(ifm),ifusfn(ifm),vt2da,vt2db,ifm,'FUS')

     endif


  enddo

  if (ngra .eq. ngrb) return

  return
end subroutine fusonest

!*************************************************************************
subroutine patch_interp_driver(icm,ifm)
   use mem_leaf
   use mem_basic
   use mem_scratch
   use mem_grid
!  use io_params
   implicit none
   integer, intent(in) :: icm, ifm

   call patch_interp(icm,ifm,nzg,nnxp(icm),nnyp(icm),npatch,nzg,nnxp(ifm),nnyp(ifm),npatch &
                    ,leaf_g(icm)%soil_water,leaf_g(ifm)%soil_water,leaf_g(icm)%patch_area  &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db)
   call patch_interp(icm,ifm,nzg,nnxp(icm),nnyp(icm),npatch,nzg,nnxp(ifm),nnyp(ifm),npatch &
                    ,leaf_g(icm)%soil_energy,leaf_g(ifm)%soil_energy                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area,scratch%vt3da           &
                    ,scratch%vt3db,scratch%vt2da,scratch%vt2db)
   call patch_interp(icm,ifm,nzs,nnxp(icm),nnyp(icm),npatch,nzs,nnxp(ifm),nnyp(ifm),npatch &
                    ,leaf_g(icm)%sfcwater_mass,leaf_g(ifm)%sfcwater_mass                   &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area,scratch%vt3da           &
                    ,scratch%vt3db,scratch%vt2da,scratch%vt2db)
   call patch_interp(icm,ifm,nzs,nnxp(icm),nnyp(icm),npatch,nzs,nnxp(ifm),nnyp(ifm),npatch &
                    ,leaf_g(icm)%sfcwater_energy,leaf_g(ifm)%sfcwater_energy               &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,nzs,nnxp(icm),nnyp(icm),npatch,nzs,nnxp(ifm),nnyp(ifm),npatch &
                    ,leaf_g(icm)%sfcwater_depth,leaf_g(ifm)%sfcwater_depth                 &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )

   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_fracarea,leaf_g(ifm)%veg_fracarea                     &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_agb,leaf_g(ifm)%veg_agb                               &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_lai,leaf_g(ifm)%veg_lai                               &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_tai,leaf_g(ifm)%veg_tai                               &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_rough,leaf_g(ifm)%veg_rough                           &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_height,leaf_g(ifm)%veg_height                         &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_albedo,leaf_g(ifm)%veg_albedo                         &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%patch_rough,leaf_g(ifm)%patch_rough                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_fracarea,leaf_g(ifm)%veg_fracarea                     &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%patch_wetind,leaf_g(ifm)%patch_wetind                     &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%soil_rough,leaf_g(ifm)%soil_rough                         &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db)
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%sfcwater_nlev,leaf_g(ifm)%sfcwater_nlev                   &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%stom_resist,leaf_g(ifm)%stom_resist                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db ) 
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%ground_rsat,leaf_g(ifm)%ground_rsat                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%ground_rvap,leaf_g(ifm)%ground_rvap                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%ground_temp,leaf_g(ifm)%ground_temp                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%ground_fliq,leaf_g(ifm)%ground_fliq                       &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_water,leaf_g(ifm)%veg_water                           &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_hcap,leaf_g(ifm)%veg_hcap                             &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_energy,leaf_g(ifm)%veg_energy                         &
                    ,leaf_g(icm)%patch_area,leaf_g(icm)%patch_area                         &
                    ,scratch%vt3da,scratch%vt3db,scratch%vt2da,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%can_prss,leaf_g(ifm)%can_prss,leaf_g(icm)%patch_area      &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db ) 
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%can_rvap,leaf_g(ifm)%can_rvap,leaf_g(icm)%patch_area      &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db ) 
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%can_co2,leaf_g(ifm)%can_co2,leaf_g(icm)%patch_area        &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db ) 
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%can_theta,leaf_g(ifm)%can_theta,leaf_g(icm)%patch_area    &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db )
   call patch_interp(icm,ifm,1,nnxp(icm),nnyp(icm),npatch,1,nnxp(ifm),nnyp(ifm),npatch     &
                    ,leaf_g(icm)%veg_ndvic,leaf_g(ifm)%veg_ndvic,leaf_g(icm)%patch_area    &
                    ,leaf_g(icm)%patch_area,scratch%vt3da,scratch%vt3db,scratch%vt2da      &
                    ,scratch%vt2db )

   return
end subroutine patch_interp_driver
