!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine grid_setup(num)

  use mem_grid, only: &
       if_adap,       &
       npatch,        &
       ngrids,        &
       nnxp,          &
       nnyp,          &
       nnzp,          &
       nzg,           &
       nnstbot,       &
       dztn,          &
       xmn,           &
       ymn,           &
       zmn,           &
       platn,         &
       plonn,         &
       grid_g,        &
       dump_mem_grid

  implicit none

  integer :: num

  integer :: ifm,ifileok

  if (num == 1) then

     call gridinit
     call gridset(1)
!!$     call dump_mem_grid()

  else  

     do ifm = 1,ngrids
        call newgrid(ifm)

        call polarst(nnxp(ifm),nnyp(ifm)              &
             ,grid_g(ifm)%glat   ,grid_g(ifm)%glon    &
             ,grid_g(ifm)%fmapu  ,grid_g(ifm)%fmapv   &
             ,grid_g(ifm)%fmapt  ,grid_g(ifm)%fmapm   &
             ,grid_g(ifm)%fmapui ,grid_g(ifm)%fmapvi  &
             ,grid_g(ifm)%fmapti ,grid_g(ifm)%fmapmi  )

        call grdspc(nnxp(ifm),nnyp(ifm)             &
             ,grid_g(ifm)%dxu   ,grid_g(ifm)%dxv    &
             ,grid_g(ifm)%dxt   ,grid_g(ifm)%dxm    &
             ,grid_g(ifm)%dyu   ,grid_g(ifm)%dyv    &
             ,grid_g(ifm)%dyt   ,grid_g(ifm)%dym    &
             ,grid_g(ifm)%fmapu ,grid_g(ifm)%fmapv  &
             ,grid_g(ifm)%fmapt ,grid_g(ifm)%fmapm  )

        ! Define transformation Jacobians for all grids

        call fill_toptuvm(nnxp(ifm),nnyp(ifm)       &
             ,grid_g(ifm)%topt  ,grid_g(ifm)%topu   &
             ,grid_g(ifm)%topv  ,grid_g(ifm)%topm   &
             ,grid_g(ifm)%topta ,grid_g(ifm)%topma  )

        call transfm(nnxp(ifm),nnyp(ifm)          &
             ,grid_g(ifm)%topt ,grid_g(ifm)%topu  &
             ,grid_g(ifm)%topv ,grid_g(ifm)%topm  &
             ,grid_g(ifm)%rtgt ,grid_g(ifm)%rtgu  &
             ,grid_g(ifm)%rtgv ,grid_g(ifm)%rtgm  &
             ,grid_g(ifm)%f13u ,grid_g(ifm)%f13v  &
             ,grid_g(ifm)%f13t ,grid_g(ifm)%f13m  &
             ,grid_g(ifm)%f23u ,grid_g(ifm)%f23v  &
             ,grid_g(ifm)%f23t ,grid_g(ifm)%f23m  &
             ,grid_g(ifm)%dxu  ,grid_g(ifm)%dxv   &
             ,grid_g(ifm)%dxt  ,grid_g(ifm)%dxm   &
             ,grid_g(ifm)%dyu  ,grid_g(ifm)%dyv   &
             ,grid_g(ifm)%dyt  ,grid_g(ifm)%dym   )

        call lpuvw_init(nnxp(ifm),nnyp(ifm),grid_g(ifm)%flpu  &
             ,grid_g(ifm)%flpv,grid_g(ifm)%flpw  )

        if (if_adap == 1) then
           call ctrlvols (nnzp(ifm),nnxp(ifm),nnyp(ifm),nnstbot(ifm)  &
                ,dztn(:,ifm),xmn(:,ifm),ymn(:,ifm),zmn(:,ifm)         &
                ,platn(ifm),plonn(ifm)                                &
                ,grid_g(ifm)%aru               ,grid_g(ifm)%arv       &
                ,grid_g(ifm)%arw               ,grid_g(ifm)%volt      &
                ,grid_g(ifm)%volu              ,grid_g(ifm)%volv      &
                ,grid_g(ifm)%volw              ,grid_g(ifm)%flpu      &
                ,grid_g(ifm)%flpv              ,grid_g(ifm)%flpw      &
                ,grid_g(ifm)%dxu               ,grid_g(ifm)%dxv       &
                ,grid_g(ifm)%dxt               ,grid_g(ifm)%dyu       &
                ,grid_g(ifm)%dyv               ,grid_g(ifm)%dyt       &
                ,grid_g(ifm)%topma             ,grid_g(ifm)%topm      &
                ,ifm, nzg, npatch  )
        end if

     enddo
  endif
  return
end subroutine grid_setup





subroutine gridinit

  ! gridinit: set up domain sizes for all grids
  !   the user has specified nnxp, nnyp, nnzp, nzg, nzs, and npatch.
  !   fill other arrays that are denote sizes and are a function of these.

  use mem_grid, only: &
       ngrids,        & ! intent(in)
       nnxp,          & ! intent(in)
       nnyp,          & ! intent(in)
       nnzp,          & ! intent(in)
       nzg,           & ! intent(in)
       nzs,           & ! intent(in)
       npatch,        & ! intent(in)
       ngrid,         & ! intent(out)
       jdim,          & ! intent(out)
       nnx,           & ! intent(out)
       nnx1,          & ! intent(out)
       nnx2,          & ! intent(out)
       nny,           & ! intent(out)
       nny1,          & ! intent(out)
       nny2,          & ! intent(out)
       nnz,           & ! intent(out)
       nnz1,          & ! intent(out)
       nnxyzp,        & ! intent(out)
       nnxysp,        & ! intent(out)
       nnxyp            ! intent(out)

  implicit none

  ! 1D or 2D horizontal simulation

  jdim = 1
  if (nnyp(1) == 1) jdim = 0

  ! set up x, y and z direction # points for all grids

  do ngrid = 1,ngrids
     nnx(ngrid) = nnxp(ngrid) - 1
     nnx1(ngrid) = nnxp(ngrid) - 2
     nnx2(ngrid) = nnxp(ngrid) - 3
     if (jdim == 1) then
        nny(ngrid) = nnyp(ngrid) - 1
        nny1(ngrid) = nnyp(ngrid) - 2
        nny2(ngrid) = nnyp(ngrid) - 3
     else
        nny(ngrid) = 1
        nny1(ngrid) = 1
        nny2(ngrid) = 1
     endif
     nnz(ngrid) = nnzp(ngrid) - 1
     nnz1(ngrid) = nnzp(ngrid) - 2
  enddo

  ! surface and volume # points

  do ngrid = 1,ngrids
     nnxyzp(ngrid) = nnxp(ngrid) * nnyp(ngrid) * nnzp(ngrid)
     nnxysp(ngrid) = nnxp(ngrid) * nnyp(ngrid) * (nzg+nzs+3) * npatch
     nnxyp(ngrid) = nnxp(ngrid) * nnyp(ngrid)
  enddo
  return
end subroutine gridinit



subroutine polarst(n2,n3,glat,glon,fmapu,fmapv,fmapt,fmapm  &
     ,fmapui,fmapvi,fmapti,fmapmi)

  use mem_grid
  use rconstants

  implicit none

  integer :: n2,n3
  real, dimension(n2,n3) :: glat,glon,fmapu,fmapv,fmapt,fmapm  &
       ,fmapui,fmapvi,fmapti,fmapmi

  integer :: i,j
  real :: c1,xm2,xt2,ym2,yt2
  !  Calculates map factors and inverse map factors at u,v,t,m-points and
  !  geographical lat/lon at t-points for a given polar stereographic grid

  c1 = (2. * erad) ** 2
  do j = 1,n3
     do i = 1,n2
        xm2 = xm(i) * xm(i)
        xt2 = xt(i) * xt(i)
        ym2 = ym(j) * ym(j)
        yt2 = yt(j) * yt(j)

        fmapt(i,j) = 1. + (xt2 + yt2) / c1
        fmapu(i,j) = 1. + (xm2 + yt2) / c1
        fmapv(i,j) = 1. + (xt2 + ym2) / c1
        fmapm(i,j) = 1. + (xm2 + ym2) / c1

        fmapui(i,j) = 1.0 / fmapu(i,j)
        fmapvi(i,j) = 1.0 / fmapv(i,j)
        fmapti(i,j) = 1.0 / fmapt(i,j)
        fmapmi(i,j) = 1.0 / fmapm(i,j)

        call xy_ll(glat(i,j),glon(i,j),platn(ngrid),plonn(ngrid)  &
             ,xt(i),yt(j))

        !       write(6,344)i,j,fmapt(i,j),fmapm(i,j),glat(i,j),glon(i,j)
        ! 344   format('polst:i,j,fmt,fmm,glt,gln',2i4,4e12.3)

     enddo
  enddo

  if (ihtran == 0) then
     call ae0(n2*n3,fmapu,1.)
     call ae0(n2*n3,fmapv,1.)
     call ae0(n2*n3,fmapt,1.)
     call ae0(n2*n3,fmapm,1.)
     call ae0(n2*n3,fmapui,1.)
     call ae0(n2*n3,fmapvi,1.)
     call ae0(n2*n3,fmapti,1.)
     call ae0(n2*n3,fmapmi,1.)
  endif
  return
end subroutine polarst





subroutine grdspc(n2,n3,dxu,dxv,dxt,dxm,dyu,dyv,dyt,dym  &
     ,fmapu,fmapv,fmapt,fmapm)
  use mem_grid
  implicit none

  integer :: n2,n3
  real, dimension(n2,n3) :: dxu,dxv,dxt,dxm,dyu,dyv,dyt,dym  &
       ,fmapu,fmapv,fmapt,fmapm

  integer :: i,j

  do j = 1,n3
     do i = 1,n2-1
        dxu(i,j) = fmapu(i,j) / (xtn(i+1,ngrid)-xtn(i,ngrid))
        dxm(i,j) = fmapm(i,j) / (xtn(i+1,ngrid)-xtn(i,ngrid))
     enddo
     dxu(nxp,j)=dxu(nx,j)*fmapu(nxp,j)/fmapu(nx,j)
     dxm(nxp,j)=dxm(nx,j)*fmapm(nxp,j)/fmapm(nx,j)
     do i = 2,n2
        dxv(i,j)=fmapv(i,j)/(xmn(i,ngrid)-xmn(i-1,ngrid))
        dxt(i,j)=fmapt(i,j)/(xmn(i,ngrid)-xmn(i-1,ngrid))
     enddo
     dxv(1,j)=dxv(2,j)*fmapv(1,j)/fmapv(2,j)
     dxt(1,j)=dxt(2,j)*fmapt(1,j)/fmapt(2,j)
  enddo

  if (jdim == 1) then
     do i = 1,n2
        do j = 1,n3-1
           dyv(i,j)=fmapv(i,j)/(ytn(j+1,ngrid)-ytn(j,ngrid))
           dym(i,j)=fmapm(i,j)/(ytn(j+1,ngrid)-ytn(j,ngrid))
        enddo
        dyv(i,nyp)=dyv(i,ny)*fmapv(i,nyp)/fmapv(i,ny)
        dym(i,nyp)=dym(i,ny)*fmapm(i,nyp)/fmapm(i,ny)
        do j = 2,n3
           dyu(i,j)=fmapu(i,j)/(ymn(j,ngrid)-ymn(j-1,ngrid))
           dyt(i,j)=fmapt(i,j)/(ymn(j,ngrid)-ymn(j-1,ngrid))
        enddo
        dyu(i,1)=dyu(i,2)*fmapu(i,1)/fmapu(i,2)
        dyt(i,1)=dyt(i,2)*fmapt(i,1)/fmapt(i,2)
     enddo
  else
     do i=1,n2
        do j=1,n3
           dyu(i,j)=1./deltayn(ngrid)
           dyv(i,j)=1./deltayn(ngrid)
           dyt(i,j)=1./deltayn(ngrid)
           dym(i,j)=1./deltayn(ngrid)
        enddo
     enddo
  endif
  return
end subroutine grdspc






subroutine fill_toptuvm(n2,n3,topt,topu,topv,topm,topta,topma)
  use mem_grid
  implicit none

  integer :: n2,n3
  real, dimension(n2,n3) :: topt,topu,topv,topm,topta,topma

  integer :: i,j
  real :: terdev

  do j = 1,n3
     do i = 1,n2
        topt(i,j) = topta(i,j)
     enddo
  enddo

  terdev = 0.
  do j = 1,n3
     do i = 1,n2-1
        topu(i,j) = topt(i,j) + (topt(i+1,j) - topt(i,j))  &
             * (xm(i) - xt(i)) / (xt(i+1) - xt(i))
        terdev = max(terdev,abs(topt(i,j)))
     enddo
     topu(n2,j) = topt(n2,j) + (topt(n2,j) - topt(n2-1,j))  &
          * (xm(n2) - xt(n2)) / (xt(n2) - xt(n2-1))
  enddo

  if (terdev < 1.e-6) then
     itopo = 0
  else
     itopo = 1
  endif

  if (jdim == 1) then
     do i = 1,n2
        do j = 1,n3-1
           topv(i,j) = topt(i,j) + (topt(i,j+1) - topt(i,j))  &
                * (ym(j) - yt(j)) / (yt(j+1) - yt(j))
           topm(i,j) = topu(i,j) + (topu(i,j+1) - topu(i,j))  &
                * (ym(j) - yt(j)) / (yt(j+1) - yt(j))
        enddo
        topv(i,n3) = topt(i,n3) + (topt(i,n3) - topt(i,n3-1))  &
             * (ym(n3) - yt(n3)) / (yt(n3) - yt(n3-1))
        topm(i,n3) = topu(i,n3) + (topu(i,n3) - topu(i,n3-1))  &
             * (ym(n3) - yt(n3)) / (yt(n3) - yt(n3-1))
     enddo
  else
     do j = 1,n3
        do i = 1,n2
           topv(i,j) = topt(i,j)
           topm(i,j) = topu(i,j)
        enddo
     enddo
  endif

  do j = 1,n3
     do i = 1,n2
        topma(i,j) = topm(i,j)
     enddo
  enddo

  if (if_adap == 1) then
     do j = 1,n3
        do i = 1,n2
           topt(i,j) = 0.
           topu(i,j) = 0.
           topv(i,j) = 0.
           topm(i,j) = 0.
        enddo
     enddo
  endif
  return
end subroutine fill_toptuvm





subroutine transfm(n2,n3,topt,topu,topv,topm,rtgt,rtgu,rtgv,rtgm  &
     ,f13u,f13v,f13t,f13m,f23u,f23v,f23t,f23m  &
     ,dxu,dxv,dxt,dxm,dyu,dyv,dyt,dym)

  use mem_grid

  implicit none

  integer :: n2,n3
  real, dimension(n2,n3) :: topt,topu,topv,topm,rtgt,rtgu,rtgv,rtgm  &
       ,f13u,f13v,f13t,f13m,f23u,f23v,f23t,f23m  &
       ,dxu,dxv,dxt,dxm,dyu,dyv,dyt,dym

  integer :: iztflag=0,i,j,k

  !     this routine computes the coordinate transformation constants
  !     based on the topographical values of TOPT.

  ztop = zmn(nnzp(1)-1,1)
  do k = 1,nzp
     htn(k,ngrid) = zt(k) / ztop - 1.
     hwn(k,ngrid) = zm(k) / ztop - 1.
  enddo
  do k = 1,nzp
     ht2n(k,ngrid) = .5 * htn(k,ngrid)
     ht4n(k,ngrid) = .25 * htn(k,ngrid)
     hw2n(k,ngrid) = .5 * hwn(k,ngrid)
     hw4n(k,ngrid) = .25 * hwn(k,ngrid)
  enddo
  do k = 1,nzp
     ht(k)  = htn(k,ngrid)
     hw(k)  = hwn(k,ngrid)
     ht2(k) = ht2n(k,ngrid)
     ht4(k) = ht4n(k,ngrid)
     hw2(k) = hw2n(k,ngrid)
     hw4(k) = hw4n(k,ngrid)
  enddo

  do j = 1,n3
     do i = 1,n2
        rtgt(i,j) = 1. - topt(i,j) / ztop
        rtgu(i,j) = 1. - topu(i,j) / ztop
        rtgv(i,j) = 1. - topv(i,j) / ztop
        rtgm(i,j) = 1. - topm(i,j) / ztop
        if (topt(i,j) > .5 * ztop) then
           print*, 'Terrain height is over half the model domain'
           print*, 'height.  Model will stop here to avoid this.'
           print*, 'ngrid, i, j, topt, ztop = ',ngrid,i,j,topt(i,j),ztop
           iztflag = 1
        endif
     enddo
  enddo

  if (iztflag == 1) stop 'topt/ztop'

  do j = 1,n3
     do i = 2,n2
        f13t(i,j) = (topu(i,j) - topu(i-1,j)) * dxt(i,j) / rtgt(i,j)
        f13v(i,j) = (topm(i,j) - topm(i-1,j)) * dxv(i,j) / rtgv(i,j)
     enddo
     do i = 1,n2-1
        f13u(i,j) = (topt(i+1,j) - topt(i,j)) * dxu(i,j) / rtgu(i,j)
        f13m(i,j) = (topv(i+1,j) - topv(i,j)) * dxm(i,j) / rtgm(i,j)
     enddo
     f13t(1,j)  = f13u(1,j)
     f13v(1,j)  = f13m(1,j)
     f13u(n2,j) = f13t(n2,j)
     f13m(n2,j) = f13v(n2,j)
  enddo

  do i = 1,n2
     do j = 2,n3
        f23t(i,j) = (topv(i,j) - topv(i,j-jdim)) * dyt(i,j) / rtgt(i,j)
        f23u(i,j) = (topm(i,j) - topm(i,j-jdim)) * dyu(i,j) / rtgu(i,j)
     enddo
     do j = 1,n3-1
        f23v(i,j) = (topt(i,j+jdim) - topt(i,j)) * dyv(i,j) / rtgv(i,j)
        f23m(i,j) = (topu(i,j+jdim) - topu(i,j)) * dym(i,j) / rtgm(i,j)
     enddo
     if (jdim == 1) then
        f23t(i,1)  = f23v(i,1)
        f23u(i,1)  = f23m(i,1)
        f23v(i,n3) = f23t(i,n3)
        f23m(i,n3) = f23u(i,n3)
     endif
  enddo
  return
end subroutine transfm



subroutine newgrid(ngr)

  use mem_grid
  use node_mod

  implicit none
  integer :: ngr

  integer :: i,j,k

  !     +----------------------------------------------------------------
  !     !    Fill the single and 1D variables that the rest of the model
  !     !      uses from the nest arrays and change grid level in the I/O.
  !     +----------------------------------------------------------------



  ngrid = ngr

  !         grid point references

  !         x - direction

  nxp=nnxp(ngr)
  nx=nnx(ngr)
  nx1=nnx1(ngr)
  nx2=nnx2(ngr)
  !
  !         y - direction
  !
  nyp=nnyp(ngr)
  ny=nny(ngr)
  ny1=nny1(ngr)
  ny2=nny2(ngr)
  !
  !         z - direction
  !
  nzp=nnzp(ngr)
  nzpp=nzp+1
  nz=nnz(ngr)
  nz1=nnz1(ngr)


  nxyzp=nnxyzp(ngr)
  nxysp=nnxysp(ngr)
  nxyp=nnxyp(ngr)

  !          grid spacings

  deltax=deltaxn(ngr)
  do i=1,nxp
     xt(i)=xtn(i,ngr)
     xm(i)=xmn(i,ngr)
  enddo
  !
  deltay=deltayn(ngr)
  do j=1,nyp
     yt(j)=ytn(j,ngr)
     ym(j)=ymn(j,ngr)
  enddo
  !
  deltaz=zmn(2,ngr)-zmn(1,ngr)
  do k=1,nzp
     zt(k)=ztn(k,ngr)
     zm(k)=zmn(k,ngr)
     dzm(k)=dzmn(k,ngr)
     dzt(k)=dztn(k,ngr)
     dzm2(k)=dzm2n(k,ngr)
     dzt2(k)=dzt2n(k,ngr)
     ht(k)=htn(k,ngr)
     ht2(k)=ht2n(k,ngr)
     ht4(k)=ht4n(k,ngr)
     hw(k)=hwn(k,ngr)
     hw2(k)=hw2n(k,ngr)
     hw4(k)=hw4n(k,ngr)
  enddo

  !         option flags

  nsttop=nnsttop(ngr)
  nstbot=nnstbot(ngr)

  !         timesteps

  dtlt=dtlongn(ngr)
  dtlv=2.*dtlt

  !        node gridpoint info

  mxp=mmxp(ngr)
  myp=mmyp(ngr)
  mzp=mmzp(ngr)
  ia=mia(ngr)
  iz=miz(ngr)
  ja=mja(ngr)
  jz=mjz(ngr)
  i0=mi0(ngr)
  j0=mj0(ngr)
  ibcon=mibcon(ngr)
  return
end subroutine newgrid
