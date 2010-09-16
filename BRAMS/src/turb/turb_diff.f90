!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine diffvel(m1,m2,m3,ia,iz,ja,jz,jd  &
     ,ia_1,ja_1,ia1,ja1,iz_1,jz_1,iz1,jz1,izu,jzv,idiffkk  &
     ,up,vp,wp,ut,vt,wt,vt3da,vt3db,vt3dc  &
     ,vt3dd,vt3de,vt3df,vt3dg,vt3dj,vt3dk  &
     ,vt3dl,vt3dm,vt3dn,vt3do,rtgu,rtgv,rtgt  &
     ,sflux_u,sflux_v,sflux_w,dn0,dn0u,dn0v,scr1,scr2,ibcon,mynum)



  use mem_grid, only     : dzm,      &     !INTENT(IN)
                           nstbot,   &     !INTENT(IN)
                           nsttop,   &     !INTENT(IN)
                           dtlv,     &     !INTENT(IN)
                           dzt,      &     !INTENT(IN)
                           grid_g,   &     !INTENT(IN)
                           zt,       &     !INTENT(IN)
                           ngrid,    &     !INTENT(IN)
                           zm              !INTENT(IN)

  use mem_scratch, only  : vctr1,    &     !INTENT(OUT)
                           vctr2,    &     !INTENT(OUT)
                           vctr3,    &     !INTENT(OUT)
                           vctr4,    &     !INTENT(OUT)
                           vctr5,    &     !INTENT(OUT)
                           vctr6,    &     !INTENT(OUT)
                           vctr7           !INTENT(OUT)

  use mem_turb, only     : ihorgrad      & !INTENT(IN)
                         , ibotflx       ! !INTENT(IN)

  implicit none

  integer, INTENT(IN) :: m1       &
                       , m2       &
                       , m3       &
                       , ia       &
                       , iz       &
                       , ja       &
                       , jz       &
                       , jd       &
                       , ia_1     &
                       , ja_1     &
                       , ia1      &
                       , ja1      &
                       , iz_1     &
                       , jz_1     &
                       , iz1      &
                       , jz1      &
                       , izu      &
                       , jzv      &
                       , idiffkk  &
                       , ibcon    &
                       , mynum

  real, dimension(m1,m2,m3), INTENT(IN)    :: up    &
                                            , vp    &
                                            , wp    &
                                            , scr1  &
                                            , scr2  &
                                            , dn0   &
                                            , dn0u  &
                                            , dn0v

  real, dimension(m1,m2,m3), INTENT(INOUT) :: ut      &
                                            , vt      &
                                            , wt      &
                                            , vt3da   &
                                            , vt3db   &
                                            , vt3dc   &
                                            , vt3dd   &
                                            , vt3de   &
                                            , vt3df   &
                                            , vt3dg   &
                                            , vt3dj   &
                                            , vt3dk   &
                                            , vt3dl   &
                                            , vt3dm   &
                                            , vt3dn   &
                                            , vt3do

  real, dimension(m2,m3), INTENT(IN)  :: sflux_u    &
                                       , sflux_v    &
                                       , sflux_w    &
                                       , rtgt       &
                                       , rtgu       &
                                       , rtgv

  !local variables:
  integer :: i,j,k

  real :: akn,ako,akp,cross,c1,c2,dtlvi


  !          compute fluxes with k's and previous strains

  if (ihorgrad .eq. 1) then

     do j = 1,jz1
        do i = 1,iz1
           do k = 1,m1-1
              vt3da(k,i,j) = -vt3da(k,i,j) * scr2(k,i,j)
              vt3dc(k,i,j) = -vt3dc(k,i,j) * scr2(k,i,j)
           enddo
        enddo
     enddo

     do j = 1,jz
        do i = 1,iz
           do k = 1,m1-1
              akn = .25 * (scr2(k,i  ,j   ) + scr2(k+1,i  ,j   )  &
                   + scr2(k,i  ,j+jd) + scr2(k+1,i  ,j+jd))
              ako = .25 * (scr2(k,i  ,j   ) + scr2(k+1,i  ,j   )  &
                   + scr2(k,i+1,j   ) + scr2(k+1,i+1,j   ))
              akp = .25 * (scr2(k,i  ,j   ) + scr2(k  ,i+1,j   )  &
                   + scr2(k,i  ,j+jd) + scr2(k  ,i+1,j+jd))

              vt3db(k,i,j) = -vt3db(k,i,j) * akp
              vt3dn(k,i,j) = -vt3dn(k,i,j) * akp
              vt3dd(k,i,j) = -vt3dd(k,i,j) * ako
              vt3df(k,i,j) = -vt3df(k,i,j) * ako
              vt3de(k,i,j) = -vt3de(k,i,j) * akn
              vt3dg(k,i,j) = -vt3dg(k,i,j) * akn
           enddo
        enddo
     enddo

  elseif(ihorgrad.eq.2)then

     do j = 1,jz
        do i = 1,iz
           do k = 1,m1-1
              akn = .25 * (scr2(k,i  ,j ) + scr2(k+1,i  ,j )  &
                   + scr2(k,i  ,j+jd) + scr2(k+1,i  ,j+jd))
              ako = .25 * (scr2(k,i  ,j ) + scr2(k+1,i  ,j )  &
                   + scr2(k,i+1,j ) + scr2(k+1,i+1,j ))

              vt3df(k,i,j) = -vt3df(k,i,j) * ako
              vt3dg(k,i,j) = -vt3dg(k,i,j) * akn
           enddo
        enddo
     enddo

  endif

  cross = 0.
  select case (idiffkk)
  case (3:6)
     cross = 1.
  case default
     cross = 0.
  end select

  !------------------------ vertical u-component diffusion ----------------

  do j = ja,jz
     do i = ia,izu
        do k = 1,m1-1
           vctr1(k) = cross * vt3df(k,i,j)
           vctr2(k) = dzm(k) * (scr1(k,i,j) + scr1(k+1,i,j)  &
                + scr1(k,i+1,j) + scr1(k+1,i+1,j))
        enddo

        if (nstbot .eq. 1) then
           vctr1(1) = .5 * (sflux_u(i,j) + sflux_u(i+1,j))
        else
           vctr1(1) = vctr1(1) + up(1,i,j)  &
                * .25 * (scr1(1,i,j) + scr1(2,i,j) + scr1(1,i+1,j)  &
                + scr1(2,i+1,j)) * dzm(1) / rtgu(i,j)
        endif

        if (nsttop .eq. 1) then
           vctr1(m1-1) = 0.
        else
           vctr1(m1-1) = vctr1(m1-1) - up(m1-1,i,j)  &
                * .25 * (scr1(m1,i,j) + scr1(m1-1,i,j)  &
                + scr1(m1,i+1,j) + scr1(m1-1,i+1,j)) * dzm(m1-1)  &
                / rtgu(i,j)
        endif

        c1 = dtlv / rtgu(i,j)
        c2 = .25 * dtlv / (rtgu(i,j) * rtgu(i,j))
        do k = 2,m1-1
           vt3dl(k,i,j) = up(k,i,j) * dn0u(k,i,j)  &
                + c1 * dzt(k) * (vctr1(k-1) - vctr1(k))
           vt3dj(k,i,j) = - c2 * dzt(k) * vctr2(k-1)
           vt3dk(k,i,j) = - c2 * dzt(k) * vctr2(k)
           vt3do(k,i,j) = dn0u(k,i,j) - vt3dj(k,i,j) - vt3dk(k,i,j)
        enddo
        vt3dj(2,i,j) = 0.
        if (nstbot .eq. 1) vt3do(2,i,j) = dn0u(2,i,j) - vt3dk(2,i,j)
        vt3dk(m1-1,i,j) = 0.
        if (nsttop .eq. 1) vt3do(m1-1,i,j) = dn0u(m1-1,i,j)  &
             - vt3dj(m1-1,i,j)

     enddo
  enddo

  call tridiff1(m1,m2,m3,ia,izu,ja,jz,m1-1  &
       ,vt3dj,vt3do,vt3dk,vt3dl,vt3do,vt3dk)

  !---------------------horizontal u-component diffusion ----------------

  if(ihorgrad.eq.1)then
     call divcart(m1,m2,m3,ia,izu,ja,jz,vt3da,vt3dj,'XDIR','TPNT',ibotflx)
     call divcart(m1,m2,m3,ia,izu,ja,jz,vt3dn,vt3dk,'YDIR','PPNT',ibotflx)

  elseif(ihorgrad.eq.2)then
     !        average the dn0*hkh to the velocity points
     call avgvel(m1,m2,m3,ia,izu,ja,jz,'xdir',jd,vt3da,scr2)

     call truhor(m1,m2,m3,ia,izu,ja,jz                   &
                ,up,vt3dj,'xdir','dxt',grid_g(ngrid)%dxt &
                ,grid_g(ngrid)%topu,grid_g(ngrid)%rtgu   &
                ,zt,vctr1,vctr2,vctr3,vctr4,vctr5        &
                ,vctr6,vctr7,jd,vt3da,dn0,dtlv,ngrid)

     call avgvel(m1,m2,m3,ia,izu,ja,jz,'ydir',jd,vt3da,scr2)
     call truhor(m1,m2,m3,ia,izu,ja,jz                    &
                ,up,vt3dk,'ydir','dym',grid_g(ngrid)%dym  &
                ,grid_g(ngrid)%topu,grid_g(ngrid)%rtgu    &
                ,zt,vctr1,vctr2,vctr3,vctr4,vctr5         &
                ,vctr6,vctr7,jd,vt3da,dn0,dtlv,ngrid)

  endif

  dtlvi = 1.0 / dtlv

  if (jd .eq. 1) then
     do j = ja,jz
        do i = ia,izu
           do k = 2,m1-1
              ut(k,i,j) = ut(k,i,j)  &
                   + dtlvi * (vt3do(k,i,j) - up(k,i,j))  &
                   - (vt3dj(k,i,j) + vt3dk(k,i,j)) / dn0u(k,i,j)
           enddo
        enddo
     enddo
  else
     do j = ja,jz
        do i = ia,izu
           do k = 2,m1-1
              ut(k,i,j) = ut(k,i,j)  &
                   + dtlvi * (vt3do(k,i,j)-up(k,i,j))  &
                   - vt3dj(k,i,j) / dn0u(k,i,j)
           enddo
        enddo
     enddo
  endif

  !------------------------ vertical v-component diffusion ----------------

  if (jd .eq. 0) go to 99

  do j = ja,jzv
     do i = ia,iz
        do k = 1,m1-1
           vctr1(k) = cross * vt3dg(k,i,j)
           vctr2(k) = dzm(k) * (scr1(k,i,j) + scr1(k+1,i,j)  &
                + scr1(k,i,j+jd) + scr1(k+1,i,j+jd))
        enddo

        if (nstbot .eq. 1) then
           vctr1(1) = .5 * (sflux_v(i,j) + sflux_v(i,j+jd))
        else
           vctr1(1) = vctr1(1) + vp(1,i,j)  &
                * .25 * (scr1(1,i,j) + scr1(2,i,j) + scr1(1,i,j+jd)  &
                + scr1(2,i,j+jd)) * dzm(1) / rtgv(i,j)
        endif
        if (nsttop .eq. 1) then
           vctr1(m1-1) = 0.
        else
           vctr1(m1-1) = vctr1(m1-1) - vp(m1-1,i,j)  &
                * .25 * (scr1(m1,i,j) + scr1(m1-1,i,j)  &
                + scr1(m1,i,j+jd) + scr1(m1-1,i,j+jd)) * dzm(m1-1)  &
                / rtgv(i,j)
        endif

        c1 = dtlv / rtgv(i,j)
        c2 = .25 * dtlv / (rtgv(i,j) * rtgv(i,j))
        do k = 2,m1-1
           vt3dm(k,i,j) = vp(k,i,j) * dn0v(k,i,j)  &
                + c1 * dzt(k) * (vctr1(k-1) - vctr1(k))
           vt3dj(k,i,j) = -c2 * dzt(k) * vctr2(k-1)
           vt3dk(k,i,j) = -c2 * dzt(k) * vctr2(k)
           vt3do(k,i,j) = dn0v(k,i,j) - vt3dj(k,i,j) - vt3dk(k,i,j)
        enddo
        vt3dj(2,i,j) = 0.
        if (nstbot .eq. 1) vt3do(2,i,j) = dn0v(2,i,j) - vt3dk(2,i,j)

        vt3dk(m1-1,i,j) = 0.
        if (nsttop .eq. 1)  &
             vt3do(m1-1,i,j) = dn0v(m1-1,i,j) - vt3dj(m1-1,i,j)
     enddo
  enddo

  call tridiff1(m1,m2,m3,ia,iz,ja,jzv,m1-1  &
       ,vt3dj,vt3do,vt3dk,vt3dm,vt3do,vt3dk)

  !--------------------- horizontal v-component diffusion ----------------

  if (ihorgrad .eq. 1) then
     call divcart(m1,m2,m3,ia,iz,ja,jzv,vt3dc,vt3dj,'YDIR','TPNT',ibotflx)
     call divcart(m1,m2,m3,ia,iz,ja,jzv,vt3db,vt3dk,'XDIR','PPNT',ibotflx)

  elseif(ihorgrad.eq.2)then
     !        average the dn0*hkh to the velocity points
     call avgvel(m1,m2,m3,ia,iz,ja,jzv,'ydir',jd,vt3da,scr2)
     call truhor(m1,m2,m3,ia,iz,ja,jzv                    &
                ,vp,vt3dj,'ydir','dyt',grid_g(ngrid)%dyt  &
                ,grid_g(ngrid)%topv,grid_g(ngrid)%rtgv    &
                ,zt,vctr1,vctr2,vctr3,vctr4,vctr5         &
                ,vctr6,vctr7,jd,vt3da,dn0,dtlv,ngrid      )

     call avgvel(m1,m2,m3,ia,iz,ja,jzv,'xdir',jd,vt3da,scr2)
     call truhor(m1,m2,m3,ia,iz,ja,jzv  &
                ,vp,vt3dk,'xdir','dxm',grid_g(ngrid)%dxm  &
                ,grid_g(ngrid)%topv,grid_g(ngrid)%rtgv    &
                ,zt,vctr1,vctr2,vctr3,vctr4,vctr5         &
                ,vctr6,vctr7,jd,vt3da,dn0,dtlv,ngrid      )

  endif

  do j = ja,jzv
     do i = ia,iz
        do k = 2,m1-1
           vt(k,i,j) = vt(k,i,j) + dtlvi * (vt3do(k,i,j)-vp(k,i,j))  &
                - (vt3dj(k,i,j) + vt3dk(k,i,j)) / dn0v(k,i,j)
        enddo
     enddo
  enddo

99 continue

  !------------------------ vertical w-component diffusion ----------------

  do j = ja,jz
     do i = ia,iz
        do k = 2,m1-1
           vctr1(k) = 0.
           vctr2(k) = dzt(k) * scr1(k,i,j)
        enddo

        if (nstbot .eq. 1) then
           vctr1(2) = sflux_w(i,j)
        else
           vctr1(2) = (1.0 + cross) * wp(1,i,j)  &
                * scr1(2,i,j) * dzt(2) / rtgt(i,j)
        endif
        if (nsttop .eq. 1) then
           vctr1(m1-1) = 0.
        else
           vctr1(m1-1) = -(1.0 + cross) * wp(m1-1,i,j)  &
                * scr1(m1-1,i,j) * dzt(m1-1) / rtgt(i,j)
        endif

        c1 = dtlv / rtgt(i,j)
        c2 = 2.0 * dtlv / (rtgt(i,j) * rtgt(i,j))
        do k = 2,m1-2
           vt3dn(k,i,j) = wp(k,i,j)  &
                * (dn0(k,i,j) + dn0(k+1,i,j)) * .5  &
                + c1 * dzm(k) * (vctr1(k) - vctr1(k+1))
           vt3dj(k,i,j) = - c2 * dzm(k) * vctr2(k)
           vt3dk(k,i,j) = - c2 * dzm(k) * vctr2(k+1)
           vt3do(k,i,j) = (dn0(k,i,j) + dn0(k+1,i,j)) * .5  &
                - vt3dj(k,i,j) - vt3dk(k,i,j)
        enddo
        vt3dj(2,i,j) = 0.
        if (nstbot .eq. 1) vt3do(2,i,j) =  &
             (dn0(2,i,j) + dn0(3,i,j)) * .5 - vt3dk(2,i,j)
        vt3dk(m1-2,i,j) = 0.
        if (nsttop .eq. 1) vt3do(m1-2,i,j) =  &
             (dn0(m1-2,i,j) + dn0(m1-1,i,j)) * .5  &
             -vt3dj(m1-2,i,j)
     enddo
  enddo

  call tridiff1(m1,m2,m3,ia,iz,ja,jz,m1-2  &
       ,vt3dj,vt3do,vt3dk,vt3dn,vt3do,vt3dk)

  !--------------------- horizontal w-component diffusion ----------------

  select case (idiffkk)
  case (3:6)
     if (ihorgrad .eq. 1) then
        call divcart(m1,m2,m3,ia,iz,ja,jz,vt3dd,vt3dj,'XDIR','OPNT',ibotflx)
        call divcart(m1,m2,m3,ia,iz,ja,jz,vt3de,vt3dk,'YDIR','NPNT',ibotflx)

     elseif(ihorgrad.eq.2)then
        call truhor(m1,m2,m3,ia,iz,ja,jz                     &
                   ,wp,vt3dj,'xdir','dxu',grid_g(ngrid)%dxu  &
                   ,grid_g(ngrid)%topt,grid_g(ngrid)%rtgt    &
                   ,zm,vctr1,vctr2,vctr3,vctr4,vctr5         &
                   ,vctr6,vctr7,jd,scr2,dn0,dtlv,ngrid)

        call truhor(m1,m2,m3,ia,iz,ja,jz                    &
                   ,wp,vt3dk,'ydir','dyv',grid_g(ngrid)%dyv &
                   ,grid_g(ngrid)%topt,grid_g(ngrid)%rtgt   &
                   ,zm,vctr1,vctr2,vctr3,vctr4,vctr5        &
                   ,vctr6,vctr7,jd,scr2,dn0,dtlv,ngrid)

     endif

  case default
     if(ihorgrad.eq.1)then
        call divcart(m1,m2,m3,ia,iz,ja,jz,vt3df,vt3dj,'XDIR','OPNT',ibotflx)
        call divcart(m1,m2,m3,ia,iz,ja,jz,vt3dg,vt3dk,'YDIR','NPNT',ibotflx)

     elseif(ihorgrad.eq.2)then
        call truhor(m1,m2,m3,ia,iz,ja,jz                  &
             ,wp,vt3dj,'xdir','dxu',grid_g(ngrid)%dxu     &
             ,grid_g(ngrid)%topt,grid_g(ngrid)%rtgt       &
             ,zm,vctr1,vctr2,vctr3,vctr4,vctr5            &
             ,vctr6,vctr7,jd,scr2,dn0,dtlv,ngrid)

        call truhor(m1,m2,m3,ia,iz,ja,jz                  &
             ,wp,vt3dk,'ydir','dyv',grid_g(ngrid)%dyv     &
             ,grid_g(ngrid)%topt,grid_g(ngrid)%rtgt       &
             ,zm,vctr1,vctr2,vctr3,vctr4,vctr5            &
             ,vctr6,vctr7,jd,scr2,dn0,dtlv,ngrid)

     endif
  end select

  if (jd .eq. 1) then
     do j = ja,jz
        do i = ia,iz
           do k = 2,m1-2
              wt(k,i,j) = wt(k,i,j)  &
                   + dtlvi * (vt3do(k,i,j) - wp(k,i,j))  &
                   - (vt3dj(k,i,j) + vt3dk(k,i,j))  &
                   / ((dn0(k,i,j) + dn0(k+1,i,j)) * .5)
           enddo
        enddo
     enddo
  else
     do j = ja,jz
        do i = ia,iz
           do k = 2,m1-2
              wt(k,i,j) = wt(k,i,j)  &
                   + dtlvi * (vt3do(k,i,j) - wp(k,i,j))  &
                   - vt3dj(k,i,j)  &
                   / ((dn0(k,i,j) + dn0(k+1,i,j)) * .5)
           enddo
        enddo
     enddo
  endif

  return
end subroutine diffvel

!     ***************************************************************

subroutine diffsclr(m1,m2,m3,ia,iz,ja,jz,jd  &
     ,ia_1,ja_1,ia1,ja1,iz_1,jz_1,iz1,jz1,n,ksf  &
     ,scp,sct,vt3da,vt3db,vt3df,vt3dg  &
     ,vt3dj,vt3dk,vt3do,vt3dc,dn03i,vt3dl,vt3dm,vt2db,rtgt,sfcflx  &
     ,dn0,vkkh,hkkh) !, &
     ! CATT - srf - Large Scale Forcing for GRELL CUPAR - Not used
     ! LFR      iscalar,lsfcupar)
  !!lsfcupar,nsc)

  use mem_grid, only: dtlt      &  !INTENT(IN)
                    , dzm       &  !INTENT(IN)
                    , dzt       &  !INTENT(IN)
                    , nstbot    &  !INTENT(IN)
                    , nsttop    &  !INTENT(IN)
                    , grid_g    &  !INTENT(IN)
                    , ngrid     &  !INTENT(IN)
                    , zt           !INTENT(IN)


  use mem_scratch, only: vctr1   & !INTENT(OUT)
                       , vctr2   & !INTENT(OUT)
                       , vctr3   & !INTENT(OUT)
                       , vctr4   & !INTENT(OUT)
                       , vctr5   & !INTENT(OUT)
                       , vctr6   & !INTENT(OUT)
                       , vctr7     !INTENT(OUT)

  use mem_turb, only: ihorgrad   & !INTENT(IN)
                    , ibotflx    ! !INTENT(IN)

  use var_tables

  ! CATT
  use catt_start, only: CATT           ! intent(in)
  use mem_cuparm, only : nnqparm       ! intent(in)

  implicit none

  integer, intent(IN) :: n,ksf !,nsc

  integer, INTENT(IN) :: m1    &
                       , m2    &
                       , m3    &
                       , ia    &
                       , iz    &
                       , ja    &
                       , jz    &
                       , jd    &
                       , ia_1  &
                       , ja_1  &
                       , ia1   &
                       , ja1   &
                       , iz_1  &
                       , jz_1  &
                       , iz1   &
                       , jz1

  real, dimension(m1,m2,m3), INTENT(IN) :: scp     &
                                         , vkkh    &
                                         , hkkh    &
                                         , dn0

  real, dimension(m1,m2,m3), INTENT(INOUT) :: sct    &
                                            , vt3da  &
                                            , vt3db  &
                                            , vt3dc  &
                                            , vt3df  &
                                            , vt3dg  &
                                            , vt3dj  &
                                            , vt3dk  &
                                            , vt3do  &
                                            , dn03i  &
                                            , vt3dl  &
                                            , vt3dm

  real, dimension(m2,m3), INTENT(IN)    :: sfcflx,rtgt

  real, dimension(m2,m3), INTENT(INOUT) :: vt2db

  !<CATT-BEGIN>
  ! CATT
  !srf- Large Scale Forcing for GRELL CUPAR
  !integer ::  iscalar
  !!real, pointer :: lsfcupar(:,:,:) ! dimension(m1,m2,m3)
  !srf- ------------------------------------
  !<CATT-END>


  !local variables:
  integer :: i,j,k
  integer, save :: ksf_save = 0

  ! **(JP)** vetoriza calculo de c1
  real :: c1(m2,m3)
  ! **(JP)** fim de modificacao
  real :: dtlti
  
  real, parameter :: verysmall=1.e-20

  !          compute vertical diffusion matrix coefficients for scalars

  if (n == 1 .or. ksf /= ksf_save) then
     ksf_save = ksf
     ! **(JP)** vetoriza calculo de c1
     do j = ja,jz
        do i = ia,iz
           c1(i,j) = .5 * dtlt / (rtgt(i,j) * rtgt(i,j))
        end do
     end do
     ! **(JP)** fim de modificacao
     ! **(JP)** quebra laco para vetorizar segundo laco em ij,
     ! **(JP)** fatorando ifs para fora do laco
     do j = ja,jz
        do i = ia,iz
           do k = 1,m1-1
              vctr1(k) = dzm(k) * (vkkh(k,i,j) + vkkh(k+1,i,j))
           enddo
           ! **(JP)** vetoriza calculo de c1
           !c1 = .5 * dtlt / (rtgt(i,j) * rtgt(i,j))
           ! **(JP)** fim de modificacao
           do k = 2,m1-1
              ! **(JP)** vetoriza calculo de c1
              !vt3dj(k,i,j) = -c1 * dzt(k) * vctr1(k-1)
              !vt3dk(k,i,j) = -c1 * dzt(k) * vctr1(k)
              vt3dj(k,i,j) = -c1(i,j) * dzt(k) * vctr1(k-1)
              vt3dk(k,i,j) = -c1(i,j) * dzt(k) * vctr1(k)
              ! **(JP)** fim de modificacao
              vt3do(k,i,j) = dn0(k,i,j) - vt3dj(k,i,j) - vt3dk(k,i,j)
              dn03i(k,i,j) = 1. / dn0(k,i,j)
           enddo
        enddo
     enddo

     if ((nstbot .eq. 1) .and. (nsttop .eq. 1)) then
        do j = ja,jz
           do i = ia,iz
           vt3dj(2,i,j) = 0.
                vt3do(2,i,j) = dn0(2,i,j) - vt3dk(2,i,j)

           vt3dk(m1-1,i,j) = 0.
                vt3do(m1-1,i,j) = dn0(m1-1,i,j) - vt3dj(m1-1,i,j)

           ! Bob (10/27/00):  Evaluate ci2i, cjp1, and cji terms here
           ! (vt2db, vt3dl, vt3dm) for new tridiff2.
           ! With this, no longer need to pass ci or cip1 (vt3do, vt3dk) to tridiff2.

           vt2db(i,j) = 1. / vt3do(2,i,j)
           vt3dl(2,i,j) = vt3dk(2,i,j) * vt2db(i,j)
           enddo
        enddo

     else if (nstbot .eq. 1) then
        do j = ja,jz
           do i = ia,iz
              vt3dj(2,i,j) = 0.
              vt3do(2,i,j) = dn0(2,i,j) - vt3dk(2,i,j)
              vt3dk(m1-1,i,j) = 0.

              ! Bob (10/27/00):  Evaluate ci2i, cjp1, and cji terms here
              ! (vt2db, vt3dl, vt3dm) for new tridiff2.
              ! With this, no longer need to pass ci or cip1 (vt3do, vt3dk) to tridiff2.

              vt2db(i,j) = 1. / vt3do(2,i,j)
              vt3dl(2,i,j) = vt3dk(2,i,j) * vt2db(i,j)
           enddo
        enddo

     else if (nsttop .eq. 1) then
        do j = ja,jz
           do i = ia,iz
              vt3dj(2,i,j) = 0.
              vt3dk(m1-1,i,j) = 0.
              vt3do(m1-1,i,j) = dn0(m1-1,i,j) - vt3dj(m1-1,i,j)

              ! Bob (10/27/00):  Evaluate ci2i, cjp1, and cji terms here
              ! (vt2db, vt3dl, vt3dm) for new tridiff2.
              ! With this, no longer need to pass ci or cip1 (vt3do, vt3dk) to tridiff2.

              vt2db(i,j) = 1. / vt3do(2,i,j)
              vt3dl(2,i,j) = vt3dk(2,i,j) * vt2db(i,j)
           enddo
        enddo

     else
        do j = ja,jz
           do i = ia,iz
              vt3dj(2,i,j) = 0.
              vt3dk(m1-1,i,j) = 0.

              ! Bob (10/27/00):  Evaluate ci2i, cjp1, and cji terms here
              ! (vt2db, vt3dl, vt3dm) for new tridiff2.
              ! With this, no longer need to pass ci or cip1 (vt3do, vt3dk) to tridiff2.

              vt2db(i,j) = 1. / vt3do(2,i,j)
              vt3dl(2,i,j) = vt3dk(2,i,j) * vt2db(i,j)
           enddo
        enddo
     end if

     do j = ja,jz
        do i = ia,iz
           do k = 3,m1-1
              vt3dm(k,i,j) = 1. / (vt3do(k,i,j) - vt3dj(k,i,j) * vt3dl(k-1,i,j))
              vt3dl(k,i,j) = vt3dk(k,i,j) * vt3dm(k,i,j)

           enddo
        enddo
     enddo

     ! **(JP)** fim de modificacao

  endif

  !     compute 2 horizontal scalar gradients needed for dscp/dt

  if (ihorgrad .eq. 1) then
     call grad(m1,m2,m3,1,iz,ja,jz,scp,vt3df,'XDIR','TPNT',ibotflx)
     call grad(m1,m2,m3,ia,iz,1,jz,scp,vt3dg,'YDIR','TPNT',ibotflx)

     do j = ja,jz
        do i = 1,iz
           do k = 1,m1-1
              vt3df(k,i,j) = -vt3df(k,i,j)  &
                   * .5 * (hkkh(k,i,j) + hkkh(k,i+1,j))
           enddo
        enddo
     enddo

     do j = 1,jz
        do i = ia,iz
           do k = 1,m1-1
              vt3dg(k,i,j) = -vt3dg(k,i,j)  &
                   * .5 * (hkkh(k,i,j) + hkkh(k,i,j+jd))
           enddo
        enddo
     enddo
  endif

  !         horizontal flux divergence for scalars

  if (ihorgrad .eq. 1) then
     call divcart(m1,m2,m3,ia,iz,ja,jz,vt3df,vt3da,'XDIR','UPNT',ibotflx)
     call divcart(m1,m2,m3,ia,iz,ja,jz,vt3dg,vt3db,'YDIR','VPNT',ibotflx)

  elseif (ihorgrad .eq. 2) then
     call truhor(m1,m2,m3,ia,iz,ja,jz  &
                ,scp,vt3da,'xdir','dxu',grid_g(ngrid)%dxu  &
                ,grid_g(ngrid)%topt,grid_g(ngrid)%rtgt     &
                ,zt,vctr1,vctr2,vctr3,vctr4,vctr5  &
                ,vctr6,vctr7,jd,hkkh,dn0,dtlt,ngrid)

     call truhor(m1,m2,m3,ia,iz,ja,jz                      &
                ,scp,vt3db,'ydir','dyv',grid_g(ngrid)%dyv  &
                ,grid_g(ngrid)%topt,grid_g(ngrid)%rtgt     &
                ,zt,vctr1,vctr2,vctr3,vctr4,vctr5          &
                ,vctr6,vctr7,jd,hkkh,dn0,dtlt,ngrid)
  endif

  !         finish matrix coefficients

  do j = ja,jz
     do i = ia,iz
        do k = 2,m1-1
           vt3dc(k,i,j) = scp(k,i,j) * dn0(k,i,j)
        enddo

        if (nstbot .eq. 1) then
           vt3dc(2,i,j) = scp(2,i,j) * dn0(2,i,j)  &
                + sfcflx(i,j) * dtlt * dzt(2) / rtgt(i,j)
        else
           vt3dc(2,i,j) = scp(2,i,j) * dn0(2,i,j)  &
                + .5 * dtlt * (vkkh(1,i,j) + vkkh(2,i,j))  &
                * scp(1,i,j) * dzm(2) * dzt(2) / rtgt(i,j) ** 2
        endif

        if (nsttop .eq. 0) then
           vt3dc(m1-1,i,j) = scp(m1-1,i,j) * dn0(m1-1,i,j)  &
                - .5 * dtlt * (vkkh(m1-1,i,j) + vkkh(m1,i,j))  &
                * scp(m1,i,j) * dzm(m1-1) * dzt(m1) / rtgt(i,j) ** 2
        endif
     enddo
  enddo

  !! call tridiff2(m1,m2,m3,ia,iz,ja,jz,m1-1,vt2db,vt3dj,vt3dc,vt3dl,vt3dm,vt3df)

  do j = ja,jz
     do i = ia,iz

        vt3df(2,i,j) = vt3dc(2,i,j) * vt2db(i,j)

        do k = 3,m1-1
           vt3df(k,i,j) = (vt3dc(k,i,j) - vt3dj(k,i,j) * vt3df(k-1,i,j))  &
                * vt3dm(k,i,j)
           if (abs(vt3df(k,i,j)) < verysmall) vt3df(k,i,j) = 0. 
        enddo

        do k = m1-2,2,-1
           vt3df(k,i,j) = vt3df(k,i,j) - vt3dl(k,i,j) * vt3df(k+1,i,j)
           if (abs(vt3df(k,i,j)) < verysmall) vt3df(k,i,j) = 0. 
        enddo

     enddo
  enddo

  dtlti = 1.0 / dtlt

  if (jd .eq. 1) then
     do j = ja,jz
        do i = ia,iz
           do k = 2,m1-1
              sct(k,i,j) = sct(k,i,j) - scp(k,i,j) * dtlti  &
                   - (vt3da(k,i,j) + vt3db(k,i,j)) * dn03i(k,i,j)
           enddo
        enddo
     enddo
  else
     do j = ja,jz
        do i = ia,iz
           do k = 2,m1-1
              sct(k,i,j) = sct(k,i,j) - scp(k,i,j) * dtlti  &
                   - vt3da(k,i,j) * dn03i(k,i,j)
           enddo
        enddo
     enddo
  endif

  ! CATT
!!$  if (CATT == 1) then
!!$     !srf--------- PBL Forcing for GRELL CUPAR (only vertical term) ------
!!$     !  scalar_tab(n,ngrid)%name == 'THP'
!!$     !  if(iscalar.eq.9 .or. iscalar.eq.10) then
!!$     if(scalar_tab(nsc,ngrid)%name == 'THP' .or. &
!!$          scalar_tab(nsc,ngrid)%name == 'RTP') then
!!$
!!$        do j = ja,jz
!!$           do i = ia,iz
!!$              do k = 2,m1-1
!!$
!!$                 !print
!!$                 !          if(iscalar.eq.9) then
!!$                 !           if(j.eq.44.and.i.eq.60)then
!!$                 !
!!$                 !            if(k.eq.2) print*,'----DIFF----2-',j,i
!!$                 !
!!$                 !            print*,k,lsfcupar(k,i,j),&
!!$                 !            lsfcupar(k,i,j) + dtlti*(vt3df(k,i,j)-scp(k,i,j)) ,sct(k,i,j),dtlti*(vt3df(k,i,j)-scp(k,i,j))
!!$                 !           endif
!!$                 !          endif
!!$                 !print
!!$                 !        forcing with radiation plus only the vertical diffusion
!!$                 lsfcupar(k,i,j) = lsfcupar(k,i,j) + dtlti*(vt3df(k,i,j)-scp(k,i,j))
!!$                 !        forcing with horizontal diffusion included:
!!$                 !	 lsfcupar(k,i,j) = lsfcupar(k,i,j) + dtlti*(vt3df(k,i,j)-scp(k,i,j)) &
!!$                 !           - (vt3da(k,i,j) + vt3db(k,i,j)) / dn0(k,i,j)
!!$                 !---------------------------------------
!!$
!!$              enddo
!!$           enddo
!!$        enddo
!!$
!!$     endif
!!$     !-srf----------------
!!$  endif

  do j = ja,jz
     do i = ia,iz
        do k = 2,m1-1
           sct(k,i,j) = sct(k,i,j) + vt3df(k,i,j) * dtlti
        enddo
     enddo
  enddo

!srf--------- PBL Forcing for GRELL CUPAR (only vertical term) ------
!!$  if (CATT == 1 .and. nnqparm(ngrid) == 2) then
!!$     if(scalar_tab(nsc,ngrid)%name == 'THP' .or. &
!!$	scalar_tab(nsc,ngrid)%name == 'RTP'      ) then
!!$
!!$	do j = ja,jz
!!$	   do i = ia,iz
!!$	      do k = 2,m1-1
!!$
!!$		 !	  forcing with radiation plus only the vertical diffusion
!!$		 lsfcupar(k,i,j) = &
!!$                      lsfcupar(k,i,j) + dtlti*(vt3df(k,i,j)-scp(k,i,j))
!!$
!!$		 !if(i==49 .and. j==38) then
!!$		 !  print*,'turb k:',scalar_tab(nsc,ngrid)%name,k,lsfcupar(k,i,j)*86400.,&
!!$		 !		 dtlti*(vt3df(k,i,j)-scp(k,i,j))*86400.
!!$		 !		!if(k==m1-1) stop 4443
!!$		 !endif
!!$		
!!$		 !	  forcing with horizontal diffusion included:
!!$		 !    lsfcupar(k,i,j) = lsfcupar(k,i,j) + dtlti*(vt3df(k,i,j)-scp(k,i,j)) &
!!$		 !	             - (vt3da(k,i,j) + vt3db(k,i,j)) / dn0(k,i,j)
!!$		 !---------------------------------------
!!$
!!$	      enddo
!!$	   enddo
!!$	enddo
!!$        
!!$     endif
!!$  endif
!-srf----------------

  return
end subroutine diffsclr

!     ******************************************************************

subroutine tridiff1(m1,m2,m3,ia,iz,ja,jz,kz,cim1,ci,cip1,rhs,cj,cjp1)

  implicit none

  integer, INTENT(IN) :: m1    &
                       , m2    &
                       , m3    &
                       , ia    &
                       , iz    &
                       , ja    &
                       , jz    &
                       , kz

  real, dimension(m1,m2,m3), INTENT(IN)    :: ci      &
                                            , cip1    &
                                            , cim1

  real, dimension(m1,m2,m3), INTENT(INOUT) :: cjp1    &
                                            , rhs     &
                                            , cj

  !local variables:
  integer :: i,j,k

  !integer lpx(m2,m3)

  real cji

  do j = ja,jz
     do i = ia,iz

        cjp1(2,i,j) = cip1(2,i,j) / ci(2,i,j)
        rhs(2,i,j) = rhs(2,i,j) / ci(2,i,j)

        do k = 3,kz
           cj(k,i,j) = ci(k,i,j) - cim1(k,i,j) * cjp1(k-1,i,j)
           cji = 1. / cj(k,i,j)
           cjp1(k,i,j) = cip1(k,i,j) * cji
           rhs(k,i,j) = (rhs(k,i,j) - cim1(k,i,j) * rhs(k-1,i,j))  &
                * cji
        enddo

        cj(kz,i,j) = rhs(kz,i,j)

        do k = kz-1,2,-1
           cj(k,i,j) = rhs(k,i,j) - cjp1(k,i,j) * cj(k+1,i,j)
        enddo

     enddo
  enddo
  return
end subroutine tridiff1

!     ******************************************************************

subroutine tridiff2(m1,m2,m3,ia,iz,ja,jz,kz,ci2i,cim1,rhs,cjp1,cji,out)

  implicit none

  integer, INTENT(IN) :: m1   &
                       , m2   &
                       , m3   &
                       , ia   &
                       , iz   &
                       , ja   &
                       , jz   &
                       , kz

  real, dimension(m2,m3), INTENT(IN)       :: ci2i

  real, dimension(m1,m2,m3), INTENT(IN)    :: cim1    &
                                            , rhs     &
                                            , cjp1    &
                                            , cji

  real, dimension(m1,m2,m3), INTENT(INOUT) :: out

  !local variables:
  integer :: i,j,k
  !integer lpx(m2,m3)

  do j = ja,jz
     do i = ia,iz
        out(2,i,j) = rhs(2,i,j) * ci2i(i,j)
     enddo
  enddo

  do k = 3,kz
     do j = ja,jz
        do i = ia,iz
           out(k,i,j) = (rhs(k,i,j) - cim1(k,i,j) * out(k-1,i,j)) * cji(k,i,j)
        enddo
     enddo
  enddo

  do k = kz-1,2,-1
     do j = ja,jz
        do i = ia,iz
           out(k,i,j) = out(k,i,j) - cjp1(k,i,j) * out(k+1,i,j)
        enddo
     enddo
  enddo
  return
end subroutine tridiff2

!     ******************************************************************

subroutine tridiff2orig(m1,m2,m3,ia,iz,ja,jz,kz,ci2i,cim1,rhs,cjp1,cji,out)

  implicit none

  integer, INTENT(IN) :: m1   &
                       , m2   &
                       , m3   &
                       , ia   &
                       , iz   &
                       , ja   &
                       , jz   &
                       , kz

  real, dimension(m2,m3),    INTENT(IN)    :: ci2i

  real, dimension(m1,m2,m3), INTENT(IN)    :: cim1,rhs,cjp1,cji

  real, dimension(m1,m2,m3), INTENT(INOUT) :: out


  !local variables:
  integer :: i,j,k

  !integer lpx(m2,m3)

  do j = ja,jz
     do i = ia,iz

        out(2,i,j) = rhs(2,i,j) * ci2i(i,j)

        do k = 3,kz
           out(k,i,j) = (rhs(k,i,j) - cim1(k,i,j) * out(k-1,i,j)) * cji(k,i,j)
        enddo

        do k = kz-1,2,-1
           out(k,i,j) = out(k,i,j) - cjp1(k,i,j) * out(k+1,i,j)
        enddo

     enddo
  enddo
  return
end subroutine tridiff2orig









!     ******************************************************************

subroutine truhor(m1,m2,m3,ia,iz,ja,jz  &
     ,                 vc3da,vc3db,dir,gpnt,dxy,topo,rtg  &
     ,                 z,vctr1,vctr2,vctr3,zintl,zintr  &
     ,                 dxynu,dxytem,jd,dn0hkh,dn0,dtl,ngrid)

  implicit none

  integer, INTENT(IN) :: m1    &
                       , m2    &
                       , m3    &
                       , ia    &
                       , iz    &
                       , ja    &
                       , jz    &
                       , jd    &
                       , ngrid

  real, dimension(m1,m2,m3), INTENT(IN)    :: vc3da    &
                                            , dn0hkh   &
                                            , dn0

  real, dimension(m1,m2,m3), INTENT(INOUT) :: vc3db

  real, dimension(m2,m3),    INTENT(IN)    :: dxy    &
                                            , topo   &
                                            , rtg

  real, dimension(m1),       INTENT(IN)    :: z

  real, dimension(m1),      INTENT(OUT)    :: vctr1   &
                                            , vctr2   &
                                            , vctr3   &
                                            , dxynu   &
                                            , dxytem  &
                                            , zintl   &
                                            , zintr

  character(len=*),  INTENT(IN) :: dir,gpnt

  !local variables:
  integer :: nz,jaa,jzz,nclip,je,jf,i,j,k

  real :: delz1,delz3,distnu,distold,scale,aknu,akdiff,dtl


  nz=m1

  jaa=ja
  jzz=jz
  if(jd.eq.0) then
     jaa=1
     jzz=1
  endif

  !      if(dtl.eq.999.)then

  nclip=0
  if(dir.eq.'xdir')then
     if(gpnt.eq.'dyt'.or.gpnt.eq.'dxt')then
        je=0
        jf=1
     else
        je=-1
        jf=0
     endif

     do j=jaa,jzz
        do i=ia,iz
           do k=1,m1
              vctr1(k)=topo(i-1,j)+z(k)*rtg(i-1,j)
              vctr2(k)=topo(i  ,j)+z(k)*rtg(i ,j)
              vctr3(k)=topo(i+1,j)+z(k)*rtg(i+1,j)
              dxytem(k)=vc3da(k,i-1,j)
              dxynu(k)=vc3da(k,i+1,j)
           enddo

           !             extrapolate the gradient from 3--->2 to 2--->1.
           delz1=(vctr1(2)-vctr1(1))/(vctr1(3)-vctr1(2))
           delz3=(vctr3(2)-vctr3(1))/(vctr3(3)-vctr3(2))
           dxytem(1)=dxytem(2)-(dxytem(3)-dxytem(2))*delz1
           dxynu(1) =dxynu(2) -(dxynu(3) -dxynu(2)) *delz3

           call htint(nz,dxytem,vctr1,nz,zintl,vctr2)
           call htint(nz,dxynu, vctr3,nz,zintr,vctr2)

           !             check intersection of cartesian sfc. on left side of mtn.
           call topobnd(m1,m2,m3,i,j, 1,0,'l',vc3da,dxy(i+je,j)  &
                ,                    z,vctr1,vctr2,vctr3,zintr,dxytem,jd)
           !             check intersection of cartesian sfc. on right side of mtn.
           call topobnd(m1,m2,m3,i,j,-1,0,'r',vc3da,dxy(i+jf,j)  &
                ,                    z,vctr1,vctr2,vctr3,zintl,dxynu,jd)

           do k=2,m1
              distnu=1./dxytem(k)+1./dxynu(k)
              distold=1./dxy(i+je,j)+1./dxy(i+jf,j)
              scale=distnu**2/distold**2
              !                avgden=0.5*(dn0(k,i,j)+dn0(k,i+je+jf,j))
              !cfast                if(scale.gt.1)stop' scale incorrect-xdir--'
              !c            clp=min(dn0hkh(k,i,j)*scale,0.124/16.*avgden*distnu**2/dtl)
              !                clp=min(dn0hkh(k,i,j),0.124*avgden*distold**2/dtl)
              !cfast                if(dn0hkh(k,i,j).ne.clp)nclip=nclip+1
              !               Need to multiply hkh by 2 since going over 2dx
              aknu=dn0hkh(k,i,j)*scale*2.
              !                aknu=clp*scale
              akdiff=-aknu/distnu

              vc3db(k,i,j)=akdiff*( (zintr(k)-vc3da(k,i,j))*dxytem(k)  &
                   -(vc3da(k,i,j)-zintl(k))*dxynu(k) )

              !                print*,'distnu,old,scale,tend=',k,i,j,distnu
              !     &,                 distold,scale,vc3db(k,i,j)
              !               if(i.eq.ip.and.j.eq.jp.and.k.le.3
              !     &            .and.vc3da(k,i,j).gt.100.)
              !     &        print*,'aknu,scal,lmr=',k,i,j,akdiff,scale,zintl(k)
              !     &,                            vc3da(k,i,j),zintr(k),vc3db(k,i,j)
           enddo
           !              vc3db(2,i,j)=0.
           vc3db(1,i,j)=vc3db(2,i,j)
           if(gpnt.eq.'wpnt')then
              vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
           else
              vc3db(m1,i,j)=vc3db(m1-1,i,j)
           endif
        enddo
     enddo

  elseif(dir.eq.'ydir')then

     if(gpnt.eq.'dyt'.or.gpnt.eq.'dxt')then
        je=0
        jf=jd
     else
        je=-jd
        jf=0
     endif

     do j=jaa,jzz
        do i=ia,iz
           do k=1,m1
              vctr1(k)=topo(i,j-jd)+z(k)*rtg(i,j-jd)
              vctr2(k)=topo(i,j)   +z(k)*rtg(i,j)
              vctr3(k)=topo(i,j+jd)+z(k)*rtg(i,j+jd)
              dxytem(k)=vc3da(k,i,j-jd)
              dxynu(k)=vc3da(k,i,j+jd)
           enddo

           !             extrapolate the gradient from 3--->2 to 2--->1.
           delz1=(vctr1(2)-vctr1(1))/(vctr1(3)-vctr1(2))
           delz3=(vctr3(2)-vctr3(1))/(vctr3(3)-vctr3(2))
           dxytem(1)=dxytem(2)-(dxytem(3)-dxytem(2))*delz1
           dxynu(1) =dxynu(2) -(dxynu(3) -dxynu(2)) *delz3

           call htint(nz,dxytem,vctr1,nz,zintl,vctr2)
           call htint(nz,dxynu, vctr3,nz,zintr,vctr2)

           !             check intersection of cartesian sfc. on left side of mtn.
           call topobnd(m1,m2,m3,i,j,0, 1,'l',vc3da,dxy(i,j+je)  &
                ,                    z,vctr1,vctr2,vctr3,zintr,dxytem,jd)
           !             check intersection of cartesian sfc. on right side of mtn.
           call topobnd(m1,m2,m3,i,j,0,-1,'r',vc3da,dxy(i,j+jf)  &
                ,                    z,vctr1,vctr2,vctr3,zintl,dxynu,jd)
           do k=2,m1
              distnu=1./dxytem(k)+1./dxynu(k)
              distold=1./dxy(i,j+je)+1./dxy(i,j+jf)
              scale=distnu**2/distold**2
              !cfast                if(scale.gt.1)stop' scale incorrect-ydir--'
              !                avgden=0.5*(dn0(k,i,j)+dn0(k,i,j+je+jf))
              !c            clp=min(dn0hkh(k,i,j)*scale,0.124/16.*avgden*distnu**2/dtl)
              !                clp=min(dn0hkh(k,i,j),0.124*avgden*distold**2/dtl)
              !cfast                if(dn0hkh(k,i,j).ne.clp)nclip=nclip+1
              !               Need to multiply hkh by 2 since going over 2dx
              aknu=dn0hkh(k,i,j)*scale*2.
              !                aknu=clp*scale
              akdiff=-aknu/distnu
              vc3db(k,i,j)=akdiff*( (zintr(k)-vc3da(k,i,j))*dxytem(k)  &
                   -(vc3da(k,i,j)-zintl(k))*dxynu(k) )
           enddo
           !              vc3db(2,i,j)=0.
           vc3db(1,i,j)=vc3db(2,i,j)
           if(gpnt.eq.'npnt')then
              vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
           else
              vc3db(m1,i,j)=vc3db(m1-1,i,j)
           endif
        enddo
     enddo

  endif

  if(nclip.ge.1.and.gpnt.eq.'dxt')  &
       print*,'clip on grd',ngrid,' dir=',dir(1:1),' ',gpnt(1:3),nclip

  !      endif

  return
end subroutine truhor




!----------------------TRUHOR_OPT------------------------------------
!     ******************************************************************

subroutine truhor_opt(m1,m2,m3,ia,iz,ja,jz                   &
     ,                vc3da,vc3db,dir,gpnt,dxy,topo,rtg      &
     ,                z,vctr1,vctr2,vctr3,zintl,zintr        &
     ,                dxynu,dxytem,jd,dn0hkh,dn0,dtl,ngrid)

  use mem_opt, only : opt

  implicit none

  integer, INTENT(IN) :: m1    &
                       , m2    &
                       , m3    &
                       , ia    &
                       , iz    &
                       , ja    &
                       , jz    &
                       , jd    &
                       , ngrid

  real, dimension(m1,m2,m3),    INTENT(IN) :: vc3da     &
                                            , dn0hkh    &
                                            , dn0

  real, dimension(m1,m2,m3), INTENT(INOUT) :: vc3db

  real, dimension(m2,m3), INTENT(IN) :: dxy,topo,rtg

  real, dimension(m1), INTENT(IN)    :: z

  real, dimension(m1), INTENT(OUT)   :: vctr1   &
                                      , vctr2   &
                                      , vctr3   &
                                      , zintl   &
                                      , zintr   &
                                      , dxynu   &
                                      , dxytem

  character(len=*), INTENT(IN) ::  dir,gpnt


  integer :: nz,jaa,jzz,nclip,je,jf,i,j,k

  real :: delz1,delz3,distnu,distold,scale,aknu,akdiff,dtl

  integer      :: htint_i, htint_j

  nz=m1

  jaa=ja
  jzz=jz
  if(jd.eq.0) then
     jaa=1
     jzz=1
  endif

  !      if(dtl.eq.999.)then

  nclip=0
  if(dir.eq.'xdir')then
     if(gpnt.eq.'dyt'.or.gpnt.eq.'dxt')then
        je=0
        jf=1
     else
        je=-1
        jf=0
     endif

     do j=jaa,jzz
        do i=ia,iz
           do k=1,m1
              vctr1(k)=topo(i-1,j)+z(k)*rtg(i-1,j)
              vctr2(k)=topo(i  ,j)+z(k)*rtg(i ,j)
              vctr3(k)=topo(i+1,j)+z(k)*rtg(i+1,j)
              dxytem(k)=vc3da(k,i-1,j)
              dxynu(k)=vc3da(k,i+1,j)
           enddo

           !             extrapolate the gradient from 3--->2 to 2--->1.
           delz1=(vctr1(2)-vctr1(1))/(vctr1(3)-vctr1(2))
           delz3=(vctr3(2)-vctr3(1))/(vctr3(3)-vctr3(2))
           dxytem(1)=dxytem(2)-(dxytem(3)-dxytem(2))*delz1
           dxynu(1) =dxynu(2) -(dxynu(3) -dxynu(2)) *delz3

           htint_i = i-ia +1
           htint_j = j-jaa+1

           call htint_inter(nz, dxytem, nz, zintl,  &
                opt%ind1_x_a(:,htint_i,htint_j),    &
                opt%ind2_x_a(:,htint_i,htint_j),    &
                opt%weight_x_a(:,htint_i,htint_j))

           call htint_inter(nz, dxynu,  nz, zintr,  &
                opt%ind1_x_b(:,htint_i,htint_j),    &
                opt%ind2_x_b(:,htint_i,htint_j),    &
                opt%weight_x_b(:,htint_i,htint_j))

           !  check intersection of cartesian sfc. on left side of mtn.
           call topobnd(m1,m2,m3,i,j, 1,0,'l',vc3da,dxy(i+je,j)  &
                ,                    z,vctr1,vctr2,vctr3,zintr,dxytem,jd)
           !  check intersection of cartesian sfc. on right side of mtn.
           call topobnd(m1,m2,m3,i,j,-1,0,'r',vc3da,dxy(i+jf,j)  &
                ,                    z,vctr1,vctr2,vctr3,zintl,dxynu,jd)

           do k=2,m1
              distnu=1./dxytem(k)+1./dxynu(k)
              distold=1./dxy(i+je,j)+1./dxy(i+jf,j)
              scale=distnu**2/distold**2
              !      avgden=0.5*(dn0(k,i,j)+dn0(k,i+je+jf,j))
              !cfast if(scale.gt.1)stop' scale incorrect-xdir--'
              !c    clp=min(dn0hkh(k,i,j)*scale,0.124/16.*avgden*distnu**2/dtl)
              !     clp=min(dn0hkh(k,i,j),0.124*avgden*distold**2/dtl)
              !cfast if(dn0hkh(k,i,j).ne.clp)nclip=nclip+1
              !      Need to multiply hkh by 2 since going over 2dx
              aknu=dn0hkh(k,i,j)*scale*2.
              !                aknu=clp*scale
              akdiff=-aknu/distnu

              vc3db(k,i,j)=akdiff*( (zintr(k)-vc3da(k,i,j))*dxytem(k)  &
                   -(vc3da(k,i,j)-zintl(k))*dxynu(k) )

              !               print*,'distnu,old,scale,tend=',k,i,j,distnu
              !     &,                distold,scale,vc3db(k,i,j)
              !              if(i.eq.ip.and.j.eq.jp.and.k.le.3
              !     &           .and.vc3da(k,i,j).gt.100.)
              !     &       print*,'aknu,scal,lmr=',k,i,j,akdiff,scale,zintl(k)
              !     &,                      vc3da(k,i,j),zintr(k),vc3db(k,i,j)
           enddo
           !              vc3db(2,i,j)=0.
           vc3db(1,i,j)=vc3db(2,i,j)
           if(gpnt.eq.'wpnt')then
              vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
           else
              vc3db(m1,i,j)=vc3db(m1-1,i,j)
           endif
        enddo
     enddo

  elseif(dir.eq.'ydir')then

     if(gpnt.eq.'dyt'.or.gpnt.eq.'dxt')then
        je=0
        jf=jd
     else
        je=-jd
        jf=0
     endif

     do j=jaa,jzz
        do i=ia,iz
           do k=1,m1
              vctr1(k)=topo(i,j-jd)+z(k)*rtg(i,j-jd)
              vctr2(k)=topo(i,j)   +z(k)*rtg(i,j)
              vctr3(k)=topo(i,j+jd)+z(k)*rtg(i,j+jd)
              dxytem(k)=vc3da(k,i,j-jd)
              dxynu(k)=vc3da(k,i,j+jd)
           enddo

           !             extrapolate the gradient from 3--->2 to 2--->1.
           delz1=(vctr1(2)-vctr1(1))/(vctr1(3)-vctr1(2))
           delz3=(vctr3(2)-vctr3(1))/(vctr3(3)-vctr3(2))
           dxytem(1)=dxytem(2)-(dxytem(3)-dxytem(2))*delz1
           dxynu(1) =dxynu(2) -(dxynu(3) -dxynu(2)) *delz3

           htint_i = i-ia +1
           htint_j = j-jaa+1

           call htint_inter(nz, dxytem, nz, zintl,  &
                opt%ind1_y_a(:,htint_i,htint_j),    &
                opt%ind2_y_a(:,htint_i,htint_j),    &
                opt%weight_y_a(:,htint_i,htint_j))
           call htint_inter(nz, dxynu,  nz, zintr,  &
                opt%ind1_y_b(:,htint_i,htint_j),    &
                opt%ind2_y_b(:,htint_i,htint_j),    &
                opt%weight_y_b(:,htint_i,htint_j))

           !         check intersection of cartesian sfc. on left side of mtn.
           call topobnd(m1,m2,m3,i,j,0, 1,'l',vc3da,dxy(i,j+je)  &
                ,                    z,vctr1,vctr2,vctr3,zintr,dxytem,jd)
           !         check intersection of cartesian sfc. on right side of mtn.
           call topobnd(m1,m2,m3,i,j,0,-1,'r',vc3da,dxy(i,j+jf)  &
                ,                    z,vctr1,vctr2,vctr3,zintl,dxynu,jd)
           do k=2,m1
              distnu=1./dxytem(k)+1./dxynu(k)
              distold=1./dxy(i,j+je)+1./dxy(i,j+jf)
              scale=distnu**2/distold**2
              !cfast    if(scale.gt.1)stop' scale incorrect-ydir--'
              !     avgden=0.5*(dn0(k,i,j)+dn0(k,i,j+je+jf))
              !c    clp=min(dn0hkh(k,i,j)*scale,0.124/16.*avgden*distnu**2/dtl)
              !     clp=min(dn0hkh(k,i,j),0.124*avgden*distold**2/dtl)
              !cfast   if(dn0hkh(k,i,j).ne.clp)nclip=nclip+1
              !        Need to multiply hkh by 2 since going over 2dx
              aknu=dn0hkh(k,i,j)*scale*2.
              !                aknu=clp*scale
              akdiff=-aknu/distnu
              vc3db(k,i,j)=akdiff*( (zintr(k)-vc3da(k,i,j))*dxytem(k)  &
                   -(vc3da(k,i,j)-zintl(k))*dxynu(k) )
           enddo
           !              vc3db(2,i,j)=0.
           vc3db(1,i,j)=vc3db(2,i,j)
           if(gpnt.eq.'npnt')then
              vc3db(m1-1,i,j)=vc3db(m1-2,i,j)
           else
              vc3db(m1,i,j)=vc3db(m1-1,i,j)
           endif
        enddo
     enddo

  endif

  if(nclip.ge.1.and.gpnt.eq.'dxt')  &
       print*,'clip on grd',ngrid,' dir=',dir(1:1),' ',gpnt(1:3),nclip

  !      endif

  return
end subroutine truhor_opt

!---------------------------------------------------------------------
subroutine topobnd(m1,m2,m3,i,j,ipm1,jpm1,sid,vc3da,dxy  &
     ,                  z,vctr1,vctr2,vctr3,zint,dxynu,jd)

  implicit none

  integer, INTENT(IN)   :: m1    &
                         , m2    &
                         , m3    &
                         , ipm1  &
                         , jpm1  &
                         , jd

  real, INTENT(IN) :: dxy

  real, dimension(m1,m2,m3), INTENT(IN) :: vc3da

  real, dimension(m1), INTENT(IN)       :: vctr1   &
                                         , vctr2   &
                                         , z       &
                                         , vctr3

  real, dimension(m1), INTENT(OUT)      :: dxynu

  real, dimension(m1), INTENT(INOUT)    :: zint

  character(len=*), INTENT(IN)          :: sid

  integer :: i,j,ip,jp,ij,k

  real :: tanth,pct,delz1,delz2,delz3,vc1,vc3,vc2



  ip=6
  jp=2

  ij=ipm1+jpm1

  !     keep track of the distance between scalar points.
  call ae0(m1,dxynu,dxy)

  !                        |
  !                    ----|----
  !                   /    |    \
  !                  /     |     \
  !                 /      |      \
  !                /       |       \
  !               /        |        \
  !              /         |         \
  !             /          |          \
  !            /           |           \
  !      ------            |            -------
  !     |      |      |    |   |     |     |
  !     vctr1  |   vctr3   | vctr1   |     vctr3
  !           vctr2        |       vctr2
  !     t(i-1) t(i) t(i+1) | t(i-1) t(i)   t(i+1)
  !
  !     if the cartesian slice intersects topo.,
  !       find the slope of the topography (tanth), then calculate the
  !       horizontal cartesian distance to sigma=2 level (1/dxynu).
  !       then interpolate along the sigma surface to find the value of
  !       the variable at this intersection.

  !     check left side of hill for rhs intersection.
  if(sid.eq.'l')then

     do k=3,m1
        if(vctr3(2).le.vctr2(k))goto 10
        tanth=(vctr3(2)-vctr2(2))*dxy
        dxynu(k)=tanth/(vctr2(k)-vctr2(2))
        pct=dxy/dxynu(k)
        zint(k)=vc3da(2,i+ipm1,j+jpm1)*pct+vc3da(2,i,j)*(1.-pct)
     enddo
10   continue
     if(vctr3(1).le.vctr2(2))return

     tanth=(vctr3(1)-vctr2(1))*dxy
     dxynu(2)=tanth/(vctr2(2)-vctr2(1))
     pct=dxy/dxynu(2)

     !       extrapolate the gradient from 3--->2 to 2--->1.
     delz1=(vctr1(2)-vctr1(1))/(vctr1(3)-vctr1(2))
     delz2=(vctr2(2)-vctr2(1))/(vctr2(3)-vctr2(2))
     delz3=(vctr3(2)-vctr3(1))/(vctr3(3)-vctr3(2))
     vc1=vc3da(2,i-ipm1,j-jpm1)-(vc3da(3,i-ipm1,j-jpm1)  &
          -vc3da(2,i-ipm1,j-jpm1))*delz1
     vc3=vc3da(2,i+ipm1,j+jpm1)-(vc3da(3,i+ipm1,j+jpm1)  &
          -vc3da(2,i+ipm1,j+jpm1))*delz3
     !-----------------------------------------------------------------------
     !       use the first vc2 formulation for strictly hh runs, as this is
     !         the 'correct' value to use as an endpoint when interpolating
     !         along a sigma surface. the second formulation for vc2 is a
     !         'fix' around the excessive gridpoint cooling caused by
     !         extrapolating a locally cold gridpoint to k=1, and then using
     !         this excessively cold value as an endpoint. this then forces
     !         the gridpoint cooler at k=2 since the interpolation along the
     !         k=1 sigma surface will results in a value colder than the k=2
     !         gridpoint value. while a hh run with radiate and sfclyr
     !         commented out, lsflg=3, nudlat=0, initial=1 and us,vs=0 will
     !         not stay hh, it appears that simulations with the standard
     !         rturb.f and this version produce very similar results. this
     !         gives confidence to the solution when using the substantially
     !         higher resolution allowed by this formulation.
     !       note that this scheme does not conserve energy in its present
     !         formulation since it diffuses according to a simple 1:2:1
     !         scheme.
     !       it also appears that while vertical resolution can be as high
     !         as 20m with no winds, experiments with the real time
     !         forecasting system indicate a resolution nearer 70m is needed
     !         since problems with horizontal gradients in a sigma coordinate
     !         system near steep topography exist in other routines as well.
     !-----------------------------------------------------------------------
     !        vc2=vc3da(2,i,j)-(vc3da(3,i,j)-vc3da(2,i,j))*delz2
     vc2=0.5*(vc3+vc1)
     zint(2)=vc3*pct+vc2*(1.-pct)

     !     check right side of hill for lhs intersection.
  elseif(sid.eq.'r')then

     do k=3,m1
        if(vctr1(2).le.vctr2(k))goto 20
        tanth=(vctr1(2)-vctr2(2))*dxy
        dxynu(k)=tanth/(vctr2(k)-vctr2(2))
        pct=dxy/dxynu(k)
        zint(k)=vc3da(2,i+ipm1,j+jpm1)*pct+vc3da(2,i,j)*(1.-pct)
     enddo
20   continue
     if(vctr1(1).le.vctr2(2))return

     tanth=(vctr1(1)-vctr2(1))*dxy
     dxynu(2)=tanth/(vctr2(2)-vctr2(1))
     pct=dxy/dxynu(2)

     !       extrapolate the gradient from 3--->2 to 2--->1.
     delz1=(vctr1(2)-vctr1(1))/(vctr1(3)-vctr1(2))
     delz2=(vctr2(2)-vctr2(1))/(vctr2(3)-vctr2(2))
     delz3=(vctr3(2)-vctr3(1))/(vctr3(3)-vctr3(2))
     vc1=vc3da(2,i+ipm1,j+jpm1)-(vc3da(3,i+ipm1,j+jpm1)  &
          -vc3da(2,i+ipm1,j+jpm1))*delz1
     vc3=vc3da(2,i-ipm1,j-jpm1)-(vc3da(3,i-ipm1,j-jpm1)  &
          -vc3da(2,i-ipm1,j-jpm1))*delz3
     !        vc2=vc3da(2,i,j)-(vc3da(3,i,j)-vc3da(2,i,j))*delz2
     vc2=0.5*(vc3+vc1)
     zint(2)=vc1*pct+vc2*(1.-pct)

  endif

  return
end subroutine topobnd

!-----------------------------------------------------------------------

subroutine avgvel(m1,m2,m3,ia,iz,ja,jz,dir,jd,a,b)

  implicit none


  integer, INTENT(IN) :: m1   &
                       , m2   &
                       , m3   &
                       , ia   &
                       , iz   &
                       , ja   &
                       , jz   &
                       , jd

  real, dimension(m1,m2,m3), INTENT(IN)     :: b

  real, dimension(m1,m2,m3), INTENT(INOUT)  :: a

  character*(*), INTENT(IN) :: dir

  !local variables
  integer :: jaa,jzz,i,j,k

  jaa=ja
  jzz=jz
  if(jd.eq.0) then
     jaa=1
     jzz=1
  endif

  if(dir.eq.'xdir')then
     do j=jaa,jzz
        do i=ia,iz
           do k=1,m1
              a(k,i,j)=0.5*(b(k,i,j)+b(k,i+1,j))
           enddo
        enddo
     enddo
  elseif(dir.eq.'ydir')then
     do j=jaa,jzz
        do i=ia,iz
           do k=1,m1
              a(k,i,j)=0.5*(b(k,i,j)+b(k,i,j+jd))
           enddo
        enddo
     enddo
  endif

  return
end subroutine avgvel
