subroutine diffsclr_brams31(m1,m2,m3,ia,iz,ja,jz,jd,  &
     ia_1,ja_1,ia1,ja1,iz_1,jz_1,iz1,jz1,n,ksf,  &
     scp,sct,vt3da,vt3db,vt3df,vt3dg,  &
     vt3dj,vt3dk,vt3do,vt3dc,dn03i,vt3dl,vt3dm,vt2db,rtgt,sfcflx,  &
     dn0,vkkh,hkkh)

  use mem_grid, only:    &
       dzm,     &      !INTENT(IN)
       dtlt,    &      !INTENT(IN)
       dzt,     &      !INTENT(IN)
       nstbot,  &      !INTENT(IN)
       nsttop,  &      !INTENT(IN)
       grid_g,  &      ! %dxu(IN), %dyv(IN), %rtgt(IN), %topt(IN)
       ngrid,   &      !INTENT(IN)
       zt              !INTENT(IN)

  use mem_scratch, only :   &
       vctr1,   &   !INTENT(INOUT)
       vctr2,   &   !INTENT(INOUT)
       vctr3,   &   !INTENT(INOUT)
       vctr4,   &   !INTENT(INOUT)
       vctr5,   &   !INTENT(INOUT)
       vctr6,   &   !INTENT(INOUT)
       vctr7        !INTENT(INOUT)

  implicit none

  integer, INTENT(IN) :: m1     &
                       , m2     &
                       , m3     &                     
                       , ia     &
                       , iz     &                       
                       , ja     &
                       , jz     &
                       , jd     &
                       , ia_1   &                 ! not referenced
                       , ja_1   &                 ! not referenced
                       , ia1    &                 ! not referenced
                       , ja1    &                 ! not referenced
                       , iz_1   &                 ! not referenced
                       , jz_1   &                 ! not referenced
                       , iz1    &                 ! not referenced
                       , jz1    &                 ! not referenced
                       , n      &
                       , ksf

  real, INTENT(IN) :: scp(m1,m2,m3)    &                                                                                                
                    , vt3dg(m1,m2,m3)  &                                                                                                  
                    , rtgt(m2,m3)      &
                    , sfcflx(m2,m3)    &
                    , dn0(m1,m2,m3)    &
                    , vkkh(m1,m2,m3)   &
                    , hkkh(m1,m2,m3)


  real, INTENT(INOUT) :: vt3da(m1,m2,m3)  &  
                       , vt3dj(m1,m2,m3)  &
                       , vt3dk(m1,m2,m3)  &
                       , vt3do(m1,m2,m3)  &
                       , dn03i(m1,m2,m3)  &
                       , vt3df(m1,m2,m3)  &
                       , sct(m1,m2,m3)    &
                       , vt3db(m1,m2,m3)  &
                       , vt3dm(m1,m2,m3)  &
                       , vt3dl(m1,m2,m3)  &
                       , vt3dc(m1,m2,m3)  &
                       , vt2db(m2,m3)                            

  integer :: i
  integer :: j
  integer :: k
  real    :: c1
  real    :: dtlti
  integer, save :: ksf_save = 0

  !! For optimization
  integer      :: htint_i, htint_j


  ! compute vertical diffusion matrix coefficients for scalars

  if (n == 1 .or. ksf /= ksf_save) then
     ksf_save = ksf
     do j = ja,jz
        do i = ia,iz
           do k = 1,m1-1
              vctr1(k) = dzm(k) * (vkkh(k,i,j) + vkkh(k+1,i,j))
           enddo
           c1 = .5 * dtlt / (rtgt(i,j) * rtgt(i,j))
           do k = 2,m1-1
              vt3dj(k,i,j) = -c1 * dzt(k) * vctr1(k-1)
              vt3dk(k,i,j) = -c1 * dzt(k) * vctr1(k)
              vt3do(k,i,j) = dn0(k,i,j) - vt3dj(k,i,j) - vt3dk(k,i,j)
              dn03i(k,i,j) = 1. / dn0(k,i,j)
           enddo
           vt3dj(2,i,j) = 0.
           if (nstbot .eq. 1)  &
                vt3do(2,i,j) = dn0(2,i,j) - vt3dk(2,i,j)
           vt3dk(m1-1,i,j) = 0.
           if (nsttop .eq. 1)  &
                vt3do(m1-1,i,j) = dn0(m1-1,i,j) - vt3dj(m1-1,i,j)

           ! Bob (10/27/00):  Evaluate ci2i, cjp1, and cji terms here 
           ! (vt2db, vt3dl, vt3dm) for new tridiff2.  
           ! With this, no longer need to pass ci or cip1 (vt3do, vt3dk) 
           ! to tridiff2.

           vt2db(i,j) = 1. / vt3do(2,i,j)
           vt3dl(2,i,j) = vt3dk(2,i,j) * vt2db(i,j)
           do k = 3,m1-1
              vt3dm(k,i,j) = 1. / (vt3do(k,i,j) - vt3dj(k,i,j)*vt3dl(k-1,i,j))
              vt3dl(k,i,j) = vt3dk(k,i,j)*vt3dm(k,i,j)
           enddo
        enddo
     enddo

  endif

     !truhor_opt
     call truhor_opt(m1,m2,m3,ia,iz,ja,jz                           & ! truhor_opt
          ,          scp,vt3da,'xdir','dxu',grid_g(ngrid)%dxu(1,1)  &
          ,          grid_g(ngrid)%topt(1,1)                        &
          ,          grid_g(ngrid)%rtgt(1,1)                        &
          ,          zt,vctr1,vctr2,vctr3,vctr4,vctr5               &
          ,          vctr6,vctr7,jd,hkkh,dn0,dtlt,ngrid)

     !truhor_opt
     call truhor_opt(m1,m2,m3,ia,iz,ja,jz                           & ! truhor_opt
          ,          scp,vt3db,'ydir','dyv',grid_g(ngrid)%dyv(1,1)  &
          ,          grid_g(ngrid)%topt(1,1)                        &
          ,          grid_g(ngrid)%rtgt(1,1)                        &
          ,          zt,vctr1,vctr2,vctr3,vctr4,vctr5               &
          ,          vctr6,vctr7,jd,hkkh,dn0,dtlt,ngrid)

  ! finish matrix coefficients

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

  do j = ja,jz
     do i = ia,iz

        vt3df(2,i,j) = vt3dc(2,i,j) * vt2db(i,j)

        do k = 3,m1-1
           vt3df(k,i,j) = (vt3dc(k,i,j) - vt3dj(k,i,j) * vt3df(k-1,i,j))  &
                * vt3dm(k,i,j)
        enddo

        do k = m1-2,2,-1
           vt3df(k,i,j) = vt3df(k,i,j) - vt3dl(k,i,j) * vt3df(k+1,i,j)
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

  do j = ja,jz
     do i = ia,iz
        do k = 2,m1-1
           sct(k,i,j) = sct(k,i,j) + vt3df(k,i,j) * dtlti
        enddo
     enddo
  enddo
end subroutine diffsclr_brams31
