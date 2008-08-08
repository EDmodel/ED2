!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine strain_adap(m1,m2,m3,ia,iz,ja,jz  &
   ,ia_1,ja_1,iz1,jz1,jd,flpu,flpv,flpw  &
   ,up,vp,wp,vt3da,vt3db,vt3dc,vt3dd,vt3de  &
   ,vt3df,vt3dg,vt3dh,vt3di,vt3dn,scr2,idiffk  &
   ,dxm,dxt,dxu,dxv,dym,dyt,dyu,dyv,dzm,dzt)

implicit none
integer :: m1,m2,m3,ia,iz,ja,jz,ia_1,ja_1,iz1,jz1,jd,idiffk,i,j,k,ka
real, dimension(m2,m3) :: flpu,flpv,flpw

real, dimension(m1,m2,m3) :: up,vp,wp,vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df &
   ,vt3dg,vt3dh,vt3di,vt3dn,scr2

real, dimension(m2,m3) :: dxm,dxt,dxu,dxv,dym,dyt,dyu,dyv
real, dimension(*) :: dzm,dzt

! du/dx

do j = ja,jz
   do i = ia,iz1
      ka = max(nint(flpu(i,j)),nint(flpu(i-1,j)))
      do k = 1,ka-1
         vt3da(k,i,j) = 0.
      enddo
      do k = ka,m1
         vt3da(k,i,j) = (up(k,i,j) - up(k,i-1,j)) * dxt(i,j)
      enddo
         enddo
      enddo

! dv/dx

do j = ja_1,jz
   do i = ia_1,iz
      ka = max(nint(flpv(i+1,j)),nint(flpv(i,j)))
      do k = 1,ka-1
         vt3db(k,i,j) = 0.
      enddo
      do k = ka,m1
         vt3db(k,i,j) = (vp(k,i+1,j) - vp(k,i,j)) * dxm(i,j)
      enddo
   enddo
enddo

! dw/dx

do j = ja,jz
   do i = ia_1,iz
      ka = max(nint(flpw(i+1,j)),nint(flpw(i,j)))
      do k = 1,ka-1
         vt3df(k,i,j) = 0.
      enddo
      do k = ka,m1
         vt3df(k,i,j) = (wp(k,i+1,j) - wp(k,i,j)) * dxu(i,j)
      enddo
   enddo
enddo

! du/dy

do j = ja_1,jz
   do i = ia_1,iz
      ka = max(nint(flpu(i,j+1)),nint(flpu(i,j)))
      do k = 1,ka-1
         vt3dn(k,i,j) = 0.
      enddo
      do k = ka,m1
         vt3dn(k,i,j) = (up(k,i,j+1) - up(k,i,j)) * dym(i,j)
      enddo
   enddo
enddo

! dv/dy

do j = ja,jz1
   do i = ia,iz
      ka = max(nint(flpv(i,j)),nint(flpv(i,j-1)))
      do k = 1,ka-1
         vt3dc(k,i,j) = 0.
      enddo
      do k = ka,m1
         vt3dc(k,i,j) = (vp(k,i,j) - vp(k,i,j-1)) * dyt(i,j)
      enddo
   enddo
enddo

! dw/dy

do j = ja_1,jz
   do i = ia,iz
      ka = max(nint(flpw(i,j+1)),nint(flpw(i,j)))
      do k = 1,ka-1
         vt3dg(k,i,j) = 0.
      enddo
      do k = ka,m1
         vt3dg(k,i,j) = (wp(k,i,j+1) - wp(k,i,j)) * dyv(i,j)
      enddo
   enddo
enddo

! du/dz

do j = ja,jz
   do i = ia_1,iz
      ka = nint(flpu(i,j))
      do k = 1,ka-1
         vt3dd(k,i,j) = 0.
      enddo
      do k = ka,m1-1
         vt3dd(k,i,j) = (up(k+1,i,j) - up(k,i,j)) * dzm(k)
      enddo
   enddo
enddo

! dv/dz

do j = ja_1,jz
   do i = ia,iz
      ka = nint(flpv(i,j))
      do k = 1,ka-1
         vt3de(k,i,j) = 0.
      enddo
      do k = ka,m1-1
         vt3de(k,i,j) = (vp(k+1,i,j) - vp(k,i,j)) * dzm(k)
      enddo
   enddo
enddo

! dw/dz

if (idiffk .ge. 3 .and. idiffk /= 7) then
   do j = ja,jz
      do i = ia,iz
         ka = nint(flpw(i,j))
         do k = 1,ka
            scr2(k,i,j) = 0.
         enddo
         do k = ka+1,m1
            scr2(k,i,j) = (wp(k,i,j) - wp(k-1,i,j)) * dzt(k)
         enddo
      enddo
   enddo
endif

if (idiffk .le. 2 .or. idiffk == 7) then
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            vt3dh(k,i,j) =2. * (vt3da(k,i,j) * vt3da(k,i,j)  &
               + vt3dc(k,i,j) * vt3dc(k,i,j))                &
               + .0625 * (vt3db(k,i,j) + vt3db(k,i-1,j)      &
               + vt3db(k,i,j-jd) + vt3db(k,i-1,j-jd)         &
               + vt3dn(k,i,j) + vt3dn(k,i-1,j)               &
               + vt3dn(k,i,j-jd) + vt3dn(k,i-1,j-jd)) ** 2
            vt3di(k,i,j) = .0625 * ((vt3dd(k,i,j) + vt3dd(k-1,i,j)  &
               + vt3dd(k,i-1,j) + vt3dd(k-1,i-1,j)) ** 2            &
               + (vt3de(k,i,j) + vt3de(k-1,i,j)                     &
               + vt3de(k,i,j-jd) + vt3de(k-1,i,j-jd)) ** 2)
         enddo
      enddo
   enddo
else
   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            vt3da(k,i,j) = 2. * vt3da(k,i,j)
            vt3dc(k,i,j) = 2. * vt3dc(k,i,j)
            scr2(k,i,j) = 2. * scr2(k,i,j)
            vt3db(k,i,j) = vt3db(k,i,j) + vt3dn(k,i,j)
            vt3dn(k,i,j) = vt3db(k,i,j)
            vt3dd(k,i,j) = vt3dd(k,i,j) + vt3df(k,i,j)
            vt3de(k,i,j) = vt3de(k,i,j) + vt3dg(k,i,j)
            vt3di(k,i,j) = 0.333333  &
               * (vt3da(k,i,j) + vt3dc(k,i,j) + scr2(k,i,j))
         enddo
      enddo

      do k = 2,m1-1
         vt3da(k,iz1,j) = 2. * vt3da(k,iz1,j)
         vt3db(k,ia_1,j) = vt3db(k,ia_1,j) + vt3dn(k,ia_1,j)
         vt3dn(k,ia_1,j) = vt3db(k,ia_1,j)
         vt3dd(k,ia_1,j) = vt3dd(k,ia_1,j) + vt3df(k,ia_1,j)
      enddo
   enddo

   do i = ia_1,iz
      do k = 2,m1-1
         vt3dc(k,i,jz1) = 2. * vt3dc(k,i,jz1)
         vt3db(k,i,ja_1) = vt3db(k,i,ja_1) + vt3dn(k,i,ja_1)
         vt3dn(k,i,ja_1) = vt3db(k,i,ja_1)
         vt3de(k,i,ja_1) = vt3de(k,i,ja_1) + vt3dg(k,i,ja_1)
      enddo
   enddo

   do j = ja,jz
      do i = ia,iz
         do k = 2,m1-1
            vt3dh(k,i,j) = .5 * (  &
               (vt3da(k,i,j) - vt3di(k,i,j)) ** 2           &
               + (vt3dc(k,i,j) - vt3di(k,i,j)) ** 2         &
               + ( scr2(k,i,j) - vt3di(k,i,j)) ** 2)        &
               + .0625 * ((vt3db(k,i,j) + vt3db(k,i-1,j)    &
               + vt3db(k,i,j-jd) + vt3db(k,i-1,j-jd)) ** 2  &
               + (vt3dd(k,i,j) + vt3dd(k,i-1,j)             &
               + vt3dd(k-1,i,j) + vt3dd(k-1,i-1,j)) ** 2    &
               + (vt3de(k,i,j) + vt3de(k-1,i,j)             &
               + vt3de(k,i,j-jd) + vt3de(k-1,i,j-jd)) ** 2)
            vt3di(k,i,j) = vt3dh(k,i,j)
         enddo
      enddo
   enddo
endif

return
end
