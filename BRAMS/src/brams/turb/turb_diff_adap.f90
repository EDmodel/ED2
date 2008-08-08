!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine diffvel_adap(m1,m2,m3,ia,iz,ja,jz,jd  &
   ,iz1,jz1,izu,jzv,idiffkk  &
   ,up,vp,wp,ut,vt,wt,vt3da,vt3db,vt3dc  &
   ,vt3dd,vt3de,vt3df,vt3dg,vt3dj,vt3dk  &
   ,vt3dl,vt3dm,vt3dn,vt3do  &
   ,aru,arv,arw,volu,volv,volw,flpu,flpv,flpw  &
   ,sflux_u,sflux_v,sflux_w,dn0,dn0u,dn0v,scr1,scr2,topma,ibcon,mynum)

use mem_grid
use mem_scratch

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,jd,iz1,jz1  &
   ,izu,jzv,idiffkk,ibcon,mynum,i,j,k,ka
real, dimension(m2,m3) :: flpu,flpv,flpw

real :: akn,ako,akp,cross,c1,c2,dtlvi,volui,volvi,volwi,uw,vw,ww,wt1,topma_t
real, dimension(m1,m2,m3) :: up,vp,wp,ut,vt,wt,vt3da,vt3db,vt3dc,vt3dd  &
   ,vt3de,vt3df,vt3dg,vt3dj,vt3dk,vt3dl,vt3dm,vt3dn,vt3do  &
   ,aru,arv,arw,volu,volv,volw  &
   ,scr1,scr2,dn0,dn0u,dn0v
real, dimension(m2,m3) :: sflux_u,sflux_v,sflux_w,topma

integer :: k1,k2

!          compute fluxes with k's and previous strains

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
         akn = .25 * (scr2(k,i  ,j   ) + scr2(k+1,i  ,j   )   &
                    + scr2(k,i  ,j+jd) + scr2(k+1,i  ,j+jd))
         ako = .25 * (scr2(k,i  ,j   ) + scr2(k+1,i  ,j   )   &
                    + scr2(k,i+1,j   ) + scr2(k+1,i+1,j   ))
         akp = .25 * (scr2(k,i  ,j   ) + scr2(k  ,i+1,j   )   &
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

cross = 0.
if (idiffkk >= 3 .and. idiffkk /= 7) cross = 1.

!------------------------ vertical u-component diffusion ----------------

do j = ja,jz
   do i = ia,izu
      k2 = nint(flpu(i,j))
      k1 = k2-1
      do k = k1,m1-1
         vctr1(k) = cross * vt3df(k,i,j)
         vctr2(k) = dzm(k) * (scr1(k,i,j) + scr1(k+1,i,j)  &
                  + scr1(k,i+1,j) + scr1(k+1,i+1,j))
      enddo

      if (nstbot == 0) then
         vctr1(1) = vctr1(1) + up(1,i,j) * dzm(1) * .25  &
            * (scr1(1,i,j) + scr1(2,i,j) + scr1(1,i+1,j) + scr1(2,i+1,j))
      else
         uw = sflux_u(i,j) * arw(m1-2,i,j) + sflux_u(i+1,j) * arw(m1-2,i+1,j)
      endif

      if (nsttop == 1) then
         vctr1(m1-1) = 0.
      else
         vctr1(m1-1) = vctr1(m1-1) - up(m1-1,i,j) * dzm(m1-1) * .25  &
            * (scr1(m1,i,j) + scr1(m1-1,i,j)  &
            + scr1(m1,i+1,j) + scr1(m1-1,i+1,j))
      endif

      c1 = .5 * dtlv
      c2 = .125 * dtlv
      do k = k2,m1-1
         volui = 1. / volu(k,i,j)
         vt3dl(k,i,j) = up(k,i,j) * dn0u(k,i,j) + c1 * volui  &
            * (vctr1(k-1) * (arw(k-1,i,j) + arw(k-1,i+1,j))  &
            - vctr1(k) * (arw(k,i,j) + arw(k,i+1,j)))
         vt3dj(k,i,j) = -c2 * vctr2(k-1)  &
                      * (arw(k-1,i,j) + arw(k-1,i+1,j)) * volui
         vt3dk(k,i,j) = -c2 * vctr2(k)  &
                      * (arw(k,i,j) + arw(k,i+1,j)) * volui
         vt3do(k,i,j) = dn0u(k,i,j) - vt3dj(k,i,j) - vt3dk(k,i,j)

! Surface flux representation for normal topography.  Will require
! modification for vertical or overhanging slopes or for resolved canopy flow.

         if (nstbot == 1) then
            if (k == k2) then    
               wt1 = aru(k2,i,j) / aru(k2+1,i,j)
               vt3dl(k,i,j) = vt3dl(k,i,j) + c1 * volui * uw * wt1
            elseif (k == k2+1) then             
               vt3dl(k,i,j) = vt3dl(k,i,j) + c1 * volui * uw * (1. - wt1)
            endif
         endif
      enddo
        
      vt3dj(k2,i,j) = 0.
      vt3dk(m1-1,i,j) = 0.
      if (nstbot == 1) vt3do(k2,i,j) = dn0u(k2,i,j) - vt3dk(k2,i,j)
      if (nsttop == 1) vt3do(m1-1,i,j) = dn0u(m1-1,i,j) - vt3dj(m1-1,i,j)

   enddo
enddo

call tridiff1_adap(m1,m2,m3,ia,izu,ja,jz,m1-1,flpu  &
   ,vt3dj,vt3do,vt3dk,vt3dl,vt3do,vt3dk)

!---------------------horizontal u-component diffusion ----------------

dtlvi = 1.0 / dtlv

if (jd == 1) then
   do j = ja,jz
      do i = ia,izu
         ka = nint(flpu(i,j))
         do k = ka,m1-1
            ut(k,i,j) = ut(k,i,j) + dtlvi * (vt3do(k,i,j) - up(k,i,j))  &
               + ((aru(k,i,j) + aru(k,i-1,j)) * vt3da(k,i,j)            &
               -  (aru(k,i,j) + aru(k,i+1,j)) * vt3da(k,i+1,j)          &
               +  (arv(k,i,j-1) + arv(k,i+1,j-1)) * vt3dn(k,i,j-1)      &
               -  (arv(k,i,j) + arv(k,i+1,j)) * vt3dn(k,i,j) )          &
               / (2. * volu(k,i,j) * dn0u(k,i,j))
         enddo
      enddo
   enddo
else
   do j = ja,jz
      do i = ia,izu
         ka = nint(flpu(i,j))
         do k = ka,m1-1
            ut(k,i,j) = ut(k,i,j) + dtlvi * (vt3do(k,i,j) - up(k,i,j))  &
               + ((aru(k,i,j) + aru(k,i-1,j)) * vt3da(k,i,j)            &
               -  (aru(k,i,j) + aru(k,i+1,j)) * vt3da(k,i+1,j) )        &
               / (2. * volu(k,i,j) * dn0u(k,i,j))
         enddo
      enddo
   enddo
endif

!------------------------ vertical v-component diffusion ----------------

if (jd == 0) go to 99

do j = ja,jzv
   do i = ia,iz
      k2 = nint(flpv(i,j))
      k1 = k2-1
      do k = k1,m1-1
         vctr1(k) = cross * vt3dg(k,i,j)
         vctr2(k) = dzm(k) * (scr1(k,i,j) + scr1(k+1,i,j)  &
            + scr1(k,i,j+jd) + scr1(k+1,i,j+jd))
      enddo

      if (nstbot == 0) then
         vctr1(1) = vctr1(1) + vp(1,i,j) * dzm(1) * .25  &
            * (scr1(1,i,j) + scr1(2,i,j) + scr1(1,i,j+jd) + scr1(2,i,j+jd))
      else
         vw = sflux_v(i,j) * arw(m1-2,i,j) + sflux_v(i,j+jd) * arw(m1-2,i,j+jd)   
      endif

      if (nsttop == 1) then
         vctr1(m1-1) = 0.
      else
         vctr1(m1-1) = vctr1(m1-1) - vp(m1-1,i,j)  &
            * .25 * (scr1(m1,i,j) + scr1(m1-1,i,j)  &
            + scr1(m1,i,j+jd) + scr1(m1-1,i,j+jd)) * dzm(m1-1)
      endif

      c1 = .5 * dtlv
      c2 = .125 * dtlv
      do k = k2,m1-1
         volvi = 1. / volv(k,i,j)
         vt3dm(k,i,j) = vp(k,i,j) * dn0v(k,i,j) + c1 * volvi  &
            * (vctr1(k-1) * (arw(k-1,i,j) + arw(k-1,i,j+1))  &
            - vctr1(k) * (arw(k,i,j) + arw(k,i,j+1)))
         vt3dj(k,i,j) = -c2 * vctr2(k-1)  &
                      * (arw(k-1,i,j) + arw(k-1,i,j+1)) * volvi
         vt3dk(k,i,j) = -c2 * vctr2(k)  &
                      * (arw(k,i,j) + arw(k,i,j+1)) * volvi
         vt3do(k,i,j) = dn0v(k,i,j)-vt3dj(k,i,j)-vt3dk(k,i,j)

! Surface flux representation for normal topography.  Will require
! modification for vertical or overhanging slopes or for resolved canopy flow.

         if (nstbot == 1) then
            if (k == k2) then    
               wt1 = arv(k2,i,j) / arv(k2+1,i,j)
               vt3dm(k,i,j) = vt3dm(k,i,j) + c1 * volvi * vw * wt1
            elseif (k == k2+1) then             
               vt3dm(k,i,j) = vt3dm(k,i,j) + c1 * volvi * vw * (1. - wt1)
            endif
         endif

      enddo

      vt3dj(k2,i,j) = 0.
      vt3dk(m1-1,i,j) = 0.
      if (nstbot == 1) vt3do(k2,i,j) = dn0v(k2,i,j) - vt3dk(k2,i,j)
      if (nsttop == 1) vt3do(m1-1,i,j) = dn0v(m1-1,i,j) - vt3dj(m1-1,i,j)
   enddo
enddo

call tridiff1_adap(m1,m2,m3,ia,iz,ja,jzv,m1-1,flpv  &
   ,vt3dj,vt3do,vt3dk,vt3dm,vt3do,vt3dk)

!--------------------- horizontal v-component diffusion ----------------

do j = ja,jzv
   do i = ia,iz
      ka = nint(flpv(i,j))
      do k = ka,m1-1
         vt(k,i,j) = vt(k,i,j) + dtlvi * (vt3do(k,i,j) - vp(k,i,j))  &
            + ((aru(k,i-1,j) + aru(k,i-1,j+1)) * vt3db(k,i-1,j)      &
            -  (aru(k,i,j) + aru(k,i,j+1)) * vt3db(k,i,j)            &
            +  (arv(k,i,j) + arv(k,i,j-jd)) * vt3dc(k,i,j)           &
            -  (arv(k,i,j) + arv(k,i,j+jd)) * vt3dc(k,i,j+jd))       &
            / (2. * volv(k,i,j) * dn0v(k,i,j))
      enddo
   enddo
enddo

99   continue

!------------------------ vertical w-component diffusion ----------------

do j = ja,jz
   do i = ia,iz
      ka = nint(flpw(i,j))
      do k = ka,m1-1
         vctr1(k) = 0.
         vctr2(k) = dzt(k) * scr1(k,i,j)
      enddo

      if (nstbot == 0) then
         vctr1(2) = (1.0 + cross) * wp(1,i,j) * scr1(2,i,j) * dzt(2)
      else
         ww = sflux_w(i,j) * arw(m1-2,i,j)
      endif
      
      if (nsttop == 1) then
         vctr1(m1-1) = 0.
      else
         vctr1(m1-1) = -(1.0 + cross) * wp(m1-1,i,j)  &
            * scr1(m1-1,i,j) * dzt(m1-1)
      endif

      c1 = .5 * dtlv
      c2 = 2. * dtlv
      do k = ka,m1-2
         volwi = 1. / volw(k,i,j)
         vt3dn(k,i,j) = wp(k,i,j) * (dn0(k,i,j) + dn0(k+1,i,j)) * .5  &
            + c1 * volwi  &
            * (vctr1(k) * (arw(k-1,i,j) + arw(k,i,j))  &
             - vctr1(k+1) * (arw(k+1,i,j) + arw(k,i,j)))
         vt3dj(k,i,j) = -c2 * vctr2(k)  &
                      * (arw(k-1,i,j) + arw(k,i,j)) * volwi
         vt3dk(k,i,j) = -c2 * vctr2(k+1)  &
                      * (arw(k+1,i,j) + arw(k,i,j)) * volwi
         vt3do(k,i,j) = (dn0(k,i,j) + dn0(k+1,i,j)) * .5  &
             - vt3dj(k,i,j) - vt3dk(k,i,j)

! Surface flux representation for normal topography.  Will require
! modification for vertical or overhanging slopes or for resolved canopy flow.

         if (nstbot == 1) then
            if (k == k2) then    
               topma_t = .25 * (topma(i,j) + topma(i-1,j)  &
                  + topma(i,j-jdim) + topma(i-1,j-jdim))
               wt1 = (zm(k2) - topma_t) * dzt(k2)

               vt3dn(k,i,j) = vt3dn(k,i,j) + dtlv * volwi * ww * wt1
            elseif (k == k2+1) then
               vt3dn(k,i,j) = vt3dn(k,i,j) + dtlv * volwi * ww * (1. - wt1)
            endif
         endif

      enddo

      vt3dj(ka,i,j) = 0.
      vt3dk(m1-2,i,j) = 0.
      if (nstbot == 1) vt3do(ka,i,j) =  &
           (dn0(ka,i,j) + dn0(ka+1,i,j)) * .5 - vt3dk(ka,i,j)
      if (nsttop == 1) vt3do(m1-2,i,j) =  &
           (dn0(m1-2,i,j) + dn0(m1-1,i,j)) * .5 - vt3dj(m1-2,i,j)

   enddo
enddo

call tridiff1_adap(m1,m2,m3,ia,iz,ja,jz,m1-2,flpw  &
   ,vt3dj,vt3do,vt3dk,vt3dn,vt3do,vt3dk)

!--------------------- horizontal w-component diffusion ----------------

if (idiffkk >= 3 .and. idiffkk /= 7) then

   if (jd == 1) then
      do j = ja,jz
        do i = ia,iz
           ka = nint(flpw(i,j))
           do k = ka,m1-2
              wt(k,i,j) = wt(k,i,j) + dtlvi * (vt3do(k,i,j) - wp(k,i,j))  &
                 + ((aru(k+1,i-1,j) + aru(k,i-1,j)) * vt3dd(k,i-1,j)      &
                 -  (aru(k+1,i,j) + aru(k,i,j)) * vt3dd(k,i,j)            &
                 +  (arv(k+1,i,j-1) + arv(k,i,j-1)) * vt3de(k,i,j-1)      &
                 -  (arv(k+1,i,j) + arv(k,i,j)) * vt3de(k,i,j))           &
                 / (volw(k,i,j) * (dn0(k,i,j) + dn0(k+1,i,j)))
           enddo
         enddo
      enddo
   else
      do j = ja,jz
         do i = ia,iz
            ka = nint(flpw(i,j))
            do k = ka,m1-2
               wt(k,i,j) = wt(k,i,j) + dtlvi * (vt3do(k,i,j) - wp(k,i,j))  &
                  + ((aru(k+1,i-1,j) + aru(k,i-1,j)) * vt3dd(k,i-1,j)      &
                  -  (aru(k+1,i,j) + aru(k,i,j)) * vt3dd(k,i,j))           &
                  / (volw(k,i,j) * (dn0(k,i,j) + dn0(k+1,i,j)))
            enddo
         enddo
      enddo
   endif

else

   if (jd == 1) then
      do j = ja,jz
        do i = ia,iz
           ka = nint(flpw(i,j))
           do k = ka,m1-2
              wt(k,i,j) = wt(k,i,j) + dtlvi * (vt3do(k,i,j) - wp(k,i,j))  &
                 + ((aru(k+1,i-1,j) + aru(k,i-1,j)) * vt3df(k,i-1,j)      &
                 -  (aru(k+1,i,j) + aru(k,i,j)) * vt3df(k,i,j)            &
                 +  (arv(k+1,i,j-1) + arv(k,i,j-1)) * vt3dg(k,i,j-1)      &
                 -  (arv(k+1,i,j) + arv(k,i,j)) * vt3dg(k,i,j))           &
                 / (volw(k,i,j) * (dn0(k,i,j) + dn0(k+1,i,j)))
           enddo
         enddo
      enddo
   else
      do j = ja,jz
         do i = ia,iz
            ka = nint(flpw(i,j))
            do k = ka,m1-2
               wt(k,i,j) = wt(k,i,j) + dtlvi * (vt3do(k,i,j) - wp(k,i,j))  &
                  + ((aru(k+1,i-1,j) + aru(k,i-1,j)) * vt3df(k,i-1,j)      &
                  -  (aru(k+1,i,j) + aru(k,i,j)) * vt3df(k,i,j))           &
                  / (volw(k,i,j) * (dn0(k,i,j) + dn0(k+1,i,j)))
            enddo
         enddo
      enddo
   endif

endif
return
end

!     ***************************************************************

subroutine diffsclr_adap(m1,m2,m3,ia,iz,ja,jz,jd,n,ksf,flpw  &
   ,scp,sct,vt3da,vt3dc,vt3df,vt3dg,vt3dj,vt3dk,vt3dl,vt3dm,vt3do  &
   ,sfcflx,vt2db,dn0,vkkh,hkkh,aru,arv,arw,volt,dn0volti,dxu,dyv,topma)

use mem_grid
use mem_scratch

implicit none

integer :: n,ksf

integer :: m1,m2,m3,ia,iz,ja,jz,jd,i,j,k,ka
real, dimension(m2,m3) :: flpw

integer, save :: ksf_save = 0
real :: c1,dtlti,c1volti,hdxu,hdyv,tw,wt1,topma_t
real, dimension(m1,m2,m3) :: scp,sct,vt3da,vt3dc,vt3df,vt3dg,vt3dj  &
   ,vt3dk,vt3dl,vt3dm,vt3do,dn0,vkkh,hkkh,aru,arv,arw,volt  &
   ,dn0volti

real, dimension(m2,m3) :: sfcflx,vt2db,dxu,dyv,topma

!  compute vertical diffusion matrix coefficients for scalars

if (n == 1 .or. ksf /= ksf_save) then
   ksf_save = ksf
   c1 = .5 * dtlt
   do j = ja,jz
      do i = ia,iz
         ka = nint(flpw(i,j))
         do k = ka-1,m1-1
            vctr1(k) = dzm(k) * (vkkh(k,i,j) + vkkh(k+1,i,j))
         enddo
         do k = ka,m1-1
            dn0volti(k,i,j) = 1. / (dn0(k,i,j) * volt(k,i,j))
            c1volti = c1 * dn0(k,i,j) * dn0volti(k,i,j)

            vt3dj(k,i,j) = - vctr1(k-1) * arw(k-1,i,j) * c1volti
            vt3dk(k,i,j) = - vctr1(k) * arw(k,i,j) * c1volti
            vt3do(k,i,j) = dn0(k,i,j) - vt3dj(k,i,j) - vt3dk(k,i,j)
         enddo
         vt3dj(ka,i,j) = 0.
         vt3dk(m1-1,i,j) = 0.
         if (nstbot == 1) vt3do(ka,i,j) = dn0(ka,i,j) - vt3dk(ka,i,j)
         if (nsttop == 1) vt3do(m1-1,i,j) = dn0(m1-1,i,j) - vt3dj(m1-1,i,j)

! Bob (10/27/00):  Evaluate ci2i, cjp1, and cji terms here 
! (vt2db, vt3dl, vt3dm) for new inline partial tridiff.  
! With this, no longer need to pass ci or cip1 (vt3do, vt3dk) to tridiff2.

         vt2db(i,j) = 1. / vt3do(ka,i,j)
         vt3dl(ka,i,j) = vt3dk(ka,i,j) * vt2db(i,j)
         do k = ka+1,m1-1
            vt3dm(k,i,j) = 1. / (vt3do(k,i,j) - vt3dj(k,i,j) * vt3dl(k-1,i,j))
            vt3dl(k,i,j) = vt3dk(k,i,j) * vt3dm(k,i,j)
         enddo
      enddo
   enddo

endif

!         finish matrix coefficients

do j = ja,jz
   do i = ia,iz
      ka = nint(flpw(i,j))

      if (nstbot == 1) tw = sfcflx(i,j) * arw(m1-2,i,j)

      do k = ka,m1-1
         vt3dc(k,i,j) = scp(k,i,j) * dn0(k,i,j)

! Surface flux representation for normal topography.  Will require
! modification for vertical or overhanging slopes or for resolved canopy flow.

         if (nstbot == 1) then
            if (k == ka) then    
               topma_t = .25 * (topma(i,j) + topma(i-1,j)  &
                  + topma(i,j-jdim) + topma(i-1,j-jdim))
                  
! Bob (2/1/02) set minimum value on wt1 of .1                 
                  
               wt1 = max(.1,(zm(ka) - topma_t) * dzt(ka))

               vt3dc(k,i,j) = vt3dc(k,i,j) + dtlt / volt(k,i,j) * tw * wt1
            elseif (k == ka+1) then             
               vt3dc(k,i,j) = vt3dc(k,i,j) + dtlt / volt(k,i,j) * tw  &
                  * (1. - wt1)                  
            endif
         endif

      enddo

      if (nstbot == 0) then
         vt3dc(ka,i,j) = vt3dc(ka,i,j)  &
            + c1 * (vkkh(ka-1,i,j) + vkkh(ka,i,j))  &
            * scp(ka-1,i,j) * dzm(ka-1) * arw(ka-1,i,j) / volt(ka,i,j)
      endif

      if (nsttop == 0) then
         vt3dc(m1-1,i,j) = vt3dc(m1-1,i,j)  &
            - c1 * (vkkh(m1-1,i,j) + vkkh(m1,i,j))  &
            * scp(m1,i,j) * dzm(m1-1) * arw(m1-1,i,j) / volt(m1-1,i,j)
      endif
   enddo
enddo

!  Finish tridiff solution

do j = ja,jz
   do i = ia,iz

     ka = nint(flpw(i,j))
     vt3da(ka,i,j) = vt3dc(ka,i,j) * vt2db(i,j)

      do k = ka+1,m1-1
         vt3da(k,i,j) = (vt3dc(k,i,j) - vt3dj(k,i,j) * vt3da(k-1,i,j))  &
            * vt3dm(k,i,j)
      enddo

      do k = m1-2,ka,-1
         vt3da(k,i,j) = vt3da(k,i,j) - vt3dl(k,i,j) * vt3da(k+1,i,j)
      enddo

   enddo
enddo


!     compute 2 horizontal scalar gradients needed for dscp/dt

! -K dT/dx

do j = ja,jz
   do i = 1,iz
      ka = max(nint(flpw(i+1,j)),nint(flpw(i,j)))
      do k = 1,ka-1
         vt3df(k,i,j) = 0.
      enddo
      hdxu = .5 * dxu(i,j)
      do k = ka,m1
         vt3df(k,i,j) = (hkkh(k,i,j) + hkkh(k,i+1,j))  &
                      * (scp(k,i,j) - scp(k,i+1,j)) * hdxu
      enddo
   enddo
enddo

!  -K dT/dy

if (jd == 1) then

   do j = 1,jz
      do i = ia,iz
         ka = max(nint(flpw(i,j+1)),nint(flpw(i,j)))
         do k = 1,ka-1
            vt3dg(k,i,j) = 0.
         enddo
         hdyv = .5 * dyv(i,j)
         do k = ka,m1
            vt3dg(k,i,j) = (hkkh(k,i,j) + hkkh(k,i,j+jd))  &
                         * (scp(k,i,j) - scp(k,i,j+jd)) * hdyv
         enddo
      enddo
   enddo

endif

!   horizontal flux divergence and tendency update

dtlti = 1.0 / dtlt
 
if (jd == 1) then

   do j = ja,jz
      do i = ia,iz
         ka = nint(flpw(i,j))
         do k = ka,m1-1
            sct(k,i,j) = sct(k,i,j) + (vt3da(k,i,j) - scp(k,i,j)) * dtlti  &
               + (aru(k,i-1,j) * vt3df(k,i-1,j)   &
               -  aru(k,i,j) * vt3df(k,i,j)       &
               +  arv(k,i,j-1) * vt3dg(k,i,j-1)   &
               -  arv(k,i,j) * vt3dg(k,i,j))      &
               / (volt(k,i,j) * dn0(k,i,j))
         enddo
      enddo
   enddo

else

   do j = ja,jz
      do i = ia,iz
         ka = nint(flpw(i,j))
         do k = ka,m1-1
            sct(k,i,j) = sct(k,i,j) + (vt3da(k,i,j) - scp(k,i,j)) * dtlti  &
               + (aru(k,i-1,j) * vt3df(k,i-1,j)  &
               -  aru(k,i,j) * vt3df(k,i,j))     &
               / (volt(k,i,j) * dn0(k,i,j))
         enddo
      enddo
   enddo

endif

return
end

!     ******************************************************************

subroutine tridiff1_adap(m1,m2,m3,ia,iz,ja,jz,kz,flpx  &
   ,cim1,ci,cip1,rhs,cj,cjp1)

implicit none
integer :: m1,m2,m3,ia,iz,ja,jz,kz,i,j,k,ka
real :: flpx(m2,m3)

real cji
real, dimension(m1,m2,m3) :: cim1,ci,cip1,rhs,cj,cjp1

do j = ja,jz
   do i = ia,iz
      ka = nint(flpx(i,j))

      cjp1(ka,i,j) = cip1(ka,i,j) / ci(ka,i,j)
      rhs(ka,i,j) = rhs(ka,i,j) / ci(ka,i,j)

      do k = ka+1,kz
         cj(k,i,j) = ci(k,i,j) - cim1(k,i,j) * cjp1(k-1,i,j)
         cji = 1. / cj(k,i,j)
         cjp1(k,i,j) = cip1(k,i,j) * cji
         rhs(k,i,j) = (rhs(k,i,j) - cim1(k,i,j) * rhs(k-1,i,j)) * cji
      enddo

      cj(kz,i,j) = rhs(kz,i,j)

      do k = kz-1,ka,-1
         cj(k,i,j) = rhs(k,i,j) - cjp1(k,i,j) * cj(k+1,i,j)
      enddo

   enddo
enddo
return
end
