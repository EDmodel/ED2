!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine vel_advectc_adap(m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,flpu,flpv,flpw  &
   ,uc,vc,wc,ut,vt,wt,dn0,dn0u,dn0v  &
   ,aru,arv,arw,volu,volv,volw,flxu,flxv,flxw,ucb,vcb,wcb,time)

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim
real, dimension(m1,m2,m3) :: uc,vc,wc,ut,vt,wt,dn0,dn0u,dn0v  &
   ,aru,arv,arw,volu,volv,volw,flxu,flxv,flxw,ucb,vcb,wcb

integer :: j,i,k,jm,im
real, dimension(m2,m3) :: flpu,flpv,flpw
real :: fx1,fx2,fy1,fy2,fz1,fz2,cx1,cx2,cy1,cy2,cz1,cz2,time
integer :: lpu,lpv,lpw

! Compute fluxes at velocity points and initialize boundary 
! advected velocities

do j = 1,m3
   do i = 1,m2
      do k = 1,m1-1
         flxu(k,i,j) = uc(k,i,j) * aru(k,i,j) * dn0u(k,i,j)
         flxv(k,i,j) = vc(k,i,j) * arv(k,i,j) * dn0v(k,i,j)
         flxw(k,i,j) = wc(k,i,j) * arw(k,i,j)  &
                     * .5 * (dn0(k,i,j) + dn0(k+1,i,j))
         ucb(k,i,j) = uc(k,i,j)
         vcb(k,i,j) = vc(k,i,j)
         wcb(k,i,j) = wc(k,i,j)
      enddo
   enddo
enddo

! Compute advection contribution to U tendency

do j = ja,jz
   do i = ia,izu

      call advected (m1,m2,m3,uc,vc,wc,ucb,vcb,wcb,flpu,flpv,flpw,i,j,jdim)
      
      lpu = nint(flpu(i,j))

      do k = lpu,m1-1

         fx1 = flxu(k,i,j)      + flxu(k,i-1,j)
         fx2 = flxu(k,i,j)      + flxu(k,i+1,j)
         fy1 = flxv(k,i,j-jdim) + flxv(k,i+1,j-jdim)
         fy2 = flxv(k,i,j)      + flxv(k,i+1,j)
         fz1 = flxw(k-1,i,j)    + flxw(k-1,i+1,j)
         fz2 = flxw(k,i,j)      + flxw(k,i+1,j)

! The following c values are for first order fluxes

         cx1 = sign(1.,fx1)
         cx2 = sign(1.,fx2)
         cy1 = sign(1.,fy1)
         cy2 = sign(1.,fy2)
         cz1 = sign(1.,fz1)
         cz2 = sign(1.,fz2)

         ut(k,i,j) = ut(k,i,j)                       &
            + .25 / (volu(k,i,j) * dn0u(k,i,j)) * (  &
! ut flux terms
            - fx2 * (uc(k,i,j) + ucb(k,i+1,j)        &
            + cx2 * (uc(k,i,j) - ucb(k,i+1,j)))      &

            + fx1 * (uc(k,i,j) + ucb(k,i-1,j)        &
            + cx1 * (ucb(k,i-1,j) - uc(k,i,j)))      &

            - fy2 * (uc(k,i,j) + ucb(k,i,j+jdim)     &
            + cy2 * (uc(k,i,j) - ucb(k,i,j+jdim)))   &

            + fy1 * (uc(k,i,j) + ucb(k,i,j-jdim)     &
            + cy1 * (ucb(k,i,j-jdim) - uc(k,i,j)))   &

            - fz2 * (uc(k,i,j) + ucb(k+1,i,j)        &
            + cz2 * (uc(k,i,j) - ucb(k+1,i,j)))      &

            + fz1 * (uc(k,i,j) + ucb(k-1,i,j)        &
            + cz1 * (ucb(k-1,i,j) - uc(k,i,j)))      &
! ut divergence terms
            + (fx2 - fx1 + fy2 - fy1 + fz2 - fz1)    &
            * 2. * uc(k,i,j)                         )
      enddo
   enddo
enddo

! Compute advection contribution to V tendency

do j = ja,jzv
   do i = ia,iz
      lpv = nint(flpv(i,j))
      do k = lpv,m1-1
         fx1 = flxu(k,i-1,j) + flxu(k,i-1,j+jdim)
         fx2 = flxu(k,i,j)   + flxu(k,i,j+jdim)
         fy1 = flxv(k,i,j)   + flxv(k,i,j-jdim)
         fy2 = flxv(k,i,j)   + flxv(k,i,j+jdim)
         fz1 = flxw(k-1,i,j) + flxw(k-1,i,j+jdim)
         fz2 = flxw(k,i,j)   + flxw(k,i,j+jdim)

! The following c values are for first order fluxes

         cx1 = sign(1.,fx1)
         cx2 = sign(1.,fx2)
         cy1 = sign(1.,fy1)
         cy2 = sign(1.,fy2)
         cz1 = sign(1.,fz1)
         cz2 = sign(1.,fz2)

         vt(k,i,j) = vt(k,i,j)                       &
            + .25 / (volv(k,i,j) * dn0v(k,i,j)) * (  &
! vt flux terms
            - fx2 * (vc(k,i,j) + vcb(k,i+1,j)        &
            + cx2 * (vc(k,i,j) - vcb(k,i+1,j)))      &

            + fx1 * (vc(k,i,j) + vcb(k,i-1,j)        &
            + cx1 * (vcb(k,i-1,j) - vc(k,i,j)))      &

            - fy2 * (vc(k,i,j) + vcb(k,i,j+jdim)     &
            + cy2 * (vc(k,i,j) - vcb(k,i,j+jdim)))   &

            + fy1 * (vc(k,i,j) + vcb(k,i,j-jdim)     &
            + cy1 * (vcb(k,i,j-jdim) - vc(k,i,j)))   &

            - fz2 * (vc(k,i,j) + vcb(k+1,i,j)        &
            + cz2 * (vc(k,i,j) - vcb(k+1,i,j)))      &

            + fz1 * (vc(k,i,j) + vcb(k-1,i,j)        &
            + cz1 * (vcb(k-1,i,j) - vc(k,i,j)))      &
! vt divergence terms
            + (fx2 - fx1 + fy2 - fy1 + fz2 - fz1)    &
            * 2.* vc(k,i,j)                          )
      enddo
   enddo
enddo

! Compute advection contribution to W tendency

do j = ja,jz
   do i = ia,iz
      lpw = nint(flpw(i,j))
      do k = lpw,m1-1
         fx1 = flxu(k,i-1,j)    + flxu(k+1,i-1,j)
         fx2 = flxu(k,i,j)      + flxu(k+1,i,j)
         fy1 = flxv(k,i,j-jdim) + flxv(k+1,i,j-jdim)
         fy2 = flxv(k,i,j)      + flxv(k+1,i,j)
         fz1 = flxw(k,i,j)      + flxw(k-1,i,j)
         fz2 = flxw(k,i,j)      + flxw(k+1,i,j)

! The following c values are for first order fluxes

         cx1 = sign(1.0,fx1)
         cx2 = sign(1.0,fx2)
         cy1 = sign(1.0,fy1)
         cy2 = sign(1.0,fy2)
         cz1 = sign(1.0,fz1)
         cz2 = sign(1.0,fz2)

         wt(k,i,j) = wt(k,i,j) + .5 / (volw(k,i,j)     &
                   * (dn0(k,i,j) + dn0(k+1,i,j))) * (  &
! wt flux terms

            - fx2 * (wc(k,i,j) + wcb(k,i+1,j)          &
            + cx2 * (wc(k,i,j) - wcb(k,i+1,j)))        &

            + fx1 * (wc(k,i,j) + wcb(k,i-1,j)          &
            + cx1 * (wcb(k,i-1,j) - wc(k,i,j)))        &

            - fy2 * (wc(k,i,j) + wcb(k,i,j+jdim)       &
            + cy2 * (wc(k,i,j) - wcb(k,i,j+jdim)))     &

            + fy1 * (wc(k,i,j) + wcb(k,i,j-jdim)       &
            + cy1 * (wcb(k,i,j-jdim) - wc(k,i,j)))     &

            - fz2 * (wc(k,i,j) + wcb(k+1,i,j)          &
            + cz2 * (wc(k,i,j) - wcb(k+1,i,j)))        &

            + fz1 * (wc(k,i,j) + wcb(k-1,i,j)          &
            + cz1 * (wcb(k-1,i,j) - wc(k,i,j)))        &

! wt divergence terms
            + (fx2 - fx1 + fy2 - fy1 + fz2 - fz1)      &
            * 2.* wc(k,i,j)                            ) 

      enddo
   enddo
enddo
return
end subroutine vel_advectc_adap

!     *********************************************************************

subroutine fa_preptc_adap(m1,m2,m3,vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df  &
   ,vt3dh,dn0,dn0u,dn0v,aru,arv,arw,volt  &
   ,dxu,dyv,dxt,dyt,zt,zm,dzm,vctr1,vctr2,jd,mynum)

implicit none

integer :: m1,m2,m3,j,i,k,im,ip,jm,jp,jd,mynum

real, dimension(m1,m2,m3) :: vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df,vt3dh  &
   ,dn0,dn0u,dn0v,aru,arv,arw,volt

real, dimension(m2,m3) :: dxu,dyv,dxt,dyt
real, dimension(*) :: zt,zm,dzm,vctr1,vctr2

! VT3DA, VT3DB, and VT3DC are input as the velocity components (averaged
! between past and current time levels) times dtlt.

! Compute Courant numbers: VT3DD, VT3DE and half Cnum VT3DF
! Compute weight at scalar point: VT3DH
! Compute vertical advective weights for the linear term: VCTR1, and VCTR2

do j = 1,m3
   do i = 1,m2
      do k = 1,m1
         vt3dd(k,i,j) = vt3da(k,i,j) * dxu(i,j)
         vt3de(k,i,j) = vt3db(k,i,j) * dyv(i,j)
         vt3df(k,i,j) = .5 * vt3dc(k,i,j) * dzm(k)
         vt3dh(k,i,j) = 1. / (volt(k,i,j) * dn0(k,i,j))
      enddo
   enddo
enddo

do k = 1,m1-1
   vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
   vctr2(k) =  (zm(k) - zt(k)) * dzm(k)
enddo

! Convert velocity components * dtlt (VT3DA, VT3DB, VT3DC)
! into mass fluxes times dtlt.

do j = 1,m3
   do i = 1,m2
      do k = 1,m1-1
         vt3da(k,i,j) = vt3da(k,i,j) * aru(k,i,j) * dn0u(k,i,j)
         vt3db(k,i,j) = vt3db(k,i,j) * arv(k,i,j) * dn0v(k,i,j)
         vt3dc(k,i,j) = vt3dc(k,i,j) * arw(k,i,j) * .5  &
                      * (dn0(k,i,j) + dn0(k+1,i,j))
      enddo
   enddo
enddo

return
end subroutine fa_preptc_adap

!     *********************************************************************

subroutine fa_xc_adap(m1,m2,m3,ia,iz,ja,jz,flpw  &
   ,scp,scr1,vt3da,vt3dd,vt3dg,vt3dh,mynum)

implicit none
integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k,mynum
real, dimension(m2,m3) :: flpw

real :: df
real, dimension(m1,m2,m3) :: scp,scr1,vt3da,vt3dd,vt3dg,vt3dh
integer :: lpw

df = .5
do j = 1,m3
   do i = 1,iz

! Compute scalar flux times dtlt [VT3DG]

      do k = 2,m1-1
         vt3dg(k,i,j) = vt3da(k,i,j)  &
            * .5 * (scr1(k,i,j) + scr1(k,i+1,j)  &
            +  vt3dd(k,i,j) * (scr1(k,i,j) - scr1(k,i+1,j)))
      enddo

          
! Modify fluxes to retain positive-definiteness on scalar quantities.
!    If a flux will remove 1/2 quantity during a timestep,
!    reduce to first order flux. This will remain positive-definite
!    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
!    both fluxes are evacuating the box.

      do k = 2,m1-1
         if (vt3da(k,i,j) > 0.) then
            if (vt3dg(k,i,j) * vt3dh(k,i,j) > df * scr1(k,i,j)) then
               vt3dg(k,i,j) = vt3da(k,i,j) * scr1(k,i,j)
            endif
         elseif (vt3da(k,i,j) < 0.) then
            if (-vt3dg(k,i,j) * vt3dh(k,i+1,j) > df * scr1(k,i+1,j)) then
               vt3dg(k,i,j) = vt3da(k,i,j) * scr1(k,i+1,j)
            endif
         endif
      enddo
   enddo
enddo

! Compute flux divergence

do j = 1,m3
   do i = ia,iz
      lpw = nint(flpw(i,j))
      do k = lpw,m1-1
         scr1(k,i,j) = scr1(k,i,j)  &
            + vt3dh(k,i,j) * (vt3dg(k,i-1,j) - vt3dg(k,i,j)  &
            + scp(k,i,j) * (vt3da(k,i,j) - vt3da(k,i-1,j)))
      enddo
   enddo
enddo

return
end subroutine fa_xc_adap

!     *********************************************************************

subroutine fa_yc_adap(m1,m2,m3,ia,iz,ja,jz,flpw  &
   ,scp,scr1,vt3db,vt3de,vt3dg,vt3dh,jdim,mynum)

implicit none
integer :: m1,m2,m3,ia,iz,ja,jz,jdim,mynum,i,j,k
real, dimension(m2,m3) :: flpw
integer :: lpw

real :: df
real, dimension(m1,m2,m3) :: scp,scr1,vt3db,vt3de,vt3dg,vt3dh

df = .5
do j = 1,jz
   do i = ia,iz

! Compute scalar flux VT3DG

      do k = 2,m1-1
         vt3dg(k,i,j) = vt3db(k,i,j)  &
            * .5 * (scr1(k,i,j) + scr1(k,i,j+jdim)  &
            +  vt3de(k,i,j) * (scr1(k,i,j) - scr1(k,i,j+jdim)))
      enddo

!      Modify fluxes to retain positive-definiteness on scalar quantities.
!         If a flux will remove 1/2 quantity during a timestep,
!         reduce to first order flux. This will remain positive-definite
!         under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
!         both fluxes are evacuating the box.

      do k = 2,m1-1
         if (vt3db(k,i,j) > 0.) then
            if (vt3dg(k,i,j) * vt3dh(k,i,j) > df * scr1(k,i,j)) then
               vt3dg(k,i,j) = vt3db(k,i,j) * scr1(k,i,j)
            endif
         elseif (vt3db(k,i,j) < 0.) then
            if (-vt3dg(k,i,j) * vt3dh(k,i,j+jdim) >   &
               df * scr1(k,i,j+jdim)) then
               vt3dg(k,i,j) = vt3db(k,i,j) * scr1(k,i,j+jdim)
            endif
         endif
      enddo
   enddo
enddo

! Compute flux divergence

do j = ja,jz
   do i = ia,iz
      lpw = nint(flpw(i,j))
      do k = lpw,m1-1
         scr1(k,i,j) = scr1(k,i,j)  &
            + vt3dh(k,i,j) * (vt3dg(k,i,j-jdim) - vt3dg(k,i,j)  &
            + scp(k,i,j) * (vt3db(k,i,j) - vt3db(k,i,j-jdim)))
      enddo
   enddo
enddo

return
end

!     *********************************************************************

subroutine fa_zc_adap(m1,m2,m3,ia,iz,ja,jz,flpw  &
   ,scp,scr1,vt3dc,vt3df,vt3dg,vt3dh,vctr1,vctr2,mynum)

implicit none
integer :: m1,m2,m3,ia,iz,ja,jz,mynum,i,j,k
real, dimension(m2,m3) :: flpw
integer :: lpw
real :: df
real, dimension(m1,m2,m3) :: scp,scr1,vt3dc,vt3df,vt3dg,vt3dh
real, dimension(*) :: vctr1,vctr2

df = .5
do j = ja,jz
   do i = ia,iz

! Compute scalar flux VT3DG
      lpw =  nint(flpw(i,j))
      do k = lpw-1,m1-1
         vt3dg(k,i,j) = vt3dc(k,i,j)  &
            * (vctr1(k) * scr1(k,i,j)  &
            +  vctr2(k) * scr1(k+1,i,j)  &
            +  vt3df(k,i,j) * (scr1(k,i,j) - scr1(k+1,i,j)))
      enddo

! Modify fluxes to retain positive-definiteness on scalar quantities.
!    If a flux will remove 1/2 quantity during a timestep,
!    reduce to first order flux. This will remain positive-definite
!    under the assumption that ABS(CFL(i)) + ABS(CFL(i-1)) < 1.0 if
!    both fluxes are evacuating the box.
      lpw = nint(flpw(i,j))
      do k = lpw-1,m1-1
         if (vt3dc(k,i,j) < 0.) then
            if (vt3dg(k,i,j) * vt3dh(k,i,j) > df * scr1(k,i,j)) then
               vt3dg(k,i,j) = vt3dc(k,i,j) * scr1(k,i,j)
            endif
         elseif (vt3dc(k,i,j) < 0.) then
            if (-vt3dg(k,i,j) * vt3dh(k+1,i,j) > df * scr1(k+1,i,j)) then
               vt3dg(k,i,j) = vt3dc(k,i,j) * scr1(k+1,i,j)
            endif
         endif
      enddo
   enddo
enddo

! Compute flux divergence

do j = ja,jz
   do i = ia,iz
      lpw = nint(flpw(i,j))
      do k = lpw,m1-1
          scr1(k,i,j) = scr1(k,i,j)  &
            + vt3dh(k,i,j) * (vt3dg(k-1,i,j) - vt3dg(k,i,j)  &
            + scp(k,i,j) * (vt3dc(k,i,j) - vt3dc(k-1,i,j)))
    enddo
   enddo
enddo

return
end

!     ****************************************************************

subroutine advtndc_adap(m1,m2,m3,ia,iz,ja,jz,flpw,scp,sca,sct,dtl,mynum)

implicit none
integer :: m1,m2,m3,ia,iz,ja,jz,mynum,i,j,k
real, dimension(m2,m3) :: flpw
integer :: lpw
real :: dtl,dtli
real, dimension(m1,m2,m3) :: scp,sca,sct

dtli = 1. / dtl
do j = ja,jz
   do i = ia,iz
      lpw = nint(flpw(i,j))
      do k = lpw,m1-1
         sct(k,i,j) = sct(k,i,j) + (sca(k,i,j) - scp(k,i,j)) * dtli
      enddo
   enddo
enddo

return
end
