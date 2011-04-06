!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     Acoustic terms small time-step driver: this routine calls all the necessary routines !
! to march the model through the small timesteps when the adaptive coordinate ("shaved-    !
! -eta") is used.                                                                          !
!------------------------------------------------------------------------------------------!
subroutine acoust_adap(m1,m2,m3,flpu,flpv,flpw,scr1,scr2,vt3da,vt3db,vt3dc,vt3dd,vt3de     &
                      ,vt3df,vt3dg,vt3dh,vt2da,dn0,pi0,th0,up,vp,wp,pp,ut,vt,wt,pt,dxu,dyv &
                      ,fmapui,fmapvi,dxt,dyt,fmapt,aru,arv,arw,volt,volu,volv,volw)

   use mem_grid
   use mem_scratch
   use node_mod

   implicit none

   !------ Arguments ----------------------------------------------------------------------!
   integer                   :: m1
   integer                   :: m2
   integer                   :: m3
   real, dimension(m1,m2,m3) :: dn0
   real, dimension(m1,m2,m3) :: pi0
   real, dimension(m1,m2,m3) :: th0
   real, dimension(m1,m2,m3) :: up
   real, dimension(m1,m2,m3) :: vp
   real, dimension(m1,m2,m3) :: wp
   real, dimension(m1,m2,m3) :: pp
   real, dimension(m1,m2,m3) :: aru
   real, dimension(m1,m2,m3) :: arv
   real, dimension(m1,m2,m3) :: arw
   real, dimension(m1,m2,m3) :: volt
   real, dimension(m1,m2,m3) :: volu
   real, dimension(m1,m2,m3) :: volv
   real, dimension(m1,m2,m3) :: volw
   real, dimension(   m2,m3) :: dxu
   real, dimension(   m2,m3) :: dyv
   real, dimension(   m2,m3) :: fmapui
   real, dimension(   m2,m3) :: fmapvi
   real, dimension(   m2,m3) :: dxt
   real, dimension(   m2,m3) :: dyt
   real, dimension(   m2,m3) :: fmapt
   real, dimension(   m2,m3) :: flpu
   real, dimension(   m2,m3) :: flpv
   real, dimension(   m2,m3) :: flpw
   real, dimension(*)        :: scr1
   real, dimension(*)        :: scr2
   real, dimension(*)        :: vt3da
   real, dimension(*)        :: vt3db
   real, dimension(*)        :: vt3dc
   real, dimension(*)        :: vt3dd
   real, dimension(*)        :: vt3de
   real, dimension(*)        :: vt3df
   real, dimension(*)        :: vt3dg
   real, dimension(*)        :: vt3dh
   real, dimension(*)        :: vt2da
   real, dimension(*)        :: ut
   real, dimension(*)        :: vt
   real, dimension(*)        :: wt
   real, dimension(*)        :: pt
   !----- Local variables. ----------------------------------------------------------------!
   real                      :: t1
   real                      :: w1
   real                      :: a1da2
   integer                   :: iter
   !---------------------------------------------------------------------------------------!


   do iter = 1,nnacoust(ngrid)

      !----- Get coefficients for computations. -------------------------------------------!
      dts = 2. * dtlt / nnacoust(ngrid)

      if (iter == 1)  &
         call coefz_adap(mzp,mxp,myp,ia,iz,ja,jz,flpw,vt3dc,vt3dd,vt3de,dn0,pi0,th0,a1da2  &
                        ,vt3df,vt3dg,scr2,vctr1,vctr2,arw,volt,volw)

      if (iter /= 1) then
         call mpilbc_driver('getst',4)
      end if

      call prdctu_adap(mzp,mxp,myp,ia,izu,ja,jz,ibcon,flpu,up,ut,pp,vt3da,th0,vt3db,dxu    &
                      ,vt3dh,aru,volu,mynum)

      !------------------------------------------------------------------------------------!
      !     MLO.  This block has been commented out based on RAMS-6.0.  I honestly don't   !
      !           see why this could cause problems, but it seems to work better.          !
      !------------------------------------------------------------------------------------!
      ! if (iter /= nnacoust(ngrid)) then
      !    call mpilbc_driver('sendst',2)
      ! end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      if (nxtnest(ngrid) == 0 .and. ipara == 0)  then
         call cyclic_set (nzp,nxp,nyp,up,'U')
      end if
      !------------------------------------------------------------------------------------!

      call prdctv_adap(mzp,mxp,myp,ia,iz,ja,jzv,ibcon,flpv,vp,vt,pp,vt3da,th0,vt3db,dyv    &
                      ,vt3dh,arv,volv)

      !------------------------------------------------------------------------------------!
      !     MLO.  This block has been commented out based on RAMS-6.0.  I honestly don't   !
      !           see why this could cause problems, but it seems to work better.          !
      !------------------------------------------------------------------------------------!
      ! if (iter /= nnacoust(ngrid)) then
      !    call mpilbc_driver('sendst',3)
      ! else
      !    call mpilbc_driver('sendst',5)
      ! end if
      !------------------------------------------------------------------------------------!
      call mpilbc_driver('sendst',5)

      !------------------------------------------------------------------------------------!
      if (nxtnest(ngrid) == 0 .and. ipara == 0) then
         call cyclic_set (nzp,nxp,nyp,vp,'V')
      end if
      !------------------------------------------------------------------------------------!


      call prdctw1_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,flpw,wp,wt,pp,vt3dc,a1da2,vt3dh)


      !------------------------------------------------------------------------------------!
      !     MLO.  This block has been commented out based on RAMS-6.0.  I honestly don't   !
      !           see why this could cause problems, but it seems to work better.          !
      !------------------------------------------------------------------------------------!
      ! if (iter /= nnacoust(ngrid)) then
      !    call mpilbc_driver('getst',2)
      !    call mpilbc_driver('getst',3)
      ! else
      !    call mpilbc_driver('getst',5)
      ! end if
      !------------------------------------------------------------------------------------!
      call mpilbc_driver('getst',5)


      call prdctp1_adap(mzp,mxp,myp,ia,iz,ja,jz,jdim,flpw,pp,up,vp,pi0,dn0,th0,pt,vt3da    &
                       ,vt3db,vt2da,fmapui,fmapvi,dxt,dyt,fmapt,aru,arv,volt,mynum)

      call prdctw2_adap(mzp,mxp,myp,ia,iz,ja,jz,flpw,wp,pp,vt3dc,vt3dd,vt3de,vt3dg,scr1    &
                       ,scr2,vt2da)

      call prdctw3_adap(mzp,mxp,myp,ia,iz,ja,jz,flpw,wp,scr1,vt3df,vt3dg,vt3dc,vt3dd,pp)

      !------------------------------------------------------------------------------------!
      if (nxtnest(ngrid) == 0 .and. ipara == 0) then
         call cyclic_set (nzp,nxp,nyp,wp,'W')
      end if
      !------------------------------------------------------------------------------------!

      call prdctp2_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,flpw,pp,wp,vt3dd,vt3de,mynum)

      !------------------------------------------------------------------------------------!
      if (nxtnest(ngrid) == 0 .and. ipara == 0) then
         call cyclic_set (nzp,nxp,nyp,pp,'T')
      end if
      !------------------------------------------------------------------------------------!

      if (iter /= nnacoust(ngrid)) then
         call mpilbc_driver('sendst',4)
      else
         call mpilbc_driver('fullst',6)
      end if

   end do

   return
end subroutine acoust_adap
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine prdctu_adap(m1,m2,m3,ia,iz,ja,jz,ibcon,flpu  &
   ,up,ut,pp,vt3da,th0,dpdx,dxu,vt3dh,aru,volu,mynum)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k,ibcon,mynum
real, dimension(m2,m3) :: flpu

real, dimension(m1,m2,m3) :: up,ut,pp,vt3da,th0,dpdx,vt3dh,aru,volu
real, dimension(m2,m3) ::  dxu
integer :: lpu

!     U prediction

call azero(m1*m2*m3,dpdx)

!     Calculate acoustic tendency (horizontal pressure gradient)

do j = ja,jz
   do i = ia,iz
      lpu = nint(flpu(i,j))
      do k = lpu,m1-1
         dpdx(k,i,j) = -(th0(k,i,j) + th0(k,i+1,j)) * .5  &
            * aru(k,i,j) / volu(k,i,j) * (pp(k,i+1,j) - pp(k,i,j))
      enddo
   enddo
enddo

if (distim .ne. 0.) then
   call rayf_adap(1,m1,m2,m3,ia,iz,ja,jz,ibcon,flpu(1,1),up,th0,vt3dh)
endif

do j = 1,m3
   do i = 1,m2
      lpu = nint(flpu(i,j))
      do k = lpu,m1
         up(k,i,j) = up(k,i,j) + dts * (dpdx(k,i,j) + ut(k,i,j))
      enddo
   enddo
enddo

if (nstbot == 1 .and. itopo == 1)  &
     call botset_adap(m1,m2,m3,ia,iz,ja,jz,ibcon,flpu,up,'U')

return
end

!    ******************************************************************

subroutine prdctv_adap(m1,m2,m3,ia,iz,ja,jz,ibcon,flpv  &
   ,vp,vt,pp,vt3da,th0,dpdy,dyv,vt3dh,arv,volv)
   
use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j,k
integer :: lpv
real, dimension(m2,m3) :: flpv

real, dimension (m1,m2,m3) :: vp,vt,th0,dpdy,pp,vt3da,vt3dh,arv,volv
real, dimension(m2,m3) :: dyv

!     V prediction

call azero(m1*m2*m3,dpdy)

if (jdim .eq. 1) then

!       calculate acoustic tendency (horizontal pressure gradient)

   do j = ja,jz
      do i = ia,iz
         lpv = nint(flpv(i,j))
         do k =lpv,m1-1
            dpdy(k,i,j) = -(th0(k,i,j) + th0(k,i,j+1)) * .5  &
               * arv(k,i,j) / volv(k,i,j) * (pp(k,i,j+1) - pp(k,i,j))
         enddo
      enddo
   enddo

endif

if (distim .ne. 0.) then
   call rayf_adap(1,m1,m2,m3,ia,iz,ja,jz,ibcon,flpv(1,1),vp,th0,vt3dh)
endif

do j = 1,m3
   do i = 1,m2
      lpv=nint(flpv(i,j))
      do k = lpv,m1
         vp(k,i,j) = vp(k,i,j) + dts * (dpdy(k,i,j) + vt(k,i,j))
      enddo
   enddo
enddo

if (nstbot == 1 .and. itopo == 1)  &
     call botset_adap(m1,m2,m3,ia,iz,ja,jz,ibcon,flpv,vp,'V')

return
end

!----------------------------------------------------------------------

subroutine prdctw1_adap(m1,m2,m3,ia,iz,ja,jz,ibcon,flpw  &
   ,wp,wt,pp,acoc,a1da2,vt3dh)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j,k
real, dimension(m2,m3) :: flpw
integer :: lpw
real :: a1da2
real, dimension(m1,m2,m3) :: wp,wt,pp,acoc,vt3dh

!     First part of prediction at I,J point

!     Compute forward part of Crank-Nickelson scheme. This will be total
!     W prediction for explicit case.

if (distim .ne. 0.) then
   call rayf_adap(1,m1,m2,m3,ia,iz,ja,jz,ibcon,flpw(1,1),wp,vt3dh,vt3dh)
endif

do j = 1,m3
   do i = 1,m2
      lpw = nint(flpw(i,j))
      do k = lpw,m1-2
         wp(k,i,j) = wp(k,i,j) + dts * wt(k,i,j)
      enddo
   enddo
enddo

do j = ja,jz
   do i = ia,iz
      lpw = nint(flpw(i,j))
      do k = lpw,m1-2
         wp(k,i,j) = wp(k,i,j) + a1da2 * acoc(k,i,j)  &
            * (pp(k,i,j) - pp(k+1,i,j))
      enddo
   enddo
enddo

return
end

!---------------------------------------------------------------------

subroutine  prdctw2_adap(m1,m2,m3,ia,iz,ja,jz,flpw  &
   ,wp,pp,acoc,acof,acog,amof,amog,acoaa,heatfx1)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real, dimension(m2,m3) :: flpw
integer :: ka,lpw

real, dimension(m1,m2,m3) :: wp,pp,acoc,acof,acog,amof,amog,acoaa
real, dimension(m2,m3) :: heatfx1

if (nsttop == 1) then
   do j = ja,jz
      do i = ia,iz
         wp(m1-1,i,j) = 0.
      enddo
   enddo
endif

if (impl == 1) then

!         First implicit part of the w prediction

   do j = ja,jz
      do i = ia,iz
         lpw = nint(flpw(i,j))
         do k = lpw,m1-2
            wp(k,i,j) = wp(k,i,j) - (pp(k+1,i,j) - pp(k,i,j)) * acoc(k,i,j)
         enddo
      enddo
   enddo

   do j = ja,jz
      do i = ia,iz
         ka = nint(flpw(i,j))
         amog(ka-1,i,j) = -wp(ka-1,i,j) / amof(ka-1,i,j)
         do k = ka,m1-2
            amog(k,i,j) = (-wp(k,i,j) - acoaa(k,i,j) * amog(k-1,i,j))  &
               / amof(k,i,j)
         enddo
      enddo
   enddo

endif
return
end

!---------------------------------------------------------------------

subroutine prdctw3_adap(m1,m2,m3,ia,iz,ja,jz,flpw  &
   ,wp,amog,amoe,amof,acoc,acof,pp)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k
real, dimension(m2,m3) :: flpw
integer :: lpw
real, dimension(m1,m2,m3) :: wp,pp,acoc,acof,amof,amog,amoe

!     Conclusion of implicit w prediction

if (impl == 1) then
   do j = ja,jz
      do i = ia,iz
         lpw = nint(flpw(i,j))
         do k = m1-2,lpw,-1
            wp(k,i,j) = amog(k,i,j) - amoe(k,i,j) * wp(k+1,i,j)
         enddo
      enddo
   enddo
endif

if (nstbot == 1) then
   do j = ja,jz
      do i = ia,iz
         lpw = nint(flpw(i,j))
         do k = 1,lpw-1
        !    wp(k,i,j) = wp(lpw(i,j),i,j)
            wp(k,i,j) = 0.
         enddo
      enddo
   enddo
endif

return
end

!--------------------------------------------------------------------

subroutine prdctp1_adap(m1,m2,m3,ia,iz,ja,jz,jd,flpw  &
   ,pp,up,vp,pi0,dn0,th0,pt,hdv,hfx  &
   ,hfx1,fmapui,fmapvi,dxt,dyt,fmapt,aru,arv,volt,mynum)


use mem_grid
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,jd,i,j,k,mynum
real, dimension(m2,m3) :: flpw
integer :: lpw
real :: rocvpct
real, dimension(m1,m2,m3) :: pp,up,vp,pi0,hdv,pt,hfx,dn0,th0,aru,arv,volt
real, dimension(m2,m3) :: hfx1,fmapui,fmapvi,dxt,dyt,fmapt

call azero(m1*m2*m3,hdv)
rocvpct = .5 * rocv *sspct ** 2

do j = 1,m3
   do i = 1,m2
      lpw = nint(flpw(i,j))
      do k = lpw,m1
         hfx(k,i,j) = dn0(k,i,j) * th0(k,i,j)
      enddo
   enddo
enddo

do j = ja,jz
   do i = ia,iz
      lpw = nint(flpw(i,j))
      do k = lpw,m1-1

         hdv(k,i,j) = -rocvpct * pi0(k,i,j) / (hfx(k,i,j) * volt(k,i,j))   &

            * ((up(k,i,j) * aru(k,i,j)   * (hfx(k,i,j) + hfx(k,i+1,j))     &
            - up(k,i-1,j) * aru(k,i-1,j) * (hfx(k,i,j) + hfx(k,i-1,j)))    &

            + (vp(k,i,j) * arv(k,i,j)    * (hfx(k,i,j) + hfx(k,i,j+jd))    &
            - vp(k,i,j-jd) * arv(k,i,j-jd) * (hfx(k,i,j) + hfx(k,i,j-jd))) )

      enddo
   enddo
enddo

do j = ja,jz
   do i = ia,iz
      lpw = nint(flpw(i,j))
      do k = lpw,m1
         pp(k,i,j) = pp(k,i,j) + (pt(k,i,j) + hdv(k,i,j)) * dts
      enddo
   enddo
enddo

return
end

!--------------------------------------------------------------------

subroutine prdctp2_adap(m1,m2,m3,ia,iz,ja,jz,ibcon,flpw,pp,wp,acof,acog,mynum)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j,k,mynum
integer :: lpw
real, dimension(m1,m2,m3) :: pp,wp,acof,acog
real, dimension(m2,m3) :: flpw

!           Finish pressure prediction

do j = ja,jz
   do i = ia,iz
      lpw = nint(flpw(i,j))
      do k = lpw,m1-1
         pp(k,i,j) = pp(k,i,j)  &
            + (wp(k,i,j) * acof(k,i,j) + wp(k-1,i,j) * acog(k,i,j))
      enddo
   enddo
enddo

if (nstbot .eq. 1) call botset_adap(m1,m2,m3,ia,iz,ja,jz,ibcon,flpw,pp,'P')

return
end

!******************************************************************************

subroutine coefz_adap(m1,m2,m3,ia,iz,ja,jz,flpw  &
   ,acoc,acof,acog,dn0,pi0,th0,a1da2,amoe,amof,acoaa,acobb,acocc  &
   ,arw,volt,volw)

use mem_grid
use mem_scratch
use rconstants

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,i,j,k

real :: dt2al2,a1da2,rdto2cv,dt2al2r,rdtr
real, dimension(m1,m2,m3) :: th0,pi0,dn0,acoc,acof,acog,amoe,amof,acoaa  &
                            ,arw,volt,volw
real, dimension(m2,m3) :: flpw
real, dimension(*) :: acobb,acocc
integer :: ka

! +--------------------------------------------------------------------+
! \   Calculate coefficients for the vertical pressure gradient        \
! \     and divergence terms.  These will be combined later for the    \
! \     implicit computations.                                         \
! +--------------------------------------------------------------------+

if (impl .eq. 1) then
   dt2al2 = dts * .75
   a1da2 = 1. / 3.
else
   dt2al2 = dts
   a1da2 = 1.
endif
rdto2cv = sspct ** 2 * rdry * dts / (2.0 * cv)

do j = ja,jz
   do i = ia,iz
      ka = nint(flpw(i,j))

!         Coefficient for the vertical pressure gradient term

      dt2al2r = .5 * dt2al2
      do k = ka-1,m1-1
         acoc(k,i,j) = dt2al2r * arw(k,i,j) / volw(k,i,j)  &
            * (th0(k,i,j) + th0(k+1,i,j))
      enddo

!         Coefficients for the vertical divergence term

      rdtr = rdto2cv
      do k = ka,m1
         vctr12(k) = dn0(k,i,j) * th0(k,i,j)
         vctr11(k) = rdtr * pi0(k,i,j) / (vctr12(k) * volt(k,i,j))
      enddo
      vctr12(ka-1) = dn0(ka-1,i,j) * th0(ka-1,i,j)
      do k = ka,m1-1
         acof(k,i,j) = -vctr11(k) * (vctr12(k) + vctr12(k+1)) * arw(k,i,j)
         acog(k,i,j) = vctr11(k) * (vctr12(k) + vctr12(k-1)) * arw(k-1,i,j)
      enddo
      acog(m1,i,j) = vctr11(m1) * (vctr12(m1) + vctr12(m1-1)) * arw(m1-1,i,j)

      do k = ka,m1-1
         acoaa(k,i,j) = acoc(k,i,j) * acog(k,i,j)
         acobb(k) = acoc(k,i,j) * (acof(k,i,j) - acog(k+1,i,j)) - 1.
         acocc(k) = -acoc(k,i,j) * acof(k+1,i,j)
      enddo
      acobb(ka-1) = -1.
      acocc(ka-1) = 0.
      acoaa(m1,i,j) = 0.
      acobb(m1) = -1.

      amof(ka-1,i,j) = acobb(ka-1)
      amoe(ka-1,i,j) = acocc(ka-1) / amof(ka-1,i,j)
      do k = ka,m1
         amof(k,i,j) = acobb(k) - acoaa(k,i,j) * amoe(k-1,i,j)
         amoe(k,i,j) = acocc(k) / amof(k,i,j)
      enddo

   enddo
enddo
return
end


