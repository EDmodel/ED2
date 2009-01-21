!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine latbnd()

use mem_tend
use mem_basic
use mem_grid
use node_mod

implicit none

!     This routine drives the computation of the radiative lateral
!     boundary condition for normal velocity on the coarsest grid
!     and the recomputation of all boundary tendencies on nested grids
!     after the first nested grid timestep.

if (nxtnest(ngrid) .eq. 0) then

!         Radiative and/or mesoscale compensation region lateral
!            boundary conditions.

   if (ibnd .le. 3 .or. jbnd .le. 3) then

      call latnormv(mzp,mxp,myp,ia,iz,ja,jz,ibcon    &
         ,grid_g(ngrid)%flpu   ,grid_g(ngrid)%flpv   &
         ,basic_g(ngrid)%up    ,basic_g(ngrid)%uc    &
         ,tend%ut              ,basic_g(ngrid)%vp    &
         ,basic_g(ngrid)%vc    ,tend%vt              &
         ,grid_g(ngrid)%dxt    ,grid_g(ngrid)%dyt    )

   endif

endif
return
end

!     *****************************************************************

subroutine latnormv(m1,m2,m3,ia,iz,ja,jz,ibcon,flpu,flpv  &
   ,up,uc,ut,vp,vc,vt,dxt,dyt)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j,k
real, dimension(m2,m3) :: flpu,flpv
integer :: lpu,lpv

real :: dxl,dxr,cphx,cphy
real, dimension(m1,m2,m3) :: up,uc,ut,vp,vc,vt
real, dimension(m2,m3) :: dxt,dyt

!     This routine ultimately updates tendencies at lateral boundaries
!     after first diagnosing appropriate phase speeds.
!
!     IBND and JBND are flags for the radiative type in the X and Y
!     direction. Their meaning is:
!
!        IBND=1......Klemp-Wilhelmson (1978) type; phase speed given
!                    by CPHAS
!        IBND=2......Klemp-Lilly (1980) type; doppler shifted phase
!                    speed constant with height and diagnosed from
!                    average of Orlanski speeds, i.e. function of
!                    only (X,Y)
!        IBND=3......Orlanski(1974) type; Phase speeds diagnosed
!                    from local conditions and function of (x,y,z)

!     Calculate the diagnostic phase
!     speed. The Orlanski(1976) leapfrog method is to use three time
!     levels of information, namely the T-2, T-1, and T level to
!     evaluate the phase speed given by - du/dt / du/dx = u + C.
!     If this is to be an Orlanski or Klemp-Lilly type boundary in the x
!     direction then this following diagnostic procedure is necessary.

!     If this is the first call to a routine, initialize the phase
!       speed arrays if necessary.

if (ibcon.eq.0) return

!     first compute "X" boundaries.

if (iand(ibcon,1) .ne. 0) then
   do j = 1,m3
      dxl = 1. / (dtlv * dxt(2,j))
      lpu = nint(flpu(1,j))
      do k = lpu,m1

         cphx = min(0.,max(-dxl,(up(k,1,j)-cphas)))
         ut(k,1,j) = ut(k,1,j) - cphx * dxt(2,j)  &
            * (up(k,2,j) + ut(k,2,j) * dtlv - up(k,1,j))

      enddo
   enddo
endif

if (iand(ibcon,2) .ne. 0) then
   do j = 1,m3
      dxr = 1. / (dtlv * dxt(m2-1,j))
      lpu = nint(flpu(m2-1,j))
      do k = lpu,m1

         cphx = max(0.,min(dxr,(up(k,m2-1,j)+cphas)))
         ut(k,m2-1,j) = ut(k,m2-1,j) - cphx * dxt(m2-1,j)  &
            * (up(k,m2-1,j) - (up(k,m2-2,j) + ut(k,m2-2,j) * dtlv))

      enddo
   enddo
endif

!     South and north boundaries.

if (jdim .eq. 1) then

   if (iand(ibcon,4) .ne. 0) then
      do i = 1,m2
         dxl = 1. / (dtlv * dyt(i,2))
         lpv = nint(flpv(i,1))
         do k = lpv,m1
            cphy = min(0.,max(-dxl,(vp(k,i,1)-cphas)))
            vt(k,i,1) = vt(k,i,1) - cphy * dyt(i,2)  &
               * (vp(k,i,2) + vt(k,i,2) * dtlv - vp(k,i,1))
         enddo
      enddo
   endif

   if (iand(ibcon,8) .ne. 0) then
      do i = 1,m2
         dxr = 1. / (dtlv * dyt(i,m3-1))
         lpv = nint(flpv(i,m3-1))
         do k = lpv,m1
            cphy = max(0.,min(dxr,(vp(k,i,m3-1)+cphas)))
            vt(k,i,m3-1) = vt(k,i,m3-1) - cphy * dyt(i,m3-1)  &
               * (vp(k,i,m3-1) - (vp(k,i,m3-2) + vt(k,i,m3-2) * dtlv))
         enddo
      enddo
   endif

endif
return
end

!     *****************************************************************

subroutine cyclic_set(n1,n2,n3,var,vpnt)

use mem_grid
use mem_scratch

implicit none

integer :: n1,n2,n3,ncycle,nshift,i,j,k
real, dimension(n1,n2,n3) :: var
character(len=*) :: vpnt

!     This routine inputs a 3d variable VAR dimensioned by N1,N2,N3
!     and sets the boundaries to cyclic symmetry.  This version is
!     set up for second order advection only.

if (ibnd .eq. 4) then
   ncycle = n2 - 3
   nshift = 0
   if (vpnt .eq. 'U') nshift = 1
   do j = 1,n3
      do k = 1,n1
         var(k,1,j) = var(k,ncycle+1,j)
         var(k,n2-nshift,j) = var(k,n2-nshift-ncycle,j)
      enddo
   enddo

   if (vpnt .ne. 'U') then
      do j = 1,n3
         do k = 1,n1
            var(k,2,j) = 0.5 * (var(k,2,j) + var(k,n2-1,j))
            var(k,n2-1,j) = var(k,2,j)
         enddo
      enddo
   endif
endif

if (jbnd .eq. 4 .and. jdim .eq. 1) then
   ncycle = n3 - 3
   nshift = 0
   if (vpnt .eq. 'V') nshift = 1
   do i = 1,n2
      do k = 1,n1
         var(k,i,1) = var(k,i,ncycle+1)
         var(k,i,n3-nshift) = var(k,i,n3-nshift-ncycle)
      enddo
   enddo

   if (vpnt .ne. 'V') then
      do i = 1,n2
         do k = 1,n1
            var(k,i,2) = 0.5 * (var(k,i,2) + var(k,i,n3-1))
            var(k,i,n3-1) = var(k,i,2)
         enddo
      enddo
   endif

endif
return
end subroutine cyclic_set

!*****************************************************************************

subroutine vpsets()

use mem_basic
use mem_grid
use node_mod

implicit none

if (nxtnest(ngrid) .eq. 0) then
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'U'  &
      ,basic_g(ngrid)%up   ,basic_g(ngrid)%up     &
      ,basic_g(ngrid)%vp   ,grid_g(ngrid)%dxu     &
      ,grid_g(ngrid)%dxm   ,grid_g(ngrid)%dyv     &
      ,grid_g(ngrid)%dym   )
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'V'  &
      ,basic_g(ngrid)%vp   ,basic_g(ngrid)%up     &
      ,basic_g(ngrid)%vp   ,grid_g(ngrid)%dxu     &
      ,grid_g(ngrid)%dxm   ,grid_g(ngrid)%dyv     &
      ,grid_g(ngrid)%dym   )
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'W'  &
      ,basic_g(ngrid)%wp   ,basic_g(ngrid)%up     &
      ,basic_g(ngrid)%vp   ,grid_g(ngrid)%dxu     &
      ,grid_g(ngrid)%dxm   ,grid_g(ngrid)%dyv     &
      ,grid_g(ngrid)%dym   )
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'P'  &
      ,basic_g(ngrid)%pp   ,basic_g(ngrid)%up     &
      ,basic_g(ngrid)%vp   ,grid_g(ngrid)%dxu     &
      ,grid_g(ngrid)%dxm   ,grid_g(ngrid)%dyv     &
      ,grid_g(ngrid)%dym   )                      
endif

if (nsttop .eq. 1) then
   call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%up,basic_g(ngrid)%up,'U')
   call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%vp,basic_g(ngrid)%vp,'V')
endif

if (nstbot .eq. 1) then
   if (if_adap == 0) then
      call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%up,'U')
      call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%vp,'V')
   else
      call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,grid_g(ngrid)%flpu  &
                      ,basic_g(ngrid)%up,'U')
      call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,grid_g(ngrid)%flpv  &
                      ,basic_g(ngrid)%vp,'V')
   endif
end if

call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%pp,basic_g(ngrid)%pp,'P')

if (if_adap == 0) then
   call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%pp,'P')
else
   call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,grid_g(ngrid)%flpw  &
                   ,basic_g(ngrid)%pp,'P')
end if

call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%wp,'W')

call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%up,'U')
call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%vp,'V')

if (nxtnest(ngrid) .eq. 0) then
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'U'  &
      ,basic_g(ngrid)%uc   ,basic_g(ngrid)%up     &
      ,basic_g(ngrid)%vp   ,grid_g(ngrid)%dxu     &
      ,grid_g(ngrid)%dxm   ,grid_g(ngrid)%dyv     &
      ,grid_g(ngrid)%dym   )
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'V'  &
      ,basic_g(ngrid)%vc   ,basic_g(ngrid)%up     &
      ,basic_g(ngrid)%vp   ,grid_g(ngrid)%dxu     &
      ,grid_g(ngrid)%dxm   ,grid_g(ngrid)%dyv     &
      ,grid_g(ngrid)%dym   )
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon ,'W' &
      ,basic_g(ngrid)%wc   ,basic_g(ngrid)%up     &
      ,basic_g(ngrid)%vp   ,grid_g(ngrid)%dxu     &
      ,grid_g(ngrid)%dxm   ,grid_g(ngrid)%dyv     &
      ,grid_g(ngrid)%dym   )
   call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'P'  &
      ,basic_g(ngrid)%pc   ,basic_g(ngrid)%up     &
      ,basic_g(ngrid)%vp   ,grid_g(ngrid)%dxu     &
      ,grid_g(ngrid)%dxm   ,grid_g(ngrid)%dyv     &
      ,grid_g(ngrid)%dym   )
end if

if (nsttop .eq. 1) then
   call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%uc,basic_g(ngrid)%uc,'U')
   call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%vc,basic_g(ngrid)%vc,'V')
end if

if (nstbot .eq. 1) then
   if (if_adap == 0) then
      call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%uc,'U')
      call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%vc,'V')
   else
      call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,grid_g(ngrid)%flpu  &
                      ,basic_g(ngrid)%uc,'U')
      call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,grid_g(ngrid)%flpv  &
                      ,basic_g(ngrid)%vc,'V')
   end if
end if


call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%pc,basic_g(ngrid)%pc,'P')

if (if_adap == 0) then
   call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%pc,'P')
else
   call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon,grid_g(ngrid)%flpw  &
                   ,basic_g(ngrid)%pc,'P')
endif

call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%wc,'W')
call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%uc,'U')
call dumset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,basic_g(ngrid)%vc,'V')

return
end

!*****************************************************************************

subroutine trsets()

use var_tables
use mem_basic
use mem_grid
use mem_turb
use node_mod

!for catt:
!use burns, only: latset_tracer
use catt_start, only: catt

implicit none

integer :: n

!     Apply lateral, top, and bottom boundary conditions.

do n = 1,num_scalar(ngrid)

   if (nxtnest(ngrid) .eq. 0) then
      if (ngrid .eq. 1 .and. ipara .eq. 0)  &
         call cyclic_set(nzp,nxp,nyp,scalar_tab(n,ngrid)%var_p,'T')



         if( catt == 1) then
	  if(n .le. (num_scalar(ngrid) - NADDSC) ) then

            call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'TR'   &
            ,scalar_tab(n,ngrid)%var_p  ,basic_g(ngrid)%up   &
            ,basic_g(ngrid)%vp   ,grid_g(ngrid)%dxu          & 
            ,grid_g(ngrid)%dxm   ,grid_g(ngrid)%dyv          & 
            ,grid_g(ngrid)%dym   )
	  else
            call latset_tracer(mzp,mxp,myp,ia,iz,ja,jz,ibcon &
            ,scalar_tab(n,ngrid)%var_p  ,basic_g(ngrid)%up   &
            ,basic_g(ngrid)%vp          ,grid_g(ngrid)%dxu   &
            ,grid_g(ngrid)%dxm          ,grid_g(ngrid)%dyv   &
            ,grid_g(ngrid)%dym          )
	 endif
	else
            call latset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,'TR'   &
            ,scalar_tab(n,ngrid)%var_p  ,basic_g(ngrid)%up   &
            ,basic_g(ngrid)%vp          ,grid_g(ngrid)%dxu   &
            ,grid_g(ngrid)%dxm          ,grid_g(ngrid)%dyv   &
            ,grid_g(ngrid)%dym          )
        endif
       !srf


   endif
   if (nsttop .eq. 1)  &
      call topset(mzp,mxp,myp,ia,iz,ja,jz,ibcon &
                 ,scalar_tab(n,ngrid)%var_p  ,scalar_tab(n,ngrid)%var_p  ,'T')
   if (nstbot .eq. 1)  then
      if (if_adap == 0) then
         call botset(mzp,mxp,myp,ia,iz,ja,jz,ibcon,scalar_tab(n,ngrid)%var_p  ,'T')
      else
         call botset_adap(mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
            ,grid_g(ngrid)%flpw,scalar_tab(n,ngrid)%var_p  ,'T')
      endif
   endif
enddo

!       Make sure all positive definite quantities remain such.

call tkeinit(mzp,mxp,myp,ia,iz,ja,jz)

call negadj1(mzp,mxp,myp,ia,iz,ja,jz)

return
end

!******************************************************************************

subroutine latset(m1,m2,m3,ia,iz,ja,jz,ibcon,vnam,ap,uc,vc,dxu,dxm,dyv,dym)

use mem_grid
use mem_scratch

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j,k,lbw,lbe,lbs,lbn
real :: thresh,dtlx,c1,dxr,dyr
real, dimension(m1,m2,m3) :: ap,uc,vc
real, dimension(m2,m3) :: dxu,dxm,dyv,dym
character(len=*) :: vnam

if (iand(ibcon,1) .gt. 0) lbw = ia - 1
if (iand(ibcon,2) .gt. 0) lbe = iz + 1
if (iand(ibcon,4) .gt. 0) lbs = ja - 1
if (iand(ibcon,8) .gt. 0) lbn = jz + 1

thresh = 0.
if (vnam .eq. 'U' .or. vnam .eq. 'V' .or. vnam .eq. 'W' .or. vnam .eq.'P') then
   dtlx = dtlv
else
   dtlx = dtlt
endif

if (ibnd .ne. 4 .and. vnam .ne. 'U' .and. lsflg .ne. 3) then

!     Western and Eastern boundaries for zero gradient option

   if (lsflg .eq. 0) then
      if (iand(ibcon,1) .gt. 0) then
         do j = 1,m3
            do k = 1,m1
               ap(k,lbw,j) = ap(k,ia,j)
            enddo
         enddo
      endif
      if (iand(ibcon,2) .gt. 0) then
         do j = 1,m3
            do k = 1,m1
               ap(k,lbe,j) = ap(k,iz,j)
            enddo
         enddo
      endif
   else

!     Western boundary for lsflg = 1 or 2

      if (iand(ibcon,1) .gt. 0) then
         do j = 1,m3
            if (vnam .eq. 'V') then
               dxr = dxm(ia,j) / dxm(lbw,j)
               c1 = .5 * dtlx * dxm(lbw,j)
               do k = 1,m1
                  vctr17(k) = -c1 * (uc(k,lbw,j) + uc(k,lbw,j+jdim))
               enddo
            elseif (vnam .eq. 'W') then
               dxr = dxu(ia,j) / dxu(lbw,j)
               c1 = .5 * dtlx * dxu(lbw,j)
               do k = 1,m1
                  vctr17(k) = -c1 * (uc(k,lbw,j) + uc(k+1,lbw,j))
               enddo
            else
               dxr = dxu(ia,j) / dxu(lbw,j)
               c1 = dtlx * dxu(lbw,j)
               do k = 1,m1
                  vctr17(k) = -c1 * uc(k,lbw,j)
               enddo
            endif
            do k = 1,m1
               vctr18(k) = ap(k,ia,j) + dxr * (ap(k,ia,j) - ap(k,ia+1,j))
            enddo
            do k = 1,m1
               if (vctr17(k) .ge. thresh) then
                  ap(k,lbw,j) = vctr18(k)
               elseif (lsflg .eq. 1) then
                  ap(k,lbw,j) = ap(k,ia,j)
               endif
            enddo
         enddo
      endif

!     Eastern Boundary for LSFLG = 1 or 2

      if (iand(ibcon,2) .gt. 0) then
         do j = 1,m3
            if (vnam .eq. 'V') then
               dxr = dxm(iz-1,j) / dxm(iz,j)
               c1 = .5 * dtlx * dxm(iz,j)
               do k = 1,m1
                  vctr17(k) = c1 * (uc(k,iz,j) + uc(k,iz,j+jdim))
               enddo
            elseif (vnam .eq. 'W') then
               dxr = dxu(iz-1,j) / dxu(iz,j)
               c1 = .5 * dtlx * dxu(iz,j)
               do k = 1,m1
                  vctr17(k) = c1 * (uc(k,iz,j) + uc(k+1,iz,j))
               enddo
            else
               dxr = dxu(iz-1,j) / dxu(iz,j)
               c1 = dtlx * dxu(iz,j)
               do k = 1,m1
                  vctr17(k) = c1 * uc(k,iz,j)
               enddo
            endif
            do k = 1,m1
               vctr18(k) = ap(k,iz,j) + dxr * (ap(k,iz,j) - ap(k,iz-1,j))
            enddo
            do k = 1,m1
               if (vctr17(k) .ge. thresh) then
                  ap(k,lbe,j) = vctr18(k)
               elseif (lsflg .eq. 1) then
                  ap(k,lbe,j) = ap(k,iz,j)
               endif
            enddo
         enddo
      endif
   endif
endif

if(jdim.eq.1.and.jbnd.ne.4.and.vnam.ne.'V'.and.lsflg.ne.3)then

!     Southern and Northern boundaries for zero gradient option

  if (lsflg .eq. 0) then
     if (iand(ibcon,4) .gt. 0) then
        do i = 1,m2
           do k = 1,m1
              ap(k,i,lbs) = ap(k,i,ja)
           enddo
        enddo
     endif
     if (iand(ibcon,8) .gt. 0) then
        do i = 1,m2
           do k = 1,m1
              ap(k,i,lbn) = ap(k,i,jz)
           enddo
        enddo
     endif
  else

!     Southern boundary for LSFLG = 1 or 2

     if (iand(ibcon,4) .gt. 0) then
        do i = 1,m2
           if (vnam .eq. 'U') then
              dyr = dym(i,ja) / dym(i,lbs)
              c1 = .5 * dtlx * dym(i,lbs)
              do k = 1,m1
                 vctr17(k) = -c1 * (vc(k,i,lbs) + vc(k,i+1,lbs))
              enddo
           elseif (vnam .eq. 'W') then
              dyr = dyv(i,ja) / dyv(i,lbs)
              c1 = .5 * dtlx * dyv(i,lbs)
              do k = 1,m1
                 vctr17(k) = -c1 * (vc(k,i,lbs) + vc(k+1,i,lbs))
              enddo
           else
              dyr = dyv(i,ja) / dyv(i,lbs)
              c1 = dtlx * dyv(i,lbs)
              do k = 1,nz
                 vctr17(k) = -c1 * vc(k,i,lbs)
              enddo
           endif
           do k = 1,m1
              vctr18(k) = ap(k,i,ja) + dyr * (ap(k,i,ja) - ap(k,i,ja+1))
           enddo
           do k = 1,m1
              if (vctr17(k) .ge. thresh) then
                 ap(k,i,lbs) = vctr18(k)
              elseif (lsflg .eq. 1) then
                 ap(k,i,lbs) = ap(k,i,ja)
              endif
           enddo
        enddo
     endif

!     Northern Boundary for LSFLG = 1 or 2

     if (iand(ibcon,8) .gt. 0) then
        do i = 1,m2
           if (vnam .eq. 'U') then
              dyr = dym(i,jz-1) / dym(i,jz)
              c1 = .5 * dtlx * dym(i,jz)
              do k = 1,m1
                 vctr17(k) = c1 * (vc(k,i,jz) + vc(k,i+1,jz))
              enddo
           elseif (vnam .eq. 'W') then
              dyr = dyv(i,jz-1) / dyv(i,jz)
              c1 = .5 * dtlx * dyv(i,jz)
              do k = 1,m1
                 vctr17(k) = c1 * (vc(k,i,jz) + vc(k+1,i,jz))
              enddo
           else
              dyr = dyv(i,jz-1) / dyv(i,jz)
              c1 = dtlx * dyv(i,jz)
              do k = 1,m1
                 vctr17(k) = c1 * vc(k,i,jz)
              enddo
           endif
           do k = 1,m1
              vctr18(k) = ap(k,i,jz) + dyr * (ap(k,i,jz) - ap(k,i,jz-1))
           enddo
           do k = 1,m1
              if (vctr17(k) .ge. thresh) then
                 ap(k,i,lbn) = vctr18(k)
              elseif (lsflg .eq. 1) then
                 ap(k,i,lbn) = ap(k,i,jz)
              endif
           enddo
        enddo
     endif
  endif
endif

return
end subroutine latset

!     ******************************************************************

subroutine topset(m1,m2,m3,ia,iz,ja,jz,ibcon,ap,fa,vnam)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j

real :: dzmr,dztr
real, dimension(m1,m2,m3) :: ap,fa
character(len=*) :: vnam

dzmr = dzm(m1-2) / dzm(m1-1)
dztr = dzt(m1-2) / dzt(m1-1)

!     Computation of all prognostic variables (other than W) at
!       level NZP by extrapolation from below

if (vnam .eq. 'U' .or. vnam .eq. 'V' .or. vnam .eq. 'P') then
   do j = 1,m3
      do i = 1,m2
         ap(m1,i,j) = ap(m1-1,i,j) + dzmr * (ap(m1-1,i,j) - ap(m1-2,i,j))
      enddo
   enddo
endif
if (vnam .eq. 'T') then
   do j = 1,m3
      do i = 1,m2
         ap(m1,i,j) = max(0.,ap(m1-1,i,j)+dzmr*(ap(m1-1,i,j)-ap(m1-2,i,j)))
      enddo
   enddo
endif

return
end

!     ******************************************************************

subroutine botset(m1,m2,m3,ia,iz,ja,jz,ibcon,aa,vnam)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j
real :: dzmr
real, dimension(m1,m2,m3) :: aa
character(len=*) :: vnam

if (vnam .eq. 'P') then
   dzmr = dzm(2) / dzm(1)
   do i = 1,m2
      do j = 1,m3
         aa(1,i,j) = aa(2,i,j) + (aa(2,i,j) - aa(3,i,j)) * dzmr
      enddo
   enddo
else
   do i = 1,m2
      do j = 1,m3
         aa(1,i,j) = aa(2,i,j)
      enddo
   enddo
endif

return
end

!     ******************************************************************

subroutine dumset(m1,m2,m3,ia,iz,ja,jz,ibcon,aa,vnam)

use mem_grid

implicit none

integer :: m1,m2,m3,ia,iz,ja,jz,ibcon,i,j,k
real, dimension(m1,m2,m3) :: aa
character(len=*) :: vnam

if (vnam .eq. 'U' .and. iand(ibcon,2) .gt. 0) then
   do j = 1,m3
      do k = 1,m1
         aa(k,m2,j) = aa(k,m2-1,j)
      enddo
   enddo
elseif (vnam .eq. 'V' .and. iand(ibcon,8) .gt. 0) then
   do i = 1,m2
      do k = 1,m1
         aa(k,i,m3) = aa(k,i,m3-jdim)
      enddo
   enddo
elseif (vnam .eq. 'W') then
   do j = 1,m3
      do i = 1,m2
         aa(m1,i,j) = aa(m1-1,i,j)
      enddo
   enddo
endif

return
end

!     *****************************************************************

subroutine rayft()

use mem_tend
use mem_scratch
use mem_basic
use mem_grid
use node_mod
use therm_lib, only : virtt,vapour_on

implicit none

integer :: mxyzp,ii,ind,i,j,k

!     This routine is the rayleigh friction driver for the
!     theta friction and is called from the long timestep.

if (nfpt .eq. 0 .or. distim .eq. 0.) return

mxyzp = mxp * myp * mzp

!     First load past virtual theta into temporary.

if (vapour_on) then
   ind = 0
   do j = 1,mmyp(ngrid)
      do i = 1,mmxp(ngrid)
         do k = 1,mmzp(ngrid)
            ind = ind + 1
            scratch%vt3da(ind) = virtt(basic_g(ngrid)%theta(k,i,j)                         &
                                      ,basic_g(ngrid)%rv(k,i,j),basic_g(ngrid)%rtp(k,i,j))
         enddo
      enddo
   enddo
else
  call atob(mxyzp,basic_g(ngrid)%theta,scratch%vt3da)
endif

!     Now get rayleigh friction tendency

if (if_adap == 0) then

   call rayf(4,mzp,mxp,myp,ia,iz,ja,jz,ibcon      &
      ,scratch%vt3da        ,basic_g(ngrid)%th0   &
      ,tend%tht             ,grid_g(ngrid)%rtgt   &
      ,grid_g(ngrid)%topt                         )

else

   call rayf_adap(4,mzp,mxp,myp,ia,iz,ja,jz,ibcon  &
      ,grid_g(ngrid)%flpw          ,scratch%vt3da  &
      ,basic_g(ngrid)%th0         ,tend%tht        )

endif

return
end

!******************************************************************************

subroutine rayf(ifrom,m1,m2,m3,ia,iz,ja,jz,ibcon,var,th0,tht,rtgx,topx)

use mem_grid
use mem_scratch
use ref_sounding

implicit none

integer :: ifrom,m1,m2,m3,ia,iz,ja,jz,ibcon
real, dimension(m1,m2,m3) :: var,th0,tht
real, dimension(m2,m3) :: rtgx,topx

real :: zmkf,c1,c2
integer :: kf,i,j,k

!     This routine calculates rayleigh friction terms velocity and theta_il

if (nfpt .eq. 0 .or. distim .le. 0) return
kf = nnz(1) - nfpt
zmkf = zmn(kf,1)
c1 = 1. / (distim * (ztop - zmkf))
c2 = dts * c1
goto(100,200,300,400) ifrom
100   continue

!     u friction

do j = ja,jz
   do i = ia,iz
      do k = 1,nzp
         vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
      enddo
      call htint(nzp,u01dn(1,ngrid),zt,nzp,vctr5,vctr2)
      do k = nz,2,-1
         if (vctr2(k) .le. zmkf) go to 10

         var(k,i,j) = var(k,i,j) + c2 * (vctr2(k) - zmkf)  &
            * (vctr5(k) - var(k,i,j))

      enddo
10         continue
   enddo
enddo
return
200   continue

!     V friction

if (jdim .eq. 0 .and. icorflg .eq. 0) return
do j = ja,jz
   do i = ia,iz
      do k = 1,nzp
         vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
      enddo
      call htint(nzp,v01dn(1,ngrid),zt,nzp,vctr5,vctr2)
      do k = nz,2,-1
         if (vctr2(k) .le. zmkf) go to 20
         var(k,i,j) = var(k,i,j) + c2 * (vctr2(k) - zmkf)  &
            * (vctr5(k) - var(k,i,j))
      enddo
20         continue
   enddo
enddo
return
300   continue

!     W friction

do j = ja,jz
   do i = ia,iz
      do k = nz,2,-1
         vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
         if (vctr2(k) .le. zmkf) go to 30
         var(k,i,j) = var(k,i,j) - c2 * (vctr2(k) - zmkf) * var(k,i,j)
      enddo
30         continue
   enddo
enddo
return
400   continue

!     THETA FRICTION

do j = ja,jz
   do i = ia,iz
      do k = nz,2,-1
         vctr2(k) = zt(k) * rtgx(i,j) + topx(i,j)
         if (vctr2(k) .le. zmkf) go to 40
         tht(k,i,j) = tht(k,i,j) + c1 * (vctr2(k) - zmkf)  &
              * (th0(k,i,j) - var(k,i,j))
      enddo
40         continue
   enddo
enddo
return
end

