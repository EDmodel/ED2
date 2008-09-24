!############################# Change Log ##################################################
! 5.0.0
!
!###########################################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################################
!MLO - Removed from turb_k.f90. The original turb_k.f90 was taking too long
!      to compile, so I broke the original file in two. 

!------------------------------------------------------------------------------------------!

subroutine strain(m1,m2,m3,ia,iz,ja,jz,ia_1,ja_1,iz1,jz1,jd                                &
                 ,up,vp,wp,vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df,vt3dg,vt3dh,vt3di,vt3dn     &
                 ,scr2,idiffk)

  implicit none

  integer, intent(in)  :: m1   &
                        , m2   &
                        , m3   &
                        , ia   &
                        , iz   &
                        , ja   &
                        , jz   &
                        , ia_1 &
                        , ja_1 &
                        , iz1  &
                        , jz1  &
                        , jd   &
                        , idiffk

  real, dimension(m1,m2,m3), intent(in) :: up,vp,wp

  real, dimension(m1,m2,m3), intent(inout):: vt3da,      &
                                             vt3db,      &
                                             vt3dc,      &
                                             vt3dd,      &
                                             vt3de,      &
                                             vt3df,      &
                                             vt3dg,      &
                                             vt3dh,      &
                                             vt3di,      &
                                             vt3dn,      &
                                             scr2

  !local variables:
  integer              :: i,j,k


  call grad(m1, m2, m3, ia  , iz1, ja  , jz , up, vt3da, 'XDIR', 'UPNT')
  call grad(m1, m2, m3, ia_1, iz , ja_1, jz , vp, vt3db, 'XDIR', 'VPNT')
  call grad(m1, m2, m3, ia_1, iz , ja  , jz , wp, vt3df, 'XDIR', 'WPNT')

  call grad(m1, m2, m3, ia_1, iz , ja_1, jz , up, vt3dn, 'YDIR', 'UPNT')
  call grad(m1, m2, m3, ia  , iz , ja  , jz1, vp, vt3dc, 'YDIR', 'VPNT')
  call grad(m1, m2, m3, ia  , iz , ja_1, jz , wp, vt3dg, 'YDIR', 'WPNT')

  call grad(m1, m2, m3, ia_1, iz , ja  , jz , up, vt3dd, 'ZDIR', 'UPNT')
  call grad(m1, m2, m3, ia  , iz , ja_1, jz , vp, vt3de, 'ZDIR', 'VPNT')
  if(idiffk >= 3 .and. idiffk /= 7)then
     call grad(m1,m2,m3,ia,iz,ja,jz,wp,scr2,'ZDIR','WPNT')
  endif

  if (idiffk <= 2 .or. idiffk == 7) then
     do j = ja,jz
        do i = ia,iz
           do k = 2,m1-1
              vt3dh(k,i,j) =2. * (vt3da(k,i,j) * vt3da(k,i,j)  &
                   + vt3dc(k,i,j) * vt3dc(k,i,j))  &
                   + .0625 * (vt3db(k,i,j) + vt3db(k,i-1,j)  &
                   + vt3db(k,i,j-jd) + vt3db(k,i-1,j-jd)  &
                   + vt3dn(k,i,j) + vt3dn(k,i-1,j)  &
                   + vt3dn(k,i,j-jd) + vt3dn(k,i-1,j-jd)) ** 2
              vt3di(k,i,j) = .0625 * ((vt3dd(k,i,j) + vt3dd(k-1,i,j)  &
                   + vt3dd(k,i-1,j) + vt3dd(k-1,i-1,j)) ** 2  &
                   + (vt3de(k,i,j) + vt3de(k-1,i,j)  &
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
                   (vt3da(k,i,j) - vt3di(k,i,j)) ** 2  &
                   + (vt3dc(k,i,j) - vt3di(k,i,j)) ** 2  &
                   + ( scr2(k,i,j) - vt3di(k,i,j)) ** 2)  &
                   + .0625 * ((vt3db(k,i,j) + vt3db(k,i-1,j)  &
                   + vt3db(k,i,j-jd) + vt3db(k,i-1,j-jd)) ** 2  &
                   + (vt3dd(k,i,j) + vt3dd(k,i-1,j)  &
                   + vt3dd(k-1,i,j) + vt3dd(k-1,i-1,j)) ** 2  &
                   + (vt3de(k,i,j) + vt3de(k-1,i,j)  &
                   + vt3de(k,i,j-jd) + vt3de(k-1,i,j-jd)) ** 2)
              vt3di(k,i,j) = vt3dh(k,i,j)
           enddo
        enddo
     enddo
  endif

  return
end subroutine strain
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes the Brunt-Väisälä frequency squared (en2) based on the       !
! virtual temperature profile. It will consider the effect of condensed phase.             !
!------------------------------------------------------------------------------------------!
subroutine bruvais(m1,m2,m3,ia,iz,ja,jz,theta,rtp,rv,rtgt,flpw,en2)
   use mem_scratch, only :  & !
           vctr12           & ! intent(out) - Virtual temperature
          ,vctr1            & ! intent(out) - Height above ground
          ,vctr5            & ! intent(out) - Delta-z between k and k+1
          ,vctr6            & ! intent(out) - Delta-z between k-1 and k
          ,vctr10           & ! intent(out) - Ratio between z(k)-z(k-½) and z(k+½)-z(k-½)
          ,vctr25           & ! intent(out) - d(theta_v)/dz at k+½
          ,vctr26           ! ! intent(out) - d(theta_v)/dz at k-½

   use mem_grid, only    : &
           zt              & ! intent(in)
          ,zm              & ! intent(in)
          ,nzp             & ! intent(in)
          ,nz              & ! intent(in)
          ,nzpmax          ! ! intent(in)

   use rconstants, only  : &
           g               ! ! intent(in)

   use therm_lib, only   : &
           virtt           & ! function
          ,vapour_on       ! ! intent(in)

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer                  , intent(in   ) :: m1,m2,m3 ! Z,X,Y dimensions         [   ---]
   integer                  , intent(in   ) :: ia,iz    ! West-East node bound.    [   ---]
   integer                  , intent(in   ) :: ja,jz    ! South-North node bound.  [   ---]
   real, dimension(m1,m2,m3), intent(in   ) :: theta    ! Potential temperature    [     K]
   real, dimension(m1,m2,m3), intent(in   ) :: rtp      ! Total mixing ratio       [ kg/kg]
   real, dimension(m1,m2,m3), intent(in   ) :: rv       ! Vapour mixing ratio      [ kg/kg]       
   real, dimension(   m2,m3), intent(in   ) :: rtgt     ! Sigma-z correction       [   m/m]
   real, dimension(   m2,m3), intent(in   ) :: flpw     ! Lowest point in W grid   [   ---]
   real, dimension(m1,m2,m3), intent(inout) :: en2      ! (Brunt-Väisälä freq.)²   [   Hz²]
   !----- Local variables -----------------------------------------------------------------!
   integer :: i,j,k,ki,k2,k1
   !---------------------------------------------------------------------------------------!

   do j = ja,jz
      do i = ia,iz

         k2=nint(flpw(i,j))
         k1=k2-1

         !----- Calculate the virtual potential temperature profile -----------------------!
         if (vapour_on) then
            !----- Consider the vapour and possibly the condensation effect ---------------!
            do k = k1,m1
               vctr12(k) = virtt(theta(k,i,j),rv(k,i,j),rtp(k,i,j))
            end do
         else
            !----- No water substance, theta_v = theta
            vctr12(k1:m1) = theta(k1:m1,i,j)
         endif


         do k = k1,m1
            vctr1(k) = rtgt(i,j)*(zt(k)-zm(k1))
         end do
         
         !---------------------------------------------------------------------------------!
         !    This correction is necessary because delta-z is not constant in most cases.  !
         ! If it were constant, then vctr10 is 0.5 which is the original case.             !
         !---------------------------------------------------------------------------------!
         do k = k2,m1-1
            vctr5(k)  = 1. / (vctr1(k+1) - vctr1(k))        ! 1/Delta_z(k+½)
            vctr6(k)  = 1. / (vctr1(k)   - vctr1(k-1))      ! 1/Delta_z(k-½)
            !----- vctr10 is the ratio between [z(k)-z(k-½)] and [z(k+½)+z(k-½)] ----------!
            vctr10(k) = (vctr1(k)-vctr1(k-1)) / (vctr1(k+1)-vctr1(k-1))
         end do

         do k = k2,m1-1
            vctr25(k)  = vctr5(k) * (vctr12(k+1) - vctr12(k))
            vctr26(k)  = vctr6(k) * (vctr12(k) - vctr12(k-1))
            en2(k,i,j) = g * (vctr10(k)*vctr25(k) + (1.-vctr10(k))*vctr26(k)) / vctr12(k)
         end do

      end do
   end do

   !----- **(JP)** Removed from the main loop to allow vectorisation ----------------------!
   do j = ja,jz
      do i = ia,iz
         k1=nint(flpw(i,j))-1
         en2(k1,i,j)=en2(k1,i,j)
      end do
   end do
   do j = ja,jz
      do i = ia,iz
         en2(nzp,i,j)=en2(nz,i,j)
      end do
   end do

   return
end subroutine bruvais
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine mxdefm(m1,m2,m3,ia,iz,ja,jz,ibcon,jd  &
     ,vt3dh,vt3di,vt3dj,vt3dk,scr1,scr2,dn0,rtgt,dxt,dyt,flpw,mynum)

  !     +-------------------------------------------------------------+
  !     \   this routine calculates the mixing coefficients with a    \
  !     \     smagorinsky-type deformational based k with an optional \
  !     \     unstable brunt-vaisala enhancement and an optional      \
  !     \     richardson number modification.                         \
  !     +-------------------------------------------------------------+

  use mem_scratch , only   : vctr1,    &     !INTENT(INOUT)
                             vctr2           !INTENT(INOUT)

  use mem_grid, only       : ngrid,    &     !INTENT(in)
                             zm,       &     !INTENT(in)
                             zt              !INTENT(in)

  use mem_turb, only       : csx,      &     !INTENT(in)
                             idiffk,   &     !INTENT(in)
                             rmin,     &     !INTENT(out)
                             rmax,     &     !INTENT(out)
                             zkhkm,    &     !INTENT(in)
                             csz,      &     !INTENT(in)
                             akmin           !INTENT(in)

  use rconstants, only     : vonk

  implicit none

  integer, intent(in) :: m1   &     !INTENT(in)
                       , m2   &     !INTENT(in)
                       , m3   &     !INTENT(in)
                       , ia   &     !INTENT(in)
                       , iz   &     !INTENT(in)
                       , ja   &     !INTENT(in)
                       , jz         !INTENT(in)

  integer, intent(in) :: ibcon, jd, mynum   !- EHE -> nao sao usadas!!!

  real, dimension(m1,m2,m3), INTENT(INOUT) :: vt3dh    &
                                            , vt3dk    &
                                            , scr1     &
                                            , scr2

  real, dimension(m1,m2,m3), INTENT(IN)    :: dn0      &
                                            , vt3di    &
                                            , vt3dj

  real, dimension(m2,m3), INTENT(IN) :: rtgt,dxt

  real, dimension(m2,m3), INTENT(IN) :: dyt  !- EHE -> nao e' usada!!!

  real, dimension(m2,m3), INTENT(IN) :: flpw


  !local variables:
  integer :: i,j,k,irich,ienfl,lpw

  real :: csx2,sq300,enfl,rchmax,c1,c2,c3,c4,akm,ambda,vkz2

  irich = 1
  ienfl = 1

  csx2 = csx(ngrid) * csx(ngrid)
  sq300 = 90000.
  if (idiffk(ngrid) == 2 .or. idiffk(ngrid) == 3) then
     rmin = -100.
     rmax = 1. / zkhkm(ngrid)
     do j = ja,jz
        do i = ia,iz
           lpw = nint(flpw(i,j))
           do k = lpw,m1-1
              vt3dk(k,i,j) = max(min(vt3dj(k,i,j)  &
                   / max(vt3di(k,i,j),1.e-15),rmax),rmin)
           enddo
        enddo
     enddo
     enfl = float(ienfl)
     rchmax = 1.0 + 9.0 * float(irich)
     lpw = nint(flpw(i,j))
     do k = lpw,m1
        vctr1(k) = csz(ngrid) * (zm(k) - zm(k-1))
        vctr2(k) = vctr1(k) * vctr1(k)
     enddo
  endif

  if (idiffk(ngrid) == 1 .or. idiffk(ngrid) == 7) then
     do j = ja,jz
        do i = ia,iz
           c2 = 1.0 / (dxt(i,j) * dxt(i,j))
           c3 = csx2 * c2
           akm = akmin(ngrid) * 0.075 * c2 ** (0.666667)
           lpw = nint(flpw(i,j))
           do k = lpw,m1-1
              scr2(k,i,j) = dn0(k,i,j)  &
                   * max(akm,c3*sqrt(vt3dh(k,i,j)))
           enddo
        enddo
     enddo
  elseif (idiffk(ngrid) == 2) then
     do j = ja,jz
        do i = ia,iz
           c1 = rtgt(i,j) * rtgt(i,j)
           c2 = 1.0 / (dxt(i,j) * dxt(i,j))
           c3 = csx2 * c2
           akm = akmin(ngrid) * 0.075 * c2 ** (0.666667)
           c4 = vonk * vonk * c1
           lpw = nint(flpw(i,j))
           do k = lpw,m1-1
              ! old csz*dz len  scr1(k,i,j) = dn0(k,i,j) * c1 * vctr2(k)

              ! asymptotic vertical scale length from bjorn with modifications:
              ! c3 is (csx * dx)^2, c1*vctr2(k) is (csz * dz)^2, sq300 is the square
              ! of 300 meters (used as a limit for horizontal grid spacing influence
              ! on vertical scale length), ambda is (asymptotic_vertical_length_scale)^2,
              ! and vkz2 is (vonk * height_above_surface)^2.

              ambda = max(c1 * vctr2(k),min(sq300,c3))
              vkz2 = c4 * zt(k) * zt(k)
              scr1(k,i,j) = dn0(k,i,j) * vkz2 / (vkz2 / ambda + 1)  &

              * (sqrt(vt3di(k,i,j))  &
                   + enfl * sqrt(max(0.,-vt3dj(k,i,j))))*min(rchmax  &
                   ,sqrt(max(0.,(1.-zkhkm(ngrid)*vt3dk(k,i,j)))))

              scr2(k,i,j) = dn0(k,i,j)  &
                   * max(akm,c3*sqrt(vt3dh(k,i,j)))
              vt3dh(k,i,j) = scr1(k,i,j) * zkhkm(ngrid)

           enddo
        enddo
     enddo
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !     call friclyr(nzp,nxp,nyp,a(iscr1),a(iustarl),a(itstarl)
     !    +    ,a(iustarw),a(itstarw),a(ipctlnd),a(itheta),a(irtgt))
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  elseif (idiffk(ngrid) == 3) then
     do j = ja,jz
        do i = ia,iz
           c1 = rtgt(i,j) * rtgt(i,j)
           lpw = nint(flpw(i,j))
           do k = lpw,m1-1
              scr1(k,i,j) = dn0(k,i,j) * c1 * vctr2(k)  &
                   * (sqrt(vt3dh(k,i,j))  &
                   + enfl * sqrt(max(0.,-vt3dj(k,i,j))))*min(rchmax  &
                   ,sqrt(max(0.,(1.-zkhkm(ngrid)*vt3dk(k,i,j)))))
              scr2(k,i,j) = scr1(k,i,j)
              vt3dh(k,i,j) = scr1(k,i,j) * zkhkm(ngrid)
           enddo
        enddo
     enddo
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !     call friclyr(nzp,nxp,nyp,a(iscr1),a(iustarl),a(itstarl)
     !    +    ,a(iustarw),a(itstarw),a(ipctlnd),a(itheta),a(irtgt))
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  endif

  return
end subroutine mxdefm

!------------------------------------------------------------------------------------------!

subroutine klbnd(m1,m2,m3,ibcon,jd,akay,dn0,flpw)

  implicit none

  integer, INTENT(IN)       :: m1     &
                             , m2     &
                             , m3     &
                             , ibcon  &
                             , jd

  real, dimension(m1,m2,m3), INTENT(INOUT) :: akay

  real, dimension(m1,m2,m3), INTENT(IN)    :: dn0

  real, dimension(m2,m3), INTENT(IN)    :: flpw

  !local variables:
  integer ::i,j,k,k2,lpw

  !     boundary conditions on a mixing coefficient

  do j = 1,m3
     do i = 1,m2
        k2=nint(flpw(i,j))
        do k=1,k2-1
           akay(k,i,j) = akay(k2,i,j) * dn0(k,i,j) / dn0(k2,i,j)
        enddo
        akay(m1,i,j) = akay(m1-1,i,j) * dn0(m1,i,j)  &
             / dn0(m1-1,i,j)
     enddo
  enddo
  if (iand(ibcon,1) .ne. 0) then

     do j = 1,m3
        do k = 1,m1
           akay(k,1,j) = akay(k,2,j)
        enddo
     enddo
  endif

  if (iand(ibcon,2) .ne. 0) then
     do j = 1,m3
        do k = 1,m1
           akay(k,m2,j) = akay(k,m2-1,j)
        enddo
     enddo
  endif

  if (jd .eq. 1) then
     if (iand(ibcon,4) .ne. 0) then
        do i = 1,m2
           do k = 1,m1
              akay(k,i,1) = akay(k,i,2)
           enddo
        enddo
     endif

     if (iand(ibcon,8) .ne. 0) then
        do i = 1,m2
           do k = 1,m1
              akay(k,i,m3) = akay(k,i,m3-1)
           enddo
        enddo
     endif
  endif

  return
end subroutine klbnd

!------------------------------------------------------------------------------------------!

subroutine mxdefm_tracer(m1,m2,m3,ia,iz,ja,jz,ibcon,jd  &
     ,vt3dh,khtr,dn0,dxt,dyt,flpw,mynum)  !ALF

  !     +-------------------------------------------------------------+
  !     \   this routine calculates the mixing coefficients with a    \
  !     \     smagorinsky-type deformational based k                  \
  !     +-------------------------------------------------------------+
  !       khtr = coef. dif. horizontal para tracers
  !       frtr = fator de reducao do Akmin dos campos meteorologicos
  !
  use mem_grid, only:  &
       ngrid                !INTENT(IN)
  use mem_turb, only:  &
       csx,            &    !INTENT(IN)
       akmin                !INTENT(IN)

  implicit none
  ! Arguments:
  integer, intent(in)                      :: m1, m2, m3, ia, iz, ja, jz, ibcon, jd, mynum
  ! ALF
  real, dimension(m2,m3), intent(in)    :: flpw
  real, dimension(m2,m3), intent(in)       :: dxt, dyt    !dyt nao eh usada
  real, dimension(m1,m2,m3), intent(in)    :: vt3dh, dn0
  real, dimension(m1,m2,m3), intent(inout) ::  khtr

  ! local variables:
  integer :: i,j,k,lpw
  real :: csx2,c2,c3,akm,frtr

  !frtr = 0.1
  frtr = 0.05
  csx2 = csx(ngrid) * csx(ngrid)

  do j = min(1, ja), max(jz, m3)
     do i = min(1, ia), max(iz, m2)

        c2 = 1.0 / (dxt(i,j) * dxt(i,j))
        c3 = csx2 * c2
        akm = frtr * akmin(ngrid) * 0.075 * c2 ** (0.666667)

        lpw = nint(flpw(i,j))
        do k = lpw,m1-1

           khtr(k,i,j) = dn0(k,i,j)  &
                * max(akm,c3*sqrt(vt3dh(k,i,j)))

        enddo

     enddo
  enddo

  return
end subroutine mxdefm_tracer

!------------------------------------------------------------------------------------------!


