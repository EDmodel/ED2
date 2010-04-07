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
  use mem_turb, only : ibotflx
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


  call grad(m1, m2, m3, ia  , iz1, ja  , jz , up, vt3da, 'XDIR', 'UPNT',ibotflx)
  call grad(m1, m2, m3, ia_1, iz , ja_1, jz , vp, vt3db, 'XDIR', 'VPNT',ibotflx)
  call grad(m1, m2, m3, ia_1, iz , ja  , jz , wp, vt3df, 'XDIR', 'WPNT',ibotflx)

  call grad(m1, m2, m3, ia_1, iz , ja_1, jz , up, vt3dn, 'YDIR', 'UPNT',ibotflx)
  call grad(m1, m2, m3, ia  , iz , ja  , jz1, vp, vt3dc, 'YDIR', 'VPNT',ibotflx)
  call grad(m1, m2, m3, ia  , iz , ja_1, jz , wp, vt3dg, 'YDIR', 'WPNT',ibotflx)

  call grad(m1, m2, m3, ia_1, iz , ja  , jz , up, vt3dd, 'ZDIR', 'UPNT',ibotflx)
  call grad(m1, m2, m3, ia  , iz , ja_1, jz , vp, vt3de, 'ZDIR', 'VPNT',ibotflx)
  select case (idiffk)
  case (3:6)
     call grad(m1,m2,m3,ia,iz,ja,jz,wp,scr2,'ZDIR','WPNT',ibotflx)
  end select

  select case (idiffk)
  case (1,2,7,8)
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
  case default
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
  end select

  return
end subroutine strain
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes the Brunt-Väisälä frequency squared (en2) based on the       !
! virtual temperature profile. It will consider the effect of condensed phase.             !
!------------------------------------------------------------------------------------------!
subroutine bruvais(ibruvais,m1,m2,m3,ia,iz,ja,jz,pi0,pp,theta,rtp,rv,rtgt,flpw,en2)
   use mem_scratch, only :  & !
           vctr11           & ! intent(out) - Potential temperature
          ,vctr12           & ! intent(out) - Virtual potential temperature
          ,vctr1            & ! intent(out) - Height above ground
          ,vctr2            & ! intent(out) - coefficient #1 (either cl1 or ci1)
          ,vctr3            & ! intent(out) - coefficient #2 (either cl2 or ci2)
          ,vctr4            & ! intent(out) - coefficient #3 (either cl3 or ci3)
          ,vctr5            & ! intent(out) - Delta-z between k and k+1
          ,vctr6            & ! intent(out) - Delta-z between k-1 and k
          ,vctr10           & ! intent(out) - Ratio between z(k)-z(k-½) and z(k+½)-z(k-½)
          ,vctr19           & ! intent(out) - g / Height above ground
          ,vctr25           & ! intent(out) - d(theta_v)/dz at k+½
          ,vctr26           & ! intent(out) - d(theta_v)/dz at k-½
          ,vctr27           & ! intent(out) - Full Exner function                  [J/kg/K]
          ,vctr28           & ! intent(out) - Pressure                             [    Pa]
          ,vctr29           & ! intent(out) - Temperature                          [     K]
          ,vctr30           & ! intent(out) - Saturation mixing ratio              [ kg/kg]
          ,vctr31           ! ! intent(out) - Condensed  mixing ratio              [ kg/kg]

   use mem_grid, only    : &
           zt              & ! intent(in)
          ,zm              & ! intent(in)
          ,nzp             & ! intent(in)
          ,nz              & ! intent(in)
          ,nzpmax          ! ! intent(in)

   use rconstants, only  : &
           grav            & ! intent(in)
          ,t00             & ! intent(in)
          ,p00             & ! intent(in)
          ,alvl            & ! intent(in)
          ,alvi            & ! intent(in)
          ,rdry            & ! intent(in)
          ,cp              & ! intent(in)
          ,cpi             & ! intent(in)
          ,cpor            & ! intent(in)
          ,ep              ! ! intent(in)

   use therm_lib, only   : &
           virtt           & ! function
          ,rslf            & ! function
          ,rsif            & ! function
          ,vapour_on       & ! intent(in)
          ,cloud_on        & ! intent(in)
          ,bulk_on         ! ! intent(in)

   implicit none

   !----- Arguments -----------------------------------------------------------------------!
   integer                  , intent(in   ) :: ibruvais ! Method to compute N²     [   ---]
   integer                  , intent(in   ) :: m1,m2,m3 ! Z,X,Y dimensions         [   ---]
   integer                  , intent(in   ) :: ia,iz    ! West-East node bound.    [   ---]
   integer                  , intent(in   ) :: ja,jz    ! South-North node bound.  [   ---]
   real, dimension(m1,m2,m3), intent(in   ) :: pi0      ! Ref. Exner function      [J/kg/K]
   real, dimension(m1,m2,m3), intent(in   ) :: pp       ! Perturbation on Exner    [J/kg/K]
   real, dimension(m1,m2,m3), intent(in   ) :: theta    ! Potential temperature    [     K]
   real, dimension(m1,m2,m3), intent(in   ) :: rtp      ! Total mixing ratio       [ kg/kg]
   real, dimension(m1,m2,m3), intent(in   ) :: rv       ! Vapour mixing ratio      [ kg/kg]       
   real, dimension(   m2,m3), intent(in   ) :: rtgt     ! Sigma-z correction       [   m/m]
   real, dimension(   m2,m3), intent(in   ) :: flpw     ! Lowest point in W grid   [   ---]
   real, dimension(m1,m2,m3), intent(inout) :: en2      ! (Brunt-Väisälä freq.)²   [   Hz²]
   !----- Local variables -----------------------------------------------------------------!
   integer :: i,j,k,ki,k2,k1
   real :: temp,rvlsi,rvii
   !----- Local constants, for alternative method to compute N², test only ----------------!
   real, parameter :: cl1 = alvl / rdry
   real, parameter :: cl2 = ep * alvl ** 2 / (cp * rdry)
   real, parameter :: cl3 = alvl / cp
   real, parameter :: ci1 = alvi / rdry
   real, parameter :: ci2 = ep * alvi ** 2 / (cp * rdry)
   real, parameter :: ci3 = alvi / cp

   !---------------------------------------------------------------------------------------!

   jloop: do j = ja,jz
      iloop: do i = ia,iz

         k2=nint(flpw(i,j))
         k1=k2-1

         !----- Calculate the virtual potential temperature profile -----------------------!
         if (vapour_on) then
            !----- Consider the vapour and possibly the condensation effect ---------------!
            do k = k1,m1
               vctr11(k) = theta(k,i,j)
               vctr12(k) = virtt(vctr11(k),rv(k,i,j),rtp(k,i,j))
            end do
         else
            !----- No water substance, theta_v = theta ------------------------------------!
            do k = k1,m1
               vctr11(k) = theta(k,i,j)
               vctr12(k) = vctr11(k)
            end do
         end if

         !----- vctr1(k) is the height above ground ---------------------------------------!
         do k = k1,m1
            vctr1(k) = rtgt(i,j)*(zt(k)-zm(k1))
         end do
         
         
         !---------------------------------------------------------------------------------!
         !    This correction is necessary because delta-z is not constant in most         !
         ! cases. If it were constant, then vctr10 is 0.5 which is the original case.      !
         !---------------------------------------------------------------------------------!
         do k = k2,m1-1
            vctr5(k)  = 1. / (vctr1(k+1) - vctr1(k))        ! 1/Delta_z(k+½)
            vctr6(k)  = 1. / (vctr1(k)   - vctr1(k-1))      ! 1/Delta_z(k-½)
            !----- vctr10 is the ratio between [z(k)-z(k-½)] and [z(k+½)+z(k-½)] ----------!
            vctr10(k) = (vctr1(k)-vctr1(k-1)) / (vctr1(k+1)-vctr1(k-1))
         end do
         
         
         !---- Alternative way to compute, lacks citation ---------------------------------!
         if (ibruvais == 3) then
            do k=k2,m1-1
               vctr19(k) = grav / ((zt(k+1)-zt(k-1)) * rtgt(i,j))
            end do
         elseif (cloud_on .and. ibruvais == 2) then 

            do k = k2,m1-1
               vctr19(k) = grav / ((zt(k+1) - zt(k-1)) * rtgt(i,j))
            end do

            do k=k1,m1
               vctr27(k) = pi0(k,i,j) + pp(k,i,j)
               vctr28(k) = p00 * (cpi   * vctr27(k)) ** cpor
               vctr29(k) = theta(k,i,j) * vctr27(k)  * cpi

               !---------------------------------------------------------------------------!
               !    Deciding which coefficient to use. This is inconsistent with most      !
               ! thermodynamics, but I will use it anyway, just to check whether this is   !
               ! the cause of the dramatic difference.                                     !
               !---------------------------------------------------------------------------!
               if (bulk_on .and. vctr29(k) < t00-20.) then
                  vctr30(k) = rsif(vctr28(k),vctr29(k))
                  vctr2(k) = ci1
                  vctr3(k) = ci2
                  vctr4(k) = ci3
               else
                  vctr30(k) = rslf(vctr28(k),vctr29(k))
                  vctr2(k) = cl1
                  vctr3(k) = cl2
                  vctr4(k) = cl3
               end if
               vctr31(k) = max(rtp(k,i,j)/vctr30(k) -.999 ,0.)
            end do
            !----- Finding the lowest level above the -20C level and compute rsif ---------!
            if (bulk_on .and. any(vctr29(k1:m1) >= t00 -20.)) then
               ki = minloc(vctr29(k1:m1),dim=1,mask=vctr29(k1:m1) >= t00-20.) + k1-1
               rvlsi = rsif(vctr28(ki-1),vctr29(ki-1))
            else
               ki = k1
               rvlsi = rsif(vctr28(ki),vctr29(ki))
            end if
         end if

         !----- Computing the frequency according to the option. --------------------------!
         if (ibruvais == 3) then
            do k=k2,m1-1
               en2(k,i,j)=vctr19(k)*(vctr12(k+1)-vctr12(k-1))/vctr12(k)
            end do
         elseif (cloud_on .and. ibruvais == 2) then
            do k=k2,m1-1
               if (vctr31(k) > 0) then
                  if ( k == ki) then
                     rvii = rvlsi
                  else
                     rvii = vctr30(k-1)
                  end if
                  en2(k,i,j) = vctr19(k) * ((1. + vctr2(k) * vctr30(k) / vctr29(k))        &
                                          /(1. + vctr3(k) * vctr30(k) / vctr29(k) ** 2)    &
                             * ((vctr11(k+1) - vctr11(k-1)) / vctr11(k)                    &
                             + vctr4(k) / vctr29(k) * (vctr30(k+1) - rvii))                &
                             - (rtp(k+1,i,j) - rtp(k-1,i,j)))
               else
                  vctr25(k)  = vctr5(k) * (vctr12(k+1) - vctr12(k))
                  vctr26(k)  = vctr6(k) * (vctr12(k) - vctr12(k-1))
                  en2(k,i,j) = grav * ( vctr10(k)*vctr25(k) + (1. - vctr10(k))*vctr26(k))  &
                                      / vctr12(k)
               end if
            end do
         else
            do k = k2,m1-1
               vctr25(k)  = vctr5(k) * (vctr12(k+1) - vctr12(k))
               vctr26(k)  = vctr6(k) * (vctr12(k) - vctr12(k-1))
               en2(k,i,j) = grav * ( vctr10(k)*vctr25(k) + (1.-vctr10(k))*vctr26(k))       &
                                   / vctr12(k)
            end do
         end if

      end do iloop
   end do jloop

   !----- **(JP)** Removed from the main loop to allow vectorisation ----------------------!
   do j = ja,jz
      do i = ia,iz
         k1=nint(flpw(i,j))-1
         en2(k1,i,j)=en2(k2,i,j)
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
     ,vt3dh,vt3di,vt3dj,vt3dk,scr1,scr2,dn0,rtgt,dxt,dyt,flpw,akscal,mynum)

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
                             csz             !INTENT(in)

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
  real, dimension(m2,m3), INTENT(IN) :: akscal


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

  select case (idiffk(ngrid))
  case (1,7,8)
     do j = ja,jz
        do i = ia,iz
           c2 = 1.0 / (dxt(i,j) * dxt(i,j))
           c3 = csx2 * c2
           akm = akscal(i,j) * 0.075 * c2 ** (0.666667)
           lpw = nint(flpw(i,j))
           do k = lpw,m1-1
              scr2(k,i,j) = dn0(k,i,j)  &
                   * max(akm,c3*sqrt(vt3dh(k,i,j)))
           enddo
        enddo
     enddo
  case (2)
     do j = ja,jz
        do i = ia,iz
           c1 = rtgt(i,j) * rtgt(i,j)
           c2 = 1.0 / (dxt(i,j) * dxt(i,j))
           c3 = csx2 * c2
           akm = akscal(i,j) * 0.075 * c2 ** (0.666667)
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

  case (3)
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

  end select

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
     ,vt3dh,khtr,dn0,dxt,dyt,flpw,akscal,mynum)  !ALF

  !     +-------------------------------------------------------------+
  !     \   this routine calculates the mixing coefficients with a    \
  !     \     smagorinsky-type deformational based k                  \
  !     +-------------------------------------------------------------+
  !       khtr = coef. dif. horizontal para tracers
  !       frtr = fator de reducao do Akscal dos campos meteorologicos
  !
  use mem_grid, only:  &
       ngrid                !INTENT(IN)
  use mem_turb, only:  &
       csx                  !INTENT(IN)

  implicit none
  ! Arguments:
  integer, intent(in)                      :: m1, m2, m3, ia, iz, ja, jz, ibcon, jd, mynum
  ! ALF
  real, dimension(m2,m3), intent(in)    :: flpw
  real, dimension(m2,m3), intent(in)    :: akscal
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
        akm = frtr * akscal(i,j) * 0.075 * c2 ** (0.666667)

        lpw = nint(flpw(i,j))
        do k = lpw,m1-1

           khtr(k,i,j) = dn0(k,i,j)  &
                * max(akm,c3*sqrt(vt3dh(k,i,j)))

        enddo

     enddo
  enddo

  return
end subroutine mxdefm_tracer
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine will assign an appropriate value for akscal. This goes instead of      !
! former akmin in case the user wants different scales depending on topography.            !
!------------------------------------------------------------------------------------------!
subroutine akscal_init(n2,n3,ifm,topo,akscal)
   use mem_turb, only: akmin   & ! intent(in)
                     , akmax   & ! intent(in)
                     , hgtmin  & ! intent(in)
                     , hgtmax  ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer               , intent(in)  :: n2,n3,ifm
   real, dimension(n2,n3), intent(in)  :: topo
   real, dimension(n2,n3), intent(out) :: akscal
   !----- Local variables -----------------------------------------------------------------!
   integer                             :: i,j
   real                                :: hnorm
   !----- Constants -----------------------------------------------------------------------!
   real                  , parameter   :: scal = 3.
   !----- External functions --------------------------------------------------------------!
   real                  , external    :: errorfun
   !---------------------------------------------------------------------------------------!

   if (akmax(ifm) > 0.) then
      do j=1,n3
         do i=1,n2
            hnorm = scal * (2.*(topo(i,j)-hgtmin(ifm))/(hgtmax(ifm)-hgtmin(ifm)) - 1.)
            akscal(i,j) = akmin(ifm) + 0.5 * (akmax(ifm)-akmin(ifm))*(1.+errorfun(hnorm))
         end do
      end do
   else
      !----- User doesn't want different akscal, use akmin only. --------------------------!
      do j=1,n3
         do i=1,n2
            akscal(i,j) = akmin(ifm)
         end do
      end do
   end if

   return   
end subroutine akscal_init
!==========================================================================================!
!==========================================================================================!
