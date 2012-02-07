!==========================================================================================!
!==========================================================================================!
!     This subroutine initialises the vectors containing the grid spacing vectors (actual- !
! ly the inverse of the grid spacing), so we don't need to recalculate them every time we  !
! advect stuff.                                                                            !
!------------------------------------------------------------------------------------------!
subroutine init_grid_spacing(m1,m2,m3,dxt,dyt,dzt,fmapt,rtgt,dxtw,dytw,dztw)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)  :: m1
   integer                     , intent(in)  :: m2
   integer                     , intent(in)  :: m3
   real   , dimension(   m2,m3), intent(in)  :: dxt
   real   , dimension(   m2,m3), intent(in)  :: dyt
   real   , dimension(m1      ), intent(in)  :: dzt
   real   , dimension(   m2,m3), intent(in)  :: fmapt
   real   , dimension(   m2,m3), intent(in)  :: rtgt
   real   , dimension(m1,m2,m3), intent(out) :: dxtw
   real   , dimension(m1,m2,m3), intent(out) :: dytw
   real   , dimension(m1,m2,m3), intent(out) :: dztw
   !----- Local variables. ----------------------------------------------------------------!
   integer                                   :: i
   integer                                   :: j
   integer                                   :: k
   real                                      :: rtgti
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Loop over the horizontal grid points.  The Jacobian doesn't depend on Z, so dztw  !
   ! depends only on Z.                                                                    !
   !---------------------------------------------------------------------------------------!
   do j=1,m3
      do i=1,m2
         rtgti = 1.0 / rtgt(i,j)

         do k=1,m1
            dxtw(k,i,j) = 1.0 / (dxt(i,j) * fmapt(i,j) * rtgti)
            dytw(k,i,j) = 1.0 / (dyt(i,j) * fmapt(i,j) * rtgti)
            dztw(k,i,j) = 1.0 / dzt(k)
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine init_grid_spacing
!==========================================================================================!
!==========================================================================================!





!==========================================================================================!
!==========================================================================================!
!     BRAMS often uses the reference density, but in the new advection we must use the     !
! actual densities.                                                                        !
!------------------------------------------------------------------------------------------!
subroutine find_actual_densities(m1,m2,m3,rtp,rv,pp,pi0,theta,denst,densu,densv,densw)
   use therm_lib , only : virtt        & ! function
                        , exner2press  & ! function
                        , extheta2temp ! ! function
   use rconstants, only : rdry         ! ! intent(in)
   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer                     , intent(in)    :: m1
   integer                     , intent(in)    :: m2
   integer                     , intent(in)    :: m3
   real   , dimension(m1,m2,m3), intent(in)    :: rtp
   real   , dimension(m1,m2,m3), intent(in)    :: rv
   real   , dimension(m1,m2,m3), intent(in)    :: pp
   real   , dimension(m1,m2,m3), intent(in)    :: pi0
   real   , dimension(m1,m2,m3), intent(in)    :: theta
   real   , dimension(m1,m2,m3), intent(inout) :: denst
   real   , dimension(m1,m2,m3), intent(inout) :: densu
   real   , dimension(m1,m2,m3), intent(inout) :: densv
   real   , dimension(m1,m2,m3), intent(inout) :: densw
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: i
   integer                                     :: j
   integer                                     :: k
   real                                        :: exner
   real                                        :: pres
   real                                        :: temp
   real                                        :: tvir
   !---------------------------------------------------------------------------------------!


   !----- Check whether this run has water vapour or not. ---------------------------------!
   do j=1,m3
      do i=1,m2
         do k=1,m1
            !----- Find pressure and temperature. -----------------------------------------!
            exner        = pi0(k,i,j) + pp(k,i,j)
            pres         = exner2press(exner)
            temp         = extheta2temp(exner,theta(k,i,j))
            !----- Find the virtual temperature. ------------------------------------------!
            tvir         = virtt(temp,rv(k,i,j),rtp(k,i,j))
            !----- Density comes from gas law using virtual temperature. ------------------!
            denst(k,i,j) = pres / (rdry * tvir)
            !------------------------------------------------------------------------------!
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the densities at the u and v grid points.                                    !
   !---------------------------------------------------------------------------------------!
   call fill_dn0uv(m1,m2,m3,denst,densu,densv)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the densities at the w grid points.                                          !
   !---------------------------------------------------------------------------------------!
   do j=1,m3
      do i=1,m2
         do k=1,m1-1
            densw(k,i,j) = 0.5 * (denst(k,i,j) + denst(k+1,i,j))
         end do
         densw(m1,i,j) = densw(m1-1,i,j)
      end do
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine find_actual_densities
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!      This subroutine finds the average winds for the advection, and stores the values on !
! the "mid" variables.                                                                     !
!------------------------------------------------------------------------------------------!
subroutine find_avg_winds(m1,m2,m3,ia,iz,ja,jz,ka,kz,uc,up,vc,vp,wc,wp,fmapui,fmapvi       &
                         ,rtgt,rtgu,rtgv,f13t,f23t,uavg,vavg,wavg)
   use mem_grid, only : hw4  ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                     , intent(in)    :: m1
   integer                     , intent(in)    :: m2
   integer                     , intent(in)    :: m3
   integer                     , intent(in)    :: ia
   integer                     , intent(in)    :: iz
   integer                     , intent(in)    :: ja
   integer                     , intent(in)    :: jz
   integer                     , intent(in)    :: ka
   integer                     , intent(in)    :: kz
   real   , dimension(m1,m2,m3), intent(in)    :: uc
   real   , dimension(m1,m2,m3), intent(in)    :: up
   real   , dimension(m1,m2,m3), intent(in)    :: vc
   real   , dimension(m1,m2,m3), intent(in)    :: vp
   real   , dimension(m1,m2,m3), intent(in)    :: wc
   real   , dimension(m1,m2,m3), intent(in)    :: wp
   real   , dimension(   m2,m3), intent(in)    :: fmapui
   real   , dimension(   m2,m3), intent(in)    :: fmapvi
   real   , dimension(   m2,m3), intent(in)    :: rtgt
   real   , dimension(   m2,m3), intent(in)    :: rtgu
   real   , dimension(   m2,m3), intent(in)    :: rtgv
   real   , dimension(   m2,m3), intent(in)    :: f13t
   real   , dimension(   m2,m3), intent(in)    :: f23t
   real   , dimension(m1,m2,m3), intent(inout) :: uavg
   real   , dimension(m1,m2,m3), intent(inout) :: vavg
   real   , dimension(m1,m2,m3), intent(inout) :: wavg
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: i
   integer                                     :: j
   integer                                     :: k
   integer                                     :: im1
   integer                                     :: ip1
   integer                                     :: jm1
   integer                                     :: jp1
   integer                                     :: kp1
   real                                        :: rtgti
   real                                        :: fmrt_u
   real                                        :: fmrt_v
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Start by defining the ?mid terms as temporal averages using Cartesian             !
   ! coordinates.                                                                          !
   !---------------------------------------------------------------------------------------!
   do j=1,m3
      do i=1,m2
         do k=1,m1
            uavg(k,i,j) = 0.5 * (uc(k,i,j) + up(k,i,j))
            vavg(k,i,j) = 0.5 * (vc(k,i,j) + vp(k,i,j))
            wavg(k,i,j) = 0.5 * (wc(k,i,j) + wp(k,i,j))
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Terrain-following coordinate has slopes, and we must add the slope contribution   !
   ! to the vertical component so it becomes a true sigma-z velocity.                      !
   !---------------------------------------------------------------------------------------!
   do j=1,m3
      jm1 = max( 1,j-1)
      jp1 = min(m3,j+1)
      do i=1,m2
         im1 = max( 1,i-1)
         ip1 = min(m2,i+1)
         rtgti = 1.0 / rtgt(i,j)

         do k=1,kz
            kp1 = k+1
            wavg(k,i,j) = hw4(k)                                                           &
                        * ( f13t(i,j) * ( uavg(  k,  i,  j) + uavg(kp1,  i,  j)            &
                                        + uavg(  k,im1,  j) + uavg(kp1,im1,  j) )          &
                          + f23t(i,j) * ( vavg(  k,  i,  j) + vavg(kp1,  i,  j)            &
                                        + vavg(  k,  i,jm1) + vavg(kp1,  i,jm1) ) )        &
                        + wavg(k,i,j) * rtgti
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Add the map factors to the horizontal winds.                                      !
   !---------------------------------------------------------------------------------------!
   do j=1,m3
      do i=1,m2
         fmrt_u = fmapui(i,j) * rtgu(i,j)
         fmrt_v = fmapvi(i,j) * rtgv(i,j)

         do k=1,m1
            uavg(k,i,j) = uavg(k,i,j) * fmrt_u
            vavg(k,i,j) = vavg(k,i,j) * fmrt_v
         end do

      end do
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine find_avg_winds
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine finds the density terms as in Walcek.                                !
!------------------------------------------------------------------------------------------!
subroutine find_walcek_densities(dtime,m1,m2,m3,ia,iz,ja,jz,ka,kz,uavg,vavg,wavg           &
                                ,denst,densu,densv,densw,den0_wal,den1_wal,den2_wal        &
                                ,den3_wal,dxtw,dytw,dztw)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real                        , intent(in)    :: dtime
   integer                     , intent(in)    :: m1
   integer                     , intent(in)    :: m2
   integer                     , intent(in)    :: m3
   integer                     , intent(in)    :: ia
   integer                     , intent(in)    :: iz
   integer                     , intent(in)    :: ja
   integer                     , intent(in)    :: jz
   integer                     , intent(in)    :: ka
   integer                     , intent(in)    :: kz
   real   , dimension(m1,m2,m3), intent(in)    :: uavg
   real   , dimension(m1,m2,m3), intent(in)    :: vavg
   real   , dimension(m1,m2,m3), intent(in)    :: wavg
   real   , dimension(m1,m2,m3), intent(in)    :: denst
   real   , dimension(m1,m2,m3), intent(in)    :: densu
   real   , dimension(m1,m2,m3), intent(in)    :: densv
   real   , dimension(m1,m2,m3), intent(in)    :: densw
   real   , dimension(m1,m2,m3), intent(inout) :: den0_wal
   real   , dimension(m1,m2,m3), intent(inout) :: den1_wal
   real   , dimension(m1,m2,m3), intent(inout) :: den2_wal
   real   , dimension(m1,m2,m3), intent(inout) :: den3_wal
   real   , dimension(m1,m2,m3), intent(in)    :: dxtw
   real   , dimension(m1,m2,m3), intent(in)    :: dytw
   real   , dimension(m1,m2,m3), intent(in)    :: dztw
   !----- Local variables. ----------------------------------------------------------------!
   integer                                     :: i
   integer                                     :: j
   integer                                     :: k
   integer                                     :: im1
   integer                                     :: jm1
   integer                                     :: km1
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Find the four terms.                                                              !
   !---------------------------------------------------------------------------------------!
   do j=ja,jz
      jm1 = j-1
      do i=ia,iz
         im1 = i-1
         do k=ka,kz
            km1 = k-1

            den0_wal(k,i,j) = denst(k,i,j)
            den1_wal(k,i,j) = den0_wal(k,i,j) - dtime / dxtw(k,i,j)                        &
                                              * ( densu(  k,  i,  j) * uavg(  k,  i,  j)   &
                                                - densu(  k,im1,  j) * uavg(  k,im1,  j) )
            den2_wal(k,i,j) = den1_wal(k,i,j) - dtime / dytw(k,i,j)                        &
                                              * ( densv(  k,  i,  j) * vavg(  k,  i,  j)   &
                                                - densv(  k,  i,jm1) * vavg(  k,  i,jm1) )
            den3_wal(k,i,j) = den2_wal(k,i,j) - dtime / dztw(k,i,j)                        &
                                              * ( densw(  k,  i,  j) * wavg(  k,  i,  j)   &
                                                - densw(km1,  i,  j) * wavg(km1,  i,  j) )
         end do
      end do
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine find_walcek_densities
!==========================================================================================!
!==========================================================================================!
