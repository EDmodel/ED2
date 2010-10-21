!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!
!    This subroutine computes the missing terms for the prognostic Exner function          !
! equation. Reference:                                                                     !
!                                                                                          !
! MEDVIGY, D. M.; MOORCROFT, P. R.; AVISSAR, R.; WALKO, R. L.; Mass conservation and       !
!   atmospheric dynamics in the Regional Atmospheric Modeling System (RAMS). Environ.      !
!   Fluid Mech., v. 5, p. 109-134, 2005.                                                   !
!                                                                                          !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine is the main driver for the full Exner prognostic equation. All terms  !
! except the heat flux are computed through here. Heat flux is always computed             !
!------------------------------------------------------------------------------------------!
subroutine exevolve(m1,m2,m3,ifm,ia,iz,ja,jz,izu,jzv,jdim,mynum,edt,key)


   use mem_basic,   only: basic_g
   use mem_grid,    only: grid_g, itopo
   use mem_mass,    only: mass_g
   use mem_tend,    only: tend
   use mem_scratch, only: scratch
   use therm_lib,   only: vapour_on

   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   character(len=*) , intent(in) :: key
   integer          , intent(in) :: m1,m2,m3,ifm,ia,iz,ja,jz,izu,jzv,jdim,mynum
   real             , intent(in) :: edt
   !----- Local variables -----------------------------------------------------------------!
   integer :: i,j,k
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     Copying the vapour and total mixing ratio to scratch arrays if this is not a      !
   ! "dry" run. Otherwise, leave the values equal to zero, which will impose theta=theta_v !
   !---------------------------------------------------------------------------------------!
   call azero2(m1*m2*m3,scratch%vt3dp,scratch%vt3dq)
   if (vapour_on) then
      !----- If water is allowed, copy them to scratch arrays -----------------------------!
      call atob(m1*m2*m3,basic_g(ifm)%rtp,scratch%vt3dp)
      call atob(m1*m2*m3,basic_g(ifm)%rv,scratch%vt3dq)
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Now we compute the contributions on pressure from different terms: Advection, com- !
   ! pression and heating, following equation (12) of Medvigy et al. (2005). Advection and !
   ! compression can be computed at the same time, since they don't depend on any other    !
   ! derivative. Heating, on the other hand, needs to be computed in two steps.            !
   !---------------------------------------------------------------------------------------!
   select case (trim(key))
   case ('ADV')
      !------------------------------------------------------------------------------------!
      !   Advection + Compression.                                                         !
      !------------------------------------------------------------------------------------!
      !----- Initialization ---------------------------------------------------------------!
      call azero(m1*m2*m3,mass_g(ifm)%thvlast)

      !----- Compute advective term -------------------------------------------------------!
      call exadvlf(m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,itopo                                 &
                  ,grid_g(ifm)%rtgu                       ,grid_g(ifm)%fmapui              &
                  ,grid_g(ifm)%rtgv                       ,grid_g(ifm)%fmapvi              &
                  ,grid_g(ifm)%f13t                       ,grid_g(ifm)%f23t                &
                  ,grid_g(ifm)%rtgt                       ,grid_g(ifm)%fmapt               &
                  ,grid_g(ifm)%dxt                        ,grid_g(ifm)%dyt                 &
                  ,basic_g(ifm)%uc                        ,basic_g(ifm)%dn0u               &
                  ,basic_g(ifm)%vc                        ,basic_g(ifm)%dn0v               &
                  ,basic_g(ifm)%dn0                       ,basic_g(ifm)%wc                 &
                  ,basic_g(ifm)%pc                        ,tend%pt                         )
      !----- Calculate compression term ---------------------------------------------------!
      call excondiv(m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,itopo                                &
                   ,basic_g(ifm)%uc                       ,basic_g(ifm)%vc                 &
                   ,basic_g(ifm)%wc                       ,basic_g(ifm)%pc                 &
                   ,tend%pt                               ,grid_g(ifm)%dxt                 &
                   ,grid_g(ifm)%dyt                       ,grid_g(ifm)%rtgt                &
                   ,grid_g(ifm)%rtgu                      ,grid_g(ifm)%rtgv                &
                   ,grid_g(ifm)%f13t                      ,grid_g(ifm)%f23t                &
                   ,grid_g(ifm)%fmapt                     ,grid_g(ifm)%fmapui              &
                   ,grid_g(ifm)%fmapvi                    )
      !----- Put theta_v from last timestep into memory ------------------------------------!
      call fill_thvlast(m1,m2,m3,ia,iz,ja,jz                                               &
                       ,mass_g(ifm)%thvlast                    ,basic_g(ifm)%theta         &
                       ,scratch%vt3dp                          ,scratch%vt3dq              )
      
   case ('THA')
      !------------------------------------------------------------------------------------!
      !   Advection part of the heating term.                                              !
      !------------------------------------------------------------------------------------!
      call advect_theta(m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,mynum,edt                        &
                       ,basic_g(ifm)%up                  ,basic_g(ifm)%uc                  &
                       ,basic_g(ifm)%vp                  ,basic_g(ifm)%vc                  &
                       ,basic_g(ifm)%wp                  ,basic_g(ifm)%wc                  &
                       ,basic_g(ifm)%pi0                 ,basic_g(ifm)%pc                  &
                       ,tend%pt                          ,basic_g(ifm)%theta               &
                       ,scratch%vt3dp                    ,scratch%vt3dq                    &
                       ,basic_g(ifm)%dn0                 ,basic_g(ifm)%dn0u                &
                       ,basic_g(ifm)%dn0v                ,grid_g(ifm)%rtgt                 &
                       ,grid_g(ifm)%rtgu                 ,grid_g(ifm)%rtgv                 &
                       ,grid_g(ifm)%fmapt                ,grid_g(ifm)%fmapui               &
                       ,grid_g(ifm)%fmapvi               ,grid_g(ifm)%f13t                 &
                       ,grid_g(ifm)%f23t                 ,grid_g(ifm)%dxu                  &
                       ,grid_g(ifm)%dyv                  ,grid_g(ifm)%dxt                  &
                       ,grid_g(ifm)%dyt                  ,mass_g(ifm)%lnthvadv             &
                       ,mass_g(ifm)%lnthetav             )
     
   case ('THS')
      !------------------------------------------------------------------------------------!
      !   Advection part of the heating term.                                              !
      !------------------------------------------------------------------------------------!
      call storage_theta(m1,m2,m3,ifm,ia,iz,ja,jz,izu,jzv,mynum,edt                        &
                        ,basic_g(ifm)%pi0                ,basic_g(ifm)%pc                  &
                        ,scratch%vt3dp                   ,scratch%vt3dq                    &
                        ,basic_g(ifm)%theta              ,mass_g(ifm)%thvlast              &
                        ,mass_g(ifm)%lnthvtend           ,tend%pt                          )
   
   case default
      !------------------------------------------------------------------------------------!
      !   This should never happen...                                                      !
      !------------------------------------------------------------------------------------!
      call abort_run('Unexpected key'//trim(key)//'!!!','exevolve','rexev.f90')
   end select

   return
end subroutine exevolve
!==========================================================================================!
!==========================================================================================!






![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
! »»»»»»»»»»»»»»»»»»»»»»»»»» Set of subroutines for key = 'ADV' «««««««««««««««««««««««««« !
![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!



!==========================================================================================!
!==========================================================================================!
!    This subroutine finds the advection term. Note that this will work onl for terrain-   !
! following coordinates. An adaptive version could be built based on the other advection   !
! subroutines for adaptive coordinate.                                                     !
!------------------------------------------------------------------------------------------!
subroutine exadvlf(m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,itopo,rtgu,fmapui,rtgv,fmapvi,f13t    &
                  ,f23t,rtgt,fmapt,dxt,dyt,uc,dn0u,vc,dn0v,dn0,wc,pc,pt)
   use mem_grid , only : hw4 & ! intent(in)
                       , dzt ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer , intent(in)                         :: m1,m2,m3,ia,iz,ja,jz,izu,jzv,itopo,jdim
   real    , intent(in)   , dimension(m2,m3)    :: rtgu,fmapui,rtgv,fmapvi,f13t,f23t
   real    , intent(in)   , dimension(m2,m3)    :: rtgt,fmapt,dxt,dyt
   real    , intent(in)   , dimension(m1,m2,m3) :: uc,dn0u,vc,dn0v,dn0,wc,pc
   real    , intent(inout), dimension(m1,m2,m3) :: pt
   !----- Local variables -----------------------------------------------------------------!
   integer                                      :: i,j,k,im,jm
   real                                         :: c1z,c1x,c1y
   real    ,                dimension(m1,m2,m3) :: flxu,flxv,flxw
   !---------------------------------------------------------------------------------------!



   !----- Compute momentum fluxes flxu, flxv, flxw ----------------------------------------!
   do j = 1,m3
      do i = 1,m2
         do k = 1,m1
            flxu(k,i,j) = uc(k,i,j) * dn0u(k,i,j) * rtgu(i,j) * fmapui(i,j)
            flxv(k,i,j) = vc(k,i,j) * dn0v(k,i,j) * rtgv(i,j) * fmapvi(i,j)
         enddo
      enddo
   enddo
  
   if(itopo == 0) then
      do j = 1,m3
         do i = 1,m2
            do k = 1,m1-1
               flxw(k,i,j) = wc(k,i,j) * .5 * (dn0(k,i,j) + dn0(k+1,i,j))
            end do
         end do
      end do
   else
      do j = 1,m3
         jm = max(j-1,1)
         do i = 1,m2
            im = max(i-1,1)
            do k = 1,m1-1
               flxw(k,i,j) = wc(k,i,j) * .5 * (dn0(k,i,j) + dn0(k+1,i,j))                  &
                           + hw4(k) * ( ( flxu(k,i,j) + flxu(k+1,i,j)                      &
                                        + flxu(k,im,j) + flxu(k+1,im,j) ) * f13t(i,j)      &
                                      + ( flxv(k,i,j) + flxv(k+1,i,j)                      &
                                        + flxv(k,i,jm) + flxv(k+1,i,jm) ) * f23t(i,j) )
            end do
         end do
      end do
   end if
  
   !---------------------------------------------------------------------------------------!
   !  Compute advection contribution of zonal gradient to Exner function tendency.         !
   !---------------------------------------------------------------------------------------!
   do j = ja,jz
      do i = ia,izu
         c1x = 0.5 / rtgt(i,j) * fmapt(i,j) * dxt(i,j)
         do k = 2,m1-1 
            pt(k,i,j) = pt(k,i,j)                                                          &
                      - c1x / dn0(k,i,j)                                                   &
                      * ( flxu(k,i,j)   * (pc(k,i,j) + pc(k,i+1,j))                        &
                        - flxu(k,i-1,j) * (pc(k,i,j) + pc(k,i-1,j))                        &
                        - (flxu(k,i,j) - flxu(k,i-1,j)) * 2.* pc(k,i,j) )
          end do
       end do
    end do
  
   !---------------------------------------------------------------------------------------!
   !  Compute advection contribution of meridional gradient to Exner function tendency.    !
   !---------------------------------------------------------------------------------------!
   do j=ja,jzv
      do i=ia,iz
         c1y = 0.5 / rtgt(i,j) * fmapt(i,j) * dyt(i,j)
         do k=2,m1-1
            pt(k,i,j) = pt(k,i,j)                                                          &
                      - c1y /dn0(k,i,j)                                                    &
                      * ( flxv(k,i,j)      * (pc(k,i,j)+pc(k,i,j+jdim))                    &
                        - flxv(k,i,j-jdim) * (pc(k,i,j)+pc(k,i,j-jdim))                    &
                        - (flxv(k,i,j)-flxv(k,i,j-jdim)) * 2.* pc(k,i,j) )
         end do
      end do
   end do
  
   !---------------------------------------------------------------------------------------!
   !  Compute advection contribution of vertical gradient to Exner function tendency.      !
   !---------------------------------------------------------------------------------------! 
   do j=ja,jz
      do i=ia,iz
         c1z = 0.5 / rtgt(i,j)
         do k=2,m1-1
            pt(k,i,j) = pt(k,i,j)                                                          &
                      - c1z * dzt(k) /dn0(k,i,j)                                           &
                      * ( flxw(k,i,j)   * (pc(k,i,j)+pc(k+1,i,j))                          &
                        - flxw(k-1,i,j) * (pc(k,i,j)+pc(k-1,i,j))                          &
                        -  (flxw(k,i,j)-flxw(k-1,i,j)) * 2. * pc(k,i,j) )
         end do
      end do
   end do
   return
end subroutine exadvlf
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will compute the compression term of Medvigy's equation (12).         !
!------------------------------------------------------------------------------------------!
subroutine excondiv(m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,itopo,uc,vc,wc,pc,pt,dxt,dyt,rtgt    &
                   ,rtgu,rtgv,f13t,f23t,fmapt,fmapui,fmapvi )
   use rconstants , only : rocv  ! ! intent(in)
   use mem_grid   , only : hw4   & ! intent(in)
                         , dzm   ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                      , intent(in)     :: m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,itopo
   real    , dimension(   m2,m3), intent(in)     :: dxt,dyt,rtgt,rtgu,rtgv
   real    , dimension(   m2,m3), intent(in)     :: f13t,f23t,fmapt,fmapui,fmapvi
   real    , dimension(m1,m2,m3), intent(in)     :: uc,vc,wc,pc
   real    , dimension(m1,m2,m3), intent(inout)  :: pt
   !----- Local variables -----------------------------------------------------------------!
   integer                                       :: i,j,k,im,jm
   real                                          :: c1z,c1x,c1y
   real    , dimension(m1,m2,m3)                 :: flxu,flxv,flxw
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   ! Preparing fluxes. These are:  (transformed velocities) times (a) times (mapfactor)    !
   !---------------------------------------------------------------------------------------!
   do j=1,m3
      do i=1,m2
         do k=1,m1
            flxu(k,i,j) = uc(k,i,j) * rtgu(i,j) * fmapui(i,j)
            flxv(k,i,j) = vc(k,i,j) * rtgv(i,j) * fmapvi(i,j)
         end do
      end do
   end do
   !----- Vertical flux. ------------------------------------------------------------------!
   if(itopo == 0)then
      do j=1,m3
         do i=1,m2
            do k=1,m1-1
               flxw(k,i,j)=wc(k,i,j)
            end do
         end do
      end do
   else
      do j=1,m3
         jm=max(j-1,1)
         do i=1,m2
            im=max(i-1,1)
            do k=1,m1-1
               flxw(k,i,j) = wc(k,i,j)                                                     &
                           + hw4(k) * ( ( flxu(k,i,j)  + flxu(k+1,i,j)                     &
                                        + flxu(k,im,j) + flxu(k+1,im,j) ) * f13t(i,j)      &
                                      + ( flxv(k,i,j)  + flxv(k+1,i,j)                     &
                                        + flxv(k,i,jm) + flxu(k+1,i,jm) ) * f23t(i,j) )
            end do
         end do
      end do
   end if
  
   !---------------------------------------------------------------------------------------!
   !     Computing the contribution of the zonal gradient of zonal wind to Exner tendency. !
   !---------------------------------------------------------------------------------------!
   do j=ja,jz
      do i=ia,izu
         c1x=fmapt(i,j)*dxt(i,j)/rtgt(i,j)
         do k=2,m1-1
            pt(k,i,j) = pt(k,i,j) - c1x * ( flxu(k,i,j)-flxu(k,i-1,j)) * pc(k,i,j) * rocv
         end do
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !     Computing the contribution of the meridional gradient of meridional wind to Exner !
   ! tendency.                                                                             !
   !---------------------------------------------------------------------------------------!
   do j=ja,jzv
      do i=ia,iz
         c1y=fmapt(i,j)*dyt(i,j)/rtgt(i,j)
         do k=2,m1-1
            pt(k,i,j) = pt(k,i,j) - c1y * (flxv(k,i,j)-flxv(k,i,j-jdim)) * pc(k,i,j) *rocv
         end do
      end do
   end do

   !---------------------------------------------------------------------------------------!
   !     Computing the contribution of the vertical gradient of vertical velocity to Exner !
   ! tendency.                                                                             !
   !---------------------------------------------------------------------------------------!
   do j=ja,jz
      do i=ia,iz
         c1z=1.0/rtgt(i,j)
         do k=2,m1-1
            pt(k,i,j) = pt(k,i,j)                                                          &
                      - c1z * dzm(k) * (flxw(k,i,j)-flxw(k-1,i,j)) * pc(k,i,j) * rocv
         enddo
      enddo
   enddo
  
   return
end subroutine excondiv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine will save the current value of theta-v for the advection term.         !
!------------------------------------------------------------------------------------------!
subroutine fill_thvlast(m1,m2,m3,ia,iz,ja,jz,thvlast,theta,rtp,rv)
   use therm_lib, only : virtt
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                      , intent(in)    :: m1, m2, m3, ia, iz, ja, jz
   real    , dimension(m1,m2,m3), intent(in)    :: theta, rtp, rv
   real    , dimension(m1,m2,m3), intent(out)   :: thvlast
   !----- Local variables -----------------------------------------------------------------!
   integer                                      :: k, i, j
   !---------------------------------------------------------------------------------------!

   do j=ja,jz
      do i=ia,iz
         do k=1,m1
            thvlast(k,i,j)=virtt(theta(k,i,j),rv(k,i,j),rtp(k,i,j))
         end do
      end do
   end do
   return
end subroutine fill_thvlast
!==========================================================================================!
!==========================================================================================!



!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!
!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!





![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
! »»»»»»»»»»»»»»»»»»»»»»»»»» Set of subroutines for key = 'THA' «««««««««««««««««««««««««« !
![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!



!==========================================================================================!
!==========================================================================================!
!   This is the small driver to compute the advection part of the heating term.            !
!------------------------------------------------------------------------------------------!
subroutine advect_theta(m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,mynum,edt,up,uc,vp,vc,wp,wc,pi0  &
                       ,pc,pt,theta,rtp,rv,dn0,dn0u,dn0v,rtgt,rtgu,rtgv,fmapt,fmapui       &
                       ,fmapvi,f13t,f23t,dxu,dyv,dxt,dyt,lnthvadv,lnthetav)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                       , intent(in)    :: m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,mynum
   real                          , intent(in)    :: edt 
   real    , dimension(m1,m2,m3) , intent(in)    :: up,uc,vp,vc,wp,wc,pi0,pc,rtp,theta,rv
   real    , dimension(m1,m2,m3) , intent(in)    :: dn0,dn0u,dn0v
   real    , dimension(   m2,m3) , intent(in)    :: rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi
   real    , dimension(   m2,m3) , intent(in)    :: f13t,f23t,dxu,dyv,dxt,dyt
   real    , dimension(m1,m2,m3) , intent(out)   :: lnthvadv,lnthetav
   real    , dimension(m1,m2,m3) , intent(inout) :: pt
   !---------------------------------------------------------------------------------------!

   call exthvadv(m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,mynum,edt,up,uc,vp,vc,wp,wc,theta,rtp   &
                ,rv,dn0,dn0u,dn0v,rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi,f13t,f23t,dxu,dyv,dxt &
                ,dyt,lnthvadv,lnthetav)
   call exhtend_ad(m1,m2,m3,ia,iz,ja,jz,pi0,pc,pt,lnthvadv)

   return
end subroutine advect_theta
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine will compute the virtual temperature tendency due to advection        !
! (actually, its negative, since that's the contribution to the Exner function tendency).  !
!------------------------------------------------------------------------------------------!
subroutine exthvadv(m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,mynum,edt,up,uc,vp,vc,wp,wc,theta    &
                   ,rtp,rv,dn0,dn0u,dn0v,rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi,f13t,f23t,dxu  &
                   ,dyv,dxt,dyt,lnthvadv,lnthetav)
   use mem_scratch , only : scratch & ! intent(in)
                          , vctr1   & ! intent(in)
                          , vctr2   ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                      , intent(in)   :: m1,m2,m3,ia,iz,ja,jz,izu,jzv,jdim,mynum
   real                         , intent(in)   :: edt
   real    , dimension(m1,m2,m3), intent(in)   :: up,uc,vp,vc,wp,wc,theta,rtp,rv
   real    , dimension(m1,m2,m3), intent(in)   :: dn0,dn0u,dn0v
   real    , dimension(   m2,m3), intent(in)   :: rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi
   real    , dimension(   m2,m3), intent(in)   :: f13t,f23t,dxu,dyv,dxt,dyt
   real    , dimension(m1,m2,m3), intent(out)  :: lnthvadv,lnthetav
   !----- Local variables -----------------------------------------------------------------!
   integer                                     :: i,j,k,isiz
   !---------------------------------------------------------------------------------------!

   !----- Finding the total size of matrices ----------------------------------------------!
   isiz=m1*m2*m3
   call azero3(m1*m2*m3,scratch%vt3da,scratch%vt3db,scratch%vt3dc)
   call prep_timeave(m1,m2,m3,edt,up,uc,vp,vc,wp,wc,scratch%vt3da,scratch%vt3db            &
                    ,scratch%vt3dc)
  
   call prep_lnthetv(m1,m2,m3,ia,iz,ja,jz,theta,rtp,rv,lnthetav)
  
   call fa_preptc(m1,m2,m3,scratch%vt3da,scratch%vt3db,scratch%vt3dc,scratch%vt3dd         &
                 ,scratch%vt3de,scratch%vt3df,scratch%vt3dh,scratch%vt3di,scratch%vt3dj    &
                 ,scratch%vt3dk,dn0,dn0u,dn0v,rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi,f13t,f23t &
                 ,dxu,dyv,dxt,dyt,mynum)
   call atob(m1*m2*m3,lnthetav,scratch%scr1)


   call fa_xc(m1,m2,m3,ia,iz,1,m3,lnthetav,scratch%scr1,scratch%vt3da,scratch%vt3dd        &
             ,scratch%vt3dg,scratch%vt3dh,scratch%vt3di,mynum)

   if (jdim == 1)                                                                          &
      call fa_yc(m1,m2,m3,ia,iz,ja,jz,lnthetav,scratch%scr1,scratch%vt3db,scratch%vt3de    &
                ,scratch%vt3dg,scratch%vt3dj,scratch%vt3di,jdim,mynum)
  
   call fa_zc(m1,m2,m3,ia,iz,ja,jz,lnthetav,scratch%scr1,scratch%vt3dc,scratch%vt3df       &
             ,scratch%vt3dg,scratch%vt3dk,vctr1,vctr2,mynum)

   call azero(m1*m2*m3,lnthvadv)
   call advtndc(m1,m2,m3,ia,iz,ja,jz,lnthetav,scratch%scr1,lnthvadv,edt,mynum)
   
   !---------------------------------------------------------------------------------------!
   !     Switching the sign of the advection term... This is because we need to total      !
   ! derivative of theta-V, which is the local change - advection.                         !
   !---------------------------------------------------------------------------------------!
   do j=1,m3
      do i=1,m2
         do k=1,m1
            lnthvadv(k,i,j)=-1.0*lnthvadv(k,i,j)
         end do
      end do
   end do 
  
   return
end subroutine exthvadv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    Finding the mid-point between past and present for advection terms.                   !
!------------------------------------------------------------------------------------------!
subroutine prep_timeave(m1,m2,m3,edt,up,uc,vp,vc,wp,wc,vt3da,vt3db,vt3dc)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                       , intent(in)  :: m1,m2,m3
   real                          , intent(in)  :: edt
   real    , dimension(m1,m2,m3) , intent(in)  :: up,uc,vp,vc,wp,wc
   real    , dimension(m1,m2,m3) , intent(out) :: vt3da,vt3db,vt3dc
   !----- Local variables -----------------------------------------------------------------!
   integer                                     :: i,j,k
   !---------------------------------------------------------------------------------------!
   do j=1,m3
      do i=1,m2
         do k=1,m1
            vt3da(k,i,j) = 0.5 * edt * (up(k,i,j) + uc(k,i,j))
            vt3db(k,i,j) = 0.5 * edt * (vp(k,i,j) + vc(k,i,j))
            vt3dc(k,i,j) = 0.5 * edt * (wp(k,i,j) + wc(k,i,j))
         end do
      end do
   end do
   return
end subroutine prep_timeave
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    Finding the log of virtual potential temperature.                                     !
!------------------------------------------------------------------------------------------!
subroutine prep_lnthetv(m1,m2,m3,ia,iz,ja,jz,theta,rtp,rv,lnthetav)
   use therm_lib , only : virtt
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                       , intent(in)  :: m1,m2,m3,ia,iz,ja,jz
   real    , dimension(m1,m2,m3) , intent(in)  :: theta,rtp,rv
   real    , dimension(m1,m2,m3) , intent(out) :: lnthetav
   !----- Local variables -----------------------------------------------------------------!
   integer                                     :: i,j,k
   !---------------------------------------------------------------------------------------!
   do j=1,m3
      do i=1,m2
         do k=1,m1
              lnthetav(k,i,j)=log(virtt(theta(k,i,j),rv(k,i,j),rtp(k,i,j)))
         end do
      end do
   end do
   return
end subroutine prep_lnthetv
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    Adding the contribution of advection in the heating term to the Exner function        !
! tendency.                                                                                !
!------------------------------------------------------------------------------------------!
subroutine exhtend_ad(m1,m2,m3,ia,iz,ja,jz,pi0,pc,pt,lnthvadv)
   use rconstants , only : rocv  ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                       , intent(in)    :: m1,m2,m3,ia,iz,ja,jz
   real    , dimension(m1,m2,m3) , intent(in)    :: pi0,pc,lnthvadv
   real    , dimension(m1,m2,m3) , intent(inout) :: pt
   !----- Local variables -----------------------------------------------------------------!
   integer                                       :: i,j,k
   !---------------------------------------------------------------------------------------!
  
   do j=ja,jz
      do i=ia,iz
         do k=2,m1-1
            pt(k,i,j) = pt(k,i,j) + rocv * (pi0(k,i,j) + pc(k,i,j)) * lnthvadv(k,i,j)
         end do
      end do
   end do
  
   return
end subroutine exhtend_ad
!==========================================================================================!
!==========================================================================================!



!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!
!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!





![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
! »»»»»»»»»»»»»»»»»»»»»»»»»» Set of subroutines for key = 'THS' «««««««««««««««««««««««««« !
![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!



!==========================================================================================!
!==========================================================================================!
!    This is a mini-driver for the tendency part of the heating term.                      !
!------------------------------------------------------------------------------------------!
subroutine storage_theta(m1,m2,m3,ifm,ia,iz,ja,jz,izu,jzv,mynum,edt,pi0,pc,rtp,rv,theta    &
                        ,thvlast,lnthvtend,pt)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                       , intent(in)    :: m1,m2,m3,ifm,ia,iz,ja,jz,izu,jzv,mynum
   real                          , intent(in)    :: edt
   real    , dimension(m1,m2,m3) , intent(in)    :: pi0,pc,rtp,rv,theta,thvlast
   real    , dimension(m1,m2,m3) , intent(out)   :: lnthvtend
   real    , dimension(m1,m2,m3) , intent(inout) :: pt
   !---------------------------------------------------------------------------------------!
  
   !----- Computing the tendency of theta-V -----------------------------------------------!
   call prep_lnthvtend(m1,m2,m3,ifm,ia,iz,ja,jz,izu,jzv,edt,theta,thvlast,rtp,rv,lnthvtend)
   !----- Computing the tendency part of the heating term ---------------------------------!
   call exhtend_st(m1,m2,m3,ia,iz,ja,jz,pi0,pc,rtp,theta,lnthvtend,rv,pt)
   
   return
end subroutine storage_theta
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine computes the local time derivative of theta-V.                           !
!------------------------------------------------------------------------------------------!
subroutine prep_lnthvtend(m1,m2,m3,ifm,ia,iz,ja,jz,izu,jzv,edt,theta,thvlast,rtp,rv        &
                       ,lnthvtend)
   use mem_grid  , only : time     & ! intent(in)
                        , dtlongn  ! ! intent(in)
   use therm_lib , only : virtt    ! ! Function
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                       , intent(in)  :: m1,m2,m3,ifm,ia,iz,ja,jz,izu,jzv
   real                          , intent(in)  :: edt
   real    , dimension(m1,m2,m3) , intent(in)  :: theta,rv,rtp,thvlast
   real    , dimension(m1,m2,m3) , intent(out) :: lnthvtend
   !----- Local variables -----------------------------------------------------------------!
   integer                                     :: i,j,k
   real                                        :: edti
   !---------------------------------------------------------------------------------------!
  
   !----- First time the subroutine is called in this node, leave it zero -----------------!
   if (time < dtlongn(ifm)) then
      do j=ja,jz
         do i=ia,iz
            do k=2,m1-1
               lnthvtend(k,i,j) = 0.0
            end do
         end do
      end do
   else
      edti=1.0/edt
  
      do j=ja,jz
         do i=ia,iz
            do k=2,m1-1
               lnthvtend(k,i,j) = edti * log( virtt(theta(k,i,j),rv(k,i,j),rtp(k,i,j))     &
                                            / thvlast(k,i,j))
            end do
         end do
      end do
   end if
  
   return
end subroutine prep_lnthvtend
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This will compute the tendency part of the heating term.                               !
!------------------------------------------------------------------------------------------!
subroutine exhtend_st(m1,m2,m3,ia,iz,ja,jz,pi0,pc,rtp,theta,lnthvtend,rv,pt)
   use rconstants , only : rocv  ! ! intent(in)
   use therm_lib  , only : virtt ! ! Function
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer                       , intent(in)    :: m1,m2,m3,ia,iz,ja,jz
   real    , dimension(m1,m2,m3) , intent(in)    :: pi0,pc,rtp,theta,lnthvtend,rv
   real    , dimension(m1,m2,m3) , intent(inout) :: pt
   !----- Local variables -----------------------------------------------------------------!
   integer                                       :: i,j,k
   !---------------------------------------------------------------------------------------!

   do j=ja,jz
      do i=ia,iz
         do k=2,m1-1
            !----- LnThvtend is the log of the theta-v tendency ---------------------------!
            pt(k,i,j) = pt(k,i,j) + rocv * (pi0(k,i,j) + pc(k,i,j)) * lnthvtend(k,i,j)
         end do
      end do
   end do

   return
end subroutine exhtend_st
!==========================================================================================!
!==========================================================================================!
