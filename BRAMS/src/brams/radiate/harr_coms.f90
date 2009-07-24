!====================================== Change Log ========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This module contains some scratch variables used in Harrington et al. (2000)          !
! radiation code.                                                                          !
!------------------------------------------------------------------------------------------!
module harr_coms
   use mem_harr, only: mg
   implicit none

   integer                           :: nradmax ! Maximum nrad possible.

   !---------------------------------------------------------------------------------------!
   !    Set activation flags for gases of importance: H2O, CO2, and O3, respectively.      !
   ! Flag = 1: gas active;    Flag = 0: gas not active                                     !
   !---------------------------------------------------------------------------------------!
   integer, dimension(mg), parameter :: ngass=(/1, 1, 1/)  ! Flags for  shortwave
   integer, dimension(mg), parameter :: ngast=(/1, 1, 1/)  ! Flags for for longwave
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Properties derived from default settings in the microphysics. Since we only need   !
   ! the reciprocal, only this will be assigned. Parm is seldom assigned to pristine ice   !
   ! since this variable must be either off or prognosed. Anyway, this is just for a       !
   ! simple estimate for the parametrised cloud.                                           !
   !---------------------------------------------------------------------------------------!
   real                  , parameter :: parmi_cloud=1./.3e9
   real                  , parameter :: parmi_prist=1./.1e4
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !     DN limiters for consistency with fit coefficients (note that these values are not !
   !  the same as in micro).                                                               !
   !---------------------------------------------------------------------------------------!
   real,dimension(7),parameter :: dnmin= (/   1.,   10.,   1.,  125.,   10.,   10.,   10./)
   real,dimension(7),parameter :: dnmax= (/1000.,10000., 125.,10000.,10000.,10000.,10000./)

   !---------------------------------------------------------------------------------------!
   !     Coefficients for assymetry parameter.                                             !
   !---------------------------------------------------------------------------------------!
   real,dimension(7),parameter :: sacoef= (/ 1.13, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08 /)
   real,dimension(7),parameter :: lacoef= (/ 2.10, 2.90, 2.90, 2.90, 2.90, 2.90, 2.90 /)

   !---------------------------------------------------------------------------------------!
   !    Array kradcat maps RAMS/OLAM microphysics hydrometeor categories to those          !
   ! represented in Harrington radiation code according to the following numbering:        !
   !                                                                                       !
   !     Harrington radiation code             Microphysics                                !
   ! ----------------------------------------------------------------                      !
   !  1:   cloud drops                 1.  cloud drops                                     !
   !  2:   rain                        2.  rain                                            !
   !  3:   pristine ice columns        3.  pristine ice columns                            !
   !  4:   pristine ice rosettes       4.  snow columns                                    !
   !  5:   pristine ice plates         5.  aggregates                                      !
   !  6:   snow columns                6.  graupel                                         !
   !  7:   snow rosettes               7.  hail                                            !
   !  8:   snow plates                 8.  pristine ice hexagonal plates                   !
   !  9:   aggregates columns          9.  pristine ice dendrites                          !
   !  10:  aggregates rosettes        10.  pristine ice needles                            !
   !  11:  aggregates plates          11.  pristine ice rosettes                           !
   !  12:  graupel                    12.  snow hexagonal plates                           !
   !  13:  hail                       13.  snow dendrites                                  !
   !                                  14.  snow needles                                    !
   !                                  15.  snow rosettes                                   !
   !---------------------------------------------------------------------------------------!
   integer, parameter, dimension(15) :: kradcat = (/1,2,3,6,10,12,13,5,5,3,4,8,8,6,7/)



   !----- Some parameters found in swrad/lwrad... -----------------------------------------!
   real, parameter :: top=1800.,tm=1800./296.,gma=0.002,tr=296.,ccc=6.08
   !---------------------------------------------------------------------------------------!


   !----- Dimension nrad ------------------------------------------------------------------!
   real, dimension(:)  , allocatable :: rl    ! Vapour density                     [ kg/m³]
   real, dimension(:)  , allocatable :: dzl   ! Delta-z                            [     m]
   real, dimension(:)  , allocatable :: dl    ! Air density                        [ kg/m³]
   real, dimension(:)  , allocatable :: pl    ! Pressure                           [    Pa]
   real, dimension(:)  , allocatable :: co2l  ! CO2 density                        [ kg/m³]
   real, dimension(:)  , allocatable :: o3l   ! Calculated ozone profile           [  g/m³]
   real, dimension(:)  , allocatable :: vp    ! Vapour pressure                    [    Pa]
   real, dimension(:)  , allocatable :: zml   ! Heights of W points                [     m]
   real, dimension(:)  , allocatable :: ztl   ! Heights of T points                [     m]
   real, dimension(:)  , allocatable :: tl    ! Temperature                        [     K]
   real, dimension(:)  , allocatable :: flxus ! Total upwelling shortwave flux     [  W/m²]
   real, dimension(:)  , allocatable :: flxds ! Total downwelling shortwave flux   [  W/m²]
   real, dimension(:)  , allocatable :: flxul ! Total upwelling longwave flux      [  W/m²]
   real, dimension(:)  , allocatable :: flxdl ! Total downwelling longwave flux    [  W/m²]
   !---------------------------------------------------------------------------------------!



   !----- Dimension (nrad,mg) -------------------------------------------------------------!
   real, dimension(:,:), allocatable :: u     ! path-length for H2O, CO2, O3       [    Pa]
   !---------------------------------------------------------------------------------------!
  
  
  
   !----- Dimension (nrad,mb) -------------------------------------------------------------!
   real, dimension(:,:), allocatable :: tp    ! optical depth of hydrometeors      [   1/m]
   real, dimension(:,:), allocatable :: omgp  ! Single scatter albedo of hydrom.   [   ---]
   real, dimension(:,:), allocatable :: gp    ! Asymmetry factor of hydrometeors   [   ---]
   !---------------------------------------------------------------------------------------!



   !----- Dimension(nrad,mpb) -------------------------------------------------------------!
   real, allocatable, dimension(:,:) :: fu    ! Dpwelling fluxes for pseudo-bands  [  W/m²]
   real, allocatable, dimension(:,:) :: fd    ! Downwelling flx. for pseudo-bands  [  W/m²]
   !---------------------------------------------------------------------------------------!



   !----- Dimension(m1) -------------------------------------------------------------------!
   real   , allocatable, dimension(:) :: tairk    ! Air temperature               [      K]
   real   , allocatable, dimension(:) :: rhoi     ! Vapour specific volume        [  m³/kg]
   real   , allocatable, dimension(:) :: rhoe     ! Air density                   [  kg/m³]
   real   , allocatable, dimension(:) :: rhov     ! Vapour density                [  kg/m³]
   real   , allocatable, dimension(:) :: press    ! Atmospheric pressure          [     Pa]
   real   , allocatable, dimension(:) :: rcl_parm ! Cumulus liquid mixing ratio   [  kg/kg]
   real   , allocatable, dimension(:) :: rpl_parm ! Cumulus ice mixing rati       [  kg/kg]
   integer, allocatable, dimension(:) :: nsharr   ! Habit table class: moisture
   integer, allocatable, dimension(:) :: ntharr   ! Habit table class: temperature
   !---------------------------------------------------------------------------------------!



   !----- Dimension(m1,ncat) --------------------------------------------------------------!
   integer, allocatable, dimension(:,:) :: jhcatharr ! Hydrom habit category
   real   , allocatable, dimension(:,:) :: rxharr    ! hydrom bulk mixing ratio   [  kg/kg]
   real   , allocatable, dimension(:,:) :: cxharr    ! hydrom bulk number         [   #/kg]
   real   , allocatable, dimension(:,:) :: embharr   ! hydrom mean particle mass  [   kg/#]
   !---------------------------------------------------------------------------------------!


   
   contains
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine allocates some specific scratch arrays used by LEAF-3.              !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_harr_scratch(m1,ncat,nrad,mgl,mb,mpb)
      implicit none
      !----- Arguments -------------------------------------------------------------------!
      integer, intent(in) :: m1,ncat,nrad,mgl,mb,mpb
      !-----------------------------------------------------------------------------------!
      allocate (rl       (nrad    ))
      allocate (dzl      (nrad    ))
      allocate (dl       (nrad    ))
      allocate (pl       (nrad    ))
      allocate (o3l      (nrad    ))
      allocate (co2l     (nrad    ))
      allocate (vp       (nrad    ))
      allocate (zml      (nrad    ))
      allocate (ztl      (nrad    ))
      allocate (tl       (nrad    ))
      allocate (flxus    (nrad    ))
      allocate (flxds    (nrad    ))
      allocate (flxul    (nrad    ))
      allocate (flxdl    (nrad    ))
      allocate (u        (nrad,mgl))
      allocate (tp       (nrad,mb ))
      allocate (omgp     (nrad,mb ))
      allocate (gp       (nrad,mb ))
      allocate (fu       (nrad,mpb))
      allocate (fd       (nrad,mpb))
      allocate (tairk    (m1)      )
      allocate (rhoi     (m1)      )
      allocate (rhoe     (m1)      )
      allocate (rhov     (m1)      )
      allocate (press    (m1)      )
      allocate (rcl_parm (m1)      )
      allocate (rpl_parm (m1)      )
      allocate (nsharr   (m1)      )
      allocate (ntharr   (m1)      )
      allocate (jhcatharr(m1,ncat) )
      allocate (rxharr   (m1,ncat) )
      allocate (cxharr   (m1,ncat) )
      allocate (embharr  (m1,ncat) )
       
      return
   end subroutine alloc_harr_scratch
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !   This subroutine allocates some specific scratch arrays used by LEAF-3.              !
   !---------------------------------------------------------------------------------------!
   subroutine zero_harr_scratch(m1,ncat,nrad,mgl,mb,mpb)
      implicit none

      !----- Arguments -------------------------------------------------------------------!
      integer, intent(in) :: m1,ncat,nrad,mgl,mb,mpb
      !----- Local variables -------------------------------------------------------------!
      integer             :: k,icat,ig,ib,ipb
      !-----------------------------------------------------------------------------------!
      do k=1,nrad
         rl   (k) = 0.0
         dzl  (k) = 0.0
         dl   (k) = 0.0
         pl   (k) = 0.0
         o3l  (k) = 0.0
         co2l (k) = 0.0
         vp   (k) = 0.0
         zml  (k) = 0.0
         ztl  (k) = 0.0
         tl   (k) = 0.0
         flxus(k) = 0.0
         flxds(k) = 0.0
         flxul(k) = 0.0
         flxdl(k) = 0.0
      end do
      
      do ig=1,mgl
         do k=1,nrad
            u(k,ig) = 0.0
         end do
      end do
      
      do ib=1,mb
         do k=1,nrad
            tp   (k,ib)  = 0.0
            omgp (k,ib)  = 0.0
            gp   (k,ib)  = 0.0
         end do
      end do
      
      do ipb=1,mpb
         do k=1,nrad
            fu   (k,ipb)  = 0.0
            fd   (k,ipb)  = 0.0
         end do
      end do

      do k=1,m1
         tairk    (k) = 0.0
         rhoi     (k) = 0.0
         rhoe     (k) = 0.0
         rhov     (k) = 0.0
         press    (k) = 0.0
         rcl_parm (k) = 0.0
         rpl_parm (k) = 0.0
         nsharr   (k) =   1
         ntharr   (k) =   1
      end do

      do icat=1,ncat
         do k=1,m1
            jhcatharr(k,icat) = 0
            rxharr   (k,icat) = 0.
            cxharr   (k,icat) = 0.
            embharr  (k,icat) = 0.
         end do
      end do

      return
   end subroutine zero_harr_scratch
   !=======================================================================================!
   !=======================================================================================!
end module harr_coms
!==========================================================================================!
!==========================================================================================!

