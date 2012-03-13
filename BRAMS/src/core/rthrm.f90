!===================================== Change Log =========================================!
! 5.0.0                                                                                    !
!                                                                                          !
!==========================================================================================!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!==========================================================================================!
subroutine thermo_boundary_driver(time,dtlong,f_thermo_e,f_thermo_w,f_thermo_s,f_thermo_n  &
                                 ,nzp,mxp,myp,jdim)
   !---------------------------------------------------------------------------------------!
   !    This subroutine is just a wrapper that decides whether the thermodynamic calls     !
   ! should be done at the boundaries or not.                                              !
   !---------------------------------------------------------------------------------------!
   implicit none
   !----- Input arguments -----------------------------------------------------------------!
   real        , intent(in) :: dtlong     ! Long time step                         [     s]
   real(kind=8), intent(in) :: time       ! Current model time                     [     s]
   logical     , intent(in) :: f_thermo_e ! I will run thermo for Eastern  bndry   [  ----]
   logical     , intent(in) :: f_thermo_w ! I will run thermo for Western  bndry   [  ----]
   logical     , intent(in) :: f_thermo_s ! I will run thermo for Southern bndry   [  ----]
   logical     , intent(in) :: f_thermo_n ! I will run thermo for Northern bndry   [  ----]
   integer     , intent(in) :: nzp        ! # of points in Z                       [  ----]
   integer     , intent(in) :: mxp        ! # of points in X                       [  ----]
   integer     , intent(in) :: myp        ! # of points in Y                       [  ----]
   integer     , intent(in) :: jdim       ! Flag for 2-D or 3-D run                [  ----]

   !------ Checking whether it's time to call thermo, if not, return... -------------------!
   if (mod(time,dble(dtlong))/=0) return

  
   !------ Checking longitudinal borders --------------------------------------------------!
   if (f_thermo_e) call thermo(nzp, mxp, myp, 1,   1,   1, myp)
   if (f_thermo_w) call thermo(nzp, mxp, myp, mxp, mxp, 1, myp)

   !------ Checking latitudinal borders if it is a 3-D run --------------------------------!
   if (jdim==1) then
      if (f_thermo_s) call thermo(nzp, mxp, myp, 1, mxp, 1,   1  )
      if (f_thermo_n) call thermo(nzp, mxp, myp, 1, mxp, myp, myp)
   end if

   return
end subroutine thermo_boundary_driver
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This function is also a wrapper, which will find some thermodynamic variables based  !
! on the microphysics complexity level.                                                    !
!------------------------------------------------------------------------------------------!
subroutine thermo(mzp,mxp,myp,ia,iz,ja,jz)

   use mem_grid   , only : ngrid    ! ! intent(in)    - Current grid
   use mem_basic  , only : basic_g  ! ! intent(inout) - The "basic" variables
   use mem_micro  , only : micro_g  ! ! intent(inout) - The hydrometeors
   use mem_scratch, only : scratch  & ! intent(out)   - Scratch structure, for scratch...
                         , vctr5    & ! intent(out)   - Scratch vector, for scratch...
                         , vctr6    ! ! intent(out)   - Scratch vector, for scratch...
   use therm_lib  , only : level    ! ! intent(in)    - Number of H2O phases 
   use micphys    , only : availcat ! ! intent(in)    - Hydrometeor is available [T|F]
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in)  :: mzp    ! # of points in Z                              [   ----]
   integer, intent(in)  :: mxp    ! # of points in X                              [   ----]
   integer, intent(in)  :: myp    ! # of points in Y                              [   ----]
   integer, intent(in)  :: ia     ! Node Western edge                             [   ----]
   integer, intent(in)  :: iz     ! Node Eastern edge                             [   ----]
   integer, intent(in)  :: ja     ! Node Southern edge                            [   ----]
   integer, intent(in)  :: jz     ! Node Northern edge                            [   ----]
   !----- Local variables -----------------------------------------------------------------!
   integer              :: mzxyp  ! # of points                                   [   ----]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !   Calling the appropriate driver according to the microphysics complexity level (aka  !
   ! how many phases H20 will be allowed to exist).                                        !
   !---------------------------------------------------------------------------------------!
   select case (level)
   case (0)
   !----- No water substance, copying theta-il to theta. ----------------------------------!
      call drythrm(mzp,mxp,myp,ia,iz,ja,jz                                                 &
           ,basic_g(ngrid)%thp                     ,basic_g(ngrid)%theta                   )

   !----- No condensation, if super-saturation happens, then it will be supersaturated... -!
   case (1)
      call vapthrm(mzp,mxp,myp,ia,iz,ja,jz                                                 &
           ,basic_g(ngrid)%thp                     ,basic_g(ngrid)%theta                   &
           ,basic_g(ngrid)%rtp                     ,basic_g(ngrid)%rv                      )

   !----- Liquid phase only: cloud droplets can develop, but no ice -----------------------!
   case (2)
      call satadjst(mzp,mxp,myp,ia,iz,ja,jz                                                &
           ,basic_g(ngrid)%pp                      ,scratch%scr1                           &
           ,basic_g(ngrid)%thp                     ,basic_g(ngrid)%theta                   &
           ,scratch%vt3db                          ,basic_g(ngrid)%pi0                     &
           ,basic_g(ngrid)%rtp                     ,basic_g(ngrid)%rv                      &
           ,micro_g(ngrid)%rcp                     ,scratch%scr2                           )

   !----- All three phases of water are allowed -------------------------------------------!
   case (3)
      mzxyp = mzp*mxp*myp
      call azero(mzxyp,scratch%vt3da)
      call azero(mzxyp,scratch%vt3db)
      call azero(mzxyp,scratch%vt3dc)
      call azero(mzxyp,scratch%vt3dd)
      call azero(mzxyp,scratch%vt3de)
      call azero(mzxyp,scratch%vt3df)
      call azero(mzxyp,scratch%vt3dg)
      call azero(mzxyp,scratch%vt3dh)
      call azero(mzxyp,scratch%vt3di)
      if (availcat(1)) call atob(mzxyp,micro_g(ngrid)%rcp,scratch%vt3da)
      if (availcat(2)) call atob(mzxyp,micro_g(ngrid)%rrp,scratch%vt3db)
      if (availcat(3)) call atob(mzxyp,micro_g(ngrid)%rpp,scratch%vt3dc)
      if (availcat(4)) call atob(mzxyp,micro_g(ngrid)%rsp,scratch%vt3dd)
      if (availcat(5)) call atob(mzxyp,micro_g(ngrid)%rap,scratch%vt3de)
      if (availcat(6)) then
         call atob(mzxyp,micro_g(ngrid)%rgp,scratch%vt3df)
         call atob(mzxyp,micro_g(ngrid)%q6,scratch%vt3dh)
      end if
      if (availcat(7)) then
         call atob(mzxyp,micro_g(ngrid)%rhp,scratch%vt3dg)
         call atob(mzxyp,micro_g(ngrid)%q7,scratch%vt3di)
      end if
      call wetthrm3(mzp,mxp,myp,ia,iz,ja,jz,availcat                                       &
           ,basic_g(ngrid)%pi0                     ,basic_g(ngrid)%pp                      &
           ,basic_g(ngrid)%thp                     ,basic_g(ngrid)%theta                   &
           ,basic_g(ngrid)%rtp                     ,basic_g(ngrid)%rv                      &
           ,scratch%vt3da                          ,scratch%vt3db                          &
           ,scratch%vt3dc                          ,scratch%vt3dd                          &
           ,scratch%vt3de                          ,scratch%vt3df                          &
           ,scratch%vt3dg                          ,scratch%vt3dh                          &
           ,scratch%vt3di                          ,vctr5 ,vctr6                           )

   case default
      call abort_run('Thermo option not supported...LEVEL out of bounds'                   &
                    ,'thermo','rthrm.f90')

   end select

   return
end subroutine thermo
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine drythrm(m1,m2,m3,ia,iz,ja,jz,thil,theta)
  !----------------------------------------------------------------------------------------!
  !    This routine calculates theta and rv for the case where no condensate is allowed.   !
  !----------------------------------------------------------------------------------------!

   implicit none
   !----- Arguments: ----------------------------------------------------------------------!
   integer                  , intent(in)    :: m1,m2,m3 ! Grid dimensions          [  ----]
   integer                  , intent(in)    :: ia,iz    ! Lateral boundaries       [  ----]
   integer                  , intent(in)    :: ja,jz    ! Lateral boundaries       [  ----]
   real, dimension(m1,m2,m3), intent(in)    :: thil     ! Ice-liquid pot. temper.  [     K]
   real, dimension(m1,m2,m3), intent(inout) :: theta    ! Potential temperature    [     K]
   !----- Local Variables -----------------------------------------------------------------!
   integer                                  :: i,j,k    ! Counters
   !---------------------------------------------------------------------------------------!
  
   !---------------------------------------------------------------------------------------!
   !   Here it's just a matter of copying the "potentially wet" variables to the dry ones, !
   ! they are all the same...                                                              !
   !---------------------------------------------------------------------------------------!
   do j = ja,jz
      do i = ia,iz
         do k = 1,m1
            theta(k,i,j) = thil(k,i,j)
         end do
      end do
   end do
   return
end subroutine drythrm
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine vapthrm(m1,m2,m3,ia,iz,ja,jz,thil,theta,rt,rv)
  !----------------------------------------------------------------------------------------!
  !    This routine calculates theta and rv for the case where no condensate is allowed.   !
  !----------------------------------------------------------------------------------------!

   implicit none
   !----- Arguments: ----------------------------------------------------------------------!
   integer                  , intent(in)    :: m1,m2,m3 ! Grid dimensions          [  ----]
   integer                  , intent(in)    :: ia,iz    ! Lateral boundaries       [  ----]
   integer                  , intent(in)    :: ja,jz    ! Lateral boundaries       [  ----]
   real, dimension(m1,m2,m3), intent(in)    :: thil     ! Ice-liquid pot. temper.  [     K]
   real, dimension(m1,m2,m3), intent(in)    :: rt       ! Total mixing ratio       [ kg/kg]
   real, dimension(m1,m2,m3), intent(inout) :: theta    ! Potential temperature    [     K]
   real, dimension(m1,m2,m3), intent(inout) :: rv       ! Vapour mixing ratio      [ kg/kg]
   !----- Local Variables -----------------------------------------------------------------!
   integer                                  :: i,j,k    ! Counters
   !---------------------------------------------------------------------------------------!
  
   !---------------------------------------------------------------------------------------!
   !   Here it's just a matter of copying the "potentially wet" variables to the dry ones, !
   ! they are all the same...                                                              !
   !---------------------------------------------------------------------------------------!
   do j = ja,jz
      do i = ia,iz
         do k = 1,m1
            theta(k,i,j) = thil(k,i,j)
            rv(k,i,j)    = rt(k,i,j)
         end do
      end do
   end do
   return
end subroutine vapthrm
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This routine diagnoses theta, rv, and rcp using a saturation adjustment for the case !
! when water is in the liquid and vapour phases only.                                      !
!------------------------------------------------------------------------------------------!
subroutine satadjst(m1,m2,m3,ia,iz,ja,jz,pp,p,thil,theta,t,pi0,rtp,rv,rcp,rvls)
   use therm_lib , only : rslf          & ! function
                        , alvl          & ! function
                        , exner2press   & ! function
                        , extheta2temp  & ! function
                        , extemp2theta  & ! function
                        , thil2tqliq    ! ! subroutine
   implicit none
  
   !----- Input arguments -----------------------------------------------------------------!
   integer                  , intent(in)    :: m1,m2,m3 ! Grid dimensions         [   ----]
   integer                  , intent(in)    :: ia,iz    ! Lateral boundaries      [   ----]
   integer                  , intent(in)    :: ja,jz    ! Lateral boundaries      [   ----]
   real, dimension(m1,m2,m3), intent(in)    :: pp       ! Exner perturbation      [ J/kg/K]
   real, dimension(m1,m2,m3), intent(in)    :: thil     ! Ice-liquid pot. temper. [      K]
   real, dimension(m1,m2,m3), intent(in)    :: pi0      ! Ref. Exner function     [ J/kg/K]
   real, dimension(m1,m2,m3), intent(in)    :: rtp      ! Total mixing ratio      [  kg/kg]
   !----- Output arguments ----------------------------------------------------------------!
   real, dimension(m1,m2,m3), intent(inout) :: p        ! Pressure                [     Pa]
   real, dimension(m1,m2,m3), intent(inout) :: t        ! Temperature             [      K]
   real, dimension(m1,m2,m3), intent(inout) :: rvls     ! Saturation mixing ratio [  kg/kg]
   real, dimension(m1,m2,m3), intent(inout) :: rcp      ! Liquid water mix. ratio [  kg/kg]
   real, dimension(m1,m2,m3), intent(inout) :: theta    ! Potential temperature   [      K]
   real, dimension(m1,m2,m3), intent(inout) :: rv       ! Vapour mixing ratio     [  kg/kg]
   !----- Local variables -----------------------------------------------------------------!
   integer             :: i,j,k
   real                :: exner
   !---------------------------------------------------------------------------------------!
  
   do j = ja,jz
      do i = ia,iz
         do k = 1,m1
            exner = (pi0(k,i,j) + pp(k,i,j)) 
            p(k,i,j) = exner2press(exner)
            !----- First guess for temperature and liquid mixing ratio --------------------!
            t(k,i,j)   = extheta2temp(exner,thil(k,i,j))
            rcp(k,i,j) = max(0.,rtp(k,i,j) - rslf(p(k,i,j),t(k,i,j)))
            !----- Adjusting the accordingly to the saturation point ----------------------!
            call thil2tqliq(thil(k,i,j),exner,p(k,i,j),rtp(k,i,j),rcp(k,i,j),t(k,i,j)      &
                           ,rv(k,i,j),rvls(k,i,j))
            theta(k,i,j) = extemp2theta(exner,t(k,i,j))
         end do
      end do
   end do

   return
end subroutine satadjst
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This routine calculates theta and rv for "level 3 microphysics" given prognosed       !
! theta_il, cloud, rain, pristine ice, snow, aggregates, graupel, hail, q6 (internal       !
! energy of graupel), and q7 (internal energy of hail).                                    !
!------------------------------------------------------------------------------------------!
subroutine wetthrm3(m1,m2,m3,ia,iz,ja,jz,availcat,pi0,pp,thp,theta,rtp,rv,rcp,rrp,rpp,rsp  &
                   ,rap,rgp,rhp,q6,q7,rliq,rice)

   use rconstants, only : ttripoli     & ! intent(in)
                        , alvl3        & ! intent(in)
                        , alvi3        & ! intent(in)
                        , cpdryi4      & ! intent(in)
                        , htripolii    ! ! intent(in)
   use node_mod  , only : mynum        ! ! intent(in)
   use therm_lib , only : thil2temp    & ! function
                        , exner2press  & ! function
                        , extheta2temp & ! function
                        , extemp2theta ! ! function

   implicit none

   !----- Input arguments -----------------------------------------------------------------!
   integer                     , intent(in)    :: m1       ! Vertical dimension   [   ----]
   integer                     , intent(in)    :: m2       ! Zonal dimension      [   ----]
   integer                     , intent(in)    :: m3       ! Meridional dimension [   ----]
   integer                     , intent(in)    :: ia       ! Node Western edge    [   ----]
   integer                     , intent(in)    :: iz       ! Node Eastern edge    [   ----]
   integer                     , intent(in)    :: ja       ! Node Southern edge   [   ----]
   integer                     , intent(in)    :: jz       ! Node Northern edge   [   ----]
   logical, dimension(*)       , intent(in)    :: availcat ! Microphysics handle  [    T|F]
   real   , dimension(m1,m2,m3), intent(in)    :: pi0      ! Ref. Exner function  [ J/kg/K]
   real   , dimension(m1,m2,m3), intent(in)    :: pp       ! Exner funtion pert.  [ J/kg/K]
   real   , dimension(m1,m2,m3), intent(in)    :: thp      ! Ice-liq. pot. temp.  [      K]
   real   , dimension(m1,m2,m3), intent(in)    :: rtp      ! Total mixing ratio   [  kg/kg]
   real   , dimension(m1,m2,m3), intent(in)    :: rcp      ! Cloud mixing ratio   [  kg/kg]
   real   , dimension(m1,m2,m3), intent(in)    :: rrp      ! Rain mixing ratio    [  kg/kg]
   real   , dimension(m1,m2,m3), intent(in)    :: rpp      ! Prist. ice mix. rat. [  kg/kg]
   real   , dimension(m1,m2,m3), intent(in)    :: rsp      ! Snow mixing ratio    [  kg/kg]
   real   , dimension(m1,m2,m3), intent(in)    :: rap      ! Aggregates mix. rat. [  kg/kg]
   real   , dimension(m1,m2,m3), intent(in)    :: rgp      ! Graupel mixing ratio [  kg/kg]
   real   , dimension(m1,m2,m3), intent(in)    :: rhp      ! Hail mixing ratio    [  kg/kg]
   real   , dimension(m1,m2,m3), intent(in)    :: q6       ! Graupel int. en.     [   J/kg]
   real   , dimension(m1,m2,m3), intent(in)    :: q7       ! Hail internal en.    [   J/kg]
   !----- Output arguments (some are scratch) ---------------------------------------------!
   real   , dimension(m1)      , intent(inout) :: rliq  ! Liquid mixing ratio     [  kg/kg]
   real   , dimension(m1)      , intent(inout) :: rice  ! Ice mixing ratio        [  kg/kg]
   real   , dimension(m1,m2,m3), intent(inout) :: rv    ! Vapour mixing ratio     [  kg/kg]
   real   , dimension(m1,m2,m3), intent(inout) :: theta ! Potential temperature   [      K]
   !----- Local Variables -----------------------------------------------------------------!
   integer                                     :: i,j,k ! Counter.                [   ----]
   real                                        :: exner ! Exner function          [   J/kg]
   real                                        :: pres  ! Pressure                [     Pa]
   real                                        :: temp  ! Temperature             [      K]
   real                                        :: til   ! Ice-liquid temperature  [      K]
   do j = ja,jz
      do i = ia,iz
         !----- Finding the amount of liquid and ice from the hydrometeor species ---------!
         call integ_liq_ice(m1,availcat,rcp(1:m1,i,j),rrp(1:m1,i,j),rpp(1:m1,i,j)          &
                                       ,rsp(1:m1,i,j),rap(1:m1,i,j),rgp(1:m1,i,j)          &
                                       ,rhp(1:m1,i,j),q6 (1:m1,i,j),q7 (1:m1,i,j)          &
                                       ,rliq(1:m1)   ,rice(1:m1)                 )
         !----- Vapour mixing ratio as the difference b/w total and non-vapour ------------!
         do k = 1,m1
            rv(k,i,j) = rtp(k,i,j) - rliq(k) - rice(k)
         end do

         !----- Finding the potential temperature -----------------------------------------!
         do k = 1,m1
            exner        = pi0(k,i,j) + pp(k,i,j)
            pres         = exner2press(exner)
            !----- Find the first guess. --------------------------------------------------!
            til  = extheta2temp(exner,thp  (k,i,j))
            temp = extheta2temp(exner,theta(k,i,j))
            if (temp > ttripoli) then
               temp = 0.5 * ( til                                                          &
                            + sqrt( til * (til + cpdryi4 * (alvl3*rliq(k)+alvi3*rice(k)))))
            else
               temp = til * (1. + htripolii * (alvl3*rliq(k)+alvi3*rice(k)))
            endif
            !----- First guess for temperature --------------------------------------------!
            temp         = thil2temp(thp(k,i,j),exner,pres,rliq(k),rice(k),temp)
            theta(k,i,j) = extemp2theta(exner,temp)

            if (rv(k,i,j) > rtp(k,i,j) .or. rliq(k) < 0. .or. rice(k) < 0.) then
               write (unit=*,fmt='(a)') '------ MODEL THERMODYNAMIC IS NON-SENSE... ------'
               write (unit=*,fmt='(a,1x,i5,a)'  ) 'In node ',mynum,'...'
               write (unit=*,fmt='(a,1x,i5)'    ) 'I     =',i
               write (unit=*,fmt='(a,1x,i5)'    ) 'J     =',j
               write (unit=*,fmt='(a,1x,i5)'    ) 'K     =',k
               write (unit=*,fmt='(a,1x,es12.5)') 'EXNER =',exner
               write (unit=*,fmt='(a,1x,es12.5)') 'PRESS =',pres
               write (unit=*,fmt='(a,1x,es12.5)') 'THIL  =',thp(k,i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'THETA =',thp(k,i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'TEMP  =',temp
               write (unit=*,fmt='(a,1x,es12.5)') 'RTOT  =',rtp(k,i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'RVAP  =',rv(k,i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'RLIQ  =',rliq(k)
               write (unit=*,fmt='(a,1x,es12.5)') 'RICE  =',rice(k)
               write (unit=*,fmt='(a,1x,es12.5)') 'CLOUD =',rcp(k,i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'RAIN  =',rrp(k,i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'PICE  =',rpp(k,i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'SNOW  =',rsp(k,i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'AGGR  =',rap(k,i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'GRAUP =',rgp(k,i,j)
               write (unit=*,fmt='(a,1x,es12.5)') 'HAIL  =',rhp(k,i,j)
               write (unit=*,fmt='(a)') '-------------------------------------------------'
               call abort_run('Weird thermodynamic state found!','wetthrm3','rthrm.f90')
            end if
         end do

      end do
   end do
   return
end subroutine wetthrm3
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine integ_liq_ice(m1,availcat,rcp,rrp,rpp,rsp,rap,rgp,rhp,q6,q7,rliq,rice)
   !---------------------------------------------------------------------------------------!
   !     This subroutine computes the total liquid water and ice mixing ratios from all    !
   ! micrometeors. It was part of wetthrm3, the only reason it's outside now is that the   !
   ! cumulus parametrisation also needs this part to find the first guess, but not the     !
   ! temperature and potential temperature part.                                           !
   !---------------------------------------------------------------------------------------!
   use therm_lib, only : uint2tl ! ! function
   implicit none
   !----- Input variables -----------------------------------------------------------------!
   integer               , intent(in)    :: m1       ! Vertical dimension         [   ----]
   logical, dimension(*) , intent(in)    :: availcat ! Microphysics handle        [    T|F]
   real   , dimension(m1), intent(in)    :: rcp      ! Cloud mixing ratio         [  kg/kg]
   real   , dimension(m1), intent(in)    :: rrp      ! Rain mixing ratio          [  kg/kg]
   real   , dimension(m1), intent(in)    :: rpp      ! Pristine ice mix. ratio    [  kg/kg]
   real   , dimension(m1), intent(in)    :: rsp      ! Snow mixing ratio          [  kg/kg]
   real   , dimension(m1), intent(in)    :: rap      ! Aggregates mix. ratio      [  kg/kg]
   real   , dimension(m1), intent(in)    :: rgp      ! Graupel mixing ratio       [  kg/kg]
   real   , dimension(m1), intent(in)    :: rhp      ! Hail mixing ratio          [  kg/kg]
   real   , dimension(m1), intent(in)    :: q6       ! Graupel internal energy    [   J/kg]
   real   , dimension(m1), intent(in)    :: q7       ! Hail internal energy       [   J/kg]
   !----- Output variables ----------------------------------------------------------------!
   real   , dimension(m1), intent(out)   :: rliq     ! Liquid mixing ratio        [  kg/kg]
   real   , dimension(m1), intent(out)   :: rice     ! Ice mixing ratio           [  kg/kg]
   !----- Local variables -----------------------------------------------------------------!
   integer                               :: k        ! Counter                    [   ----]
   real                                  :: tcoal    ! Coal. temperature          [      K]
   real                                  :: frcliq   ! Frac. of liquid            [   ----]
   !---------------------------------------------------------------------------------------!



   !----- Initialising total liquid and ice mixing ratios ---------------------------------!
   do k = 1,m1
      rliq(k)  = 0.
      rice(k)  = 0.
   end do
   !---------------------------------------------------------------------------------------!

   !----- Adding cloud droplet mixing ratio to liquid -------------------------------------!
   if (availcat(1)) then
      do k = 1,m1
         rliq(k) = rliq(k) + rcp(k)
      end do
   end if
   !---------------------------------------------------------------------------------------!

   !----- Adding rain mixing ratio to liquid ----------------------------------------------!
   if (availcat(2)) then
      do k = 1,m1
         rliq(k) = rliq(k) + rrp(k)
      end do
   end if
   !---------------------------------------------------------------------------------------!

   !----- Adding pristine ice mixing ratio to ice -----------------------------------------!
   if (availcat(3)) then
      do k = 1,m1
         rice(k) = rice(k) + rpp(k)
      end do
   end if
   !---------------------------------------------------------------------------------------!

   !----- Adding snow mixing ratio to ice -------------------------------------------------!
   if (availcat(4)) then
      do k = 1,m1
         rice(k) = rice(k) + rsp(k)
      end do
   end if
   !---------------------------------------------------------------------------------------!

   !----- Adding aggregates to ice --------------------------------------------------------!
   if (availcat(5)) then
      do k = 1,m1
         rice(k) = rice(k) + rap(k)
      end do
   end if
   !---------------------------------------------------------------------------------------!

   !----- Graupel is mixed, find the liquid fraction and distribute to both phases --------!
   if (availcat(6)) then
      do k = 1,m1
         call uint2tl(q6(k),tcoal,frcliq)
         rliq(k) = rliq(k) + rgp(k) * frcliq
         rice(k) = rice(k) + rgp(k) * (1. - frcliq)
      end do
   end if
   !---------------------------------------------------------------------------------------!

   !----- Hail is also mixed, do the same as graupel --------------------------------------!
   if (availcat(7)) then
      do k = 1,m1
         call uint2tl(q7(k),tcoal,frcliq)
         rliq(k) = rliq(k) + rhp(k) * frcliq
         rice(k) = rice(k) + rhp(k) * (1. - frcliq)
      end do
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine integ_liq_ice
!==========================================================================================!
!==========================================================================================!
