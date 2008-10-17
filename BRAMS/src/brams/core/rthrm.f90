!##################################### Change Log #########################################!
! 5.0.0                                                                                    !
!                                                                                          !
!##########################################################################################!
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved                       !
!  Regional Atmospheric Modeling System - RAMS                                             !
!##########################################################################################!
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

   use mem_grid, only:    &
       ngrid              ! ! intent(in)    - Current grid
   use mem_basic, only:   &
       basic_g            ! ! intent(inout) - Structure with the "basic" variables
   use mem_micro, only:   &
       micro_g            ! ! intent(inout) - Structure containing the hydrometeors
   use mem_scratch, only: &
       scratch,           & ! intent(out) - Scratch structure, for scratch...
       vctr5,             & ! intent(out) - Scratch vector, for scratch...
       vctr6              ! ! intent(out) - Scratch vector, for scratch...
   use therm_lib, only:   &
       level              ! ! intent(in) - Number of H2O phases 
   use micphys, only:     &
       availcat           ! ! intent(in) - Flag: the hydrometeor is available [T|F]
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   integer, intent(in)  :: mzp  ! # of points in Z                                [   ----]
   integer, intent(in)  :: mxp  ! # of points in X                                [   ----]
   integer, intent(in)  :: myp  ! # of points in Y                                [   ----]
   integer, intent(in)  :: ia   ! Node Western edge                               [   ----]
   integer, intent(in)  :: iz   ! Node Eastern edge                               [   ----]
   integer, intent(in)  :: ja   ! Node Southern edge                              [   ----]
   integer, intent(in)  :: jz   ! Node Northern edge                              [   ----]
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !   Calling the appropriate driver according to the microphysics complexity level (aka  !
   ! how many phases H20 will be allowed to exist).                                        !
   !---------------------------------------------------------------------------------------!
   select case (level)
   !----- No condensation, if super-saturation happens, then it will be supersaturated... -!
   case (1)
      call drythrm(mzp,mxp,myp,ia,iz,ja,jz                                                 &
           ,basic_g(ngrid)%thp                     ,basic_g(ngrid)%theta                   &
           ,basic_g(ngrid)%rtp                     ,basic_g(ngrid)%rv                      &
           ,level                                  )

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
      call wetthrm3(mzp,mxp,myp,ia,iz,ja,jz,availcat                                       &
           ,basic_g(ngrid)%pi0                     ,basic_g(ngrid)%pp                      &
           ,basic_g(ngrid)%thp                     ,basic_g(ngrid)%theta                   &
           ,basic_g(ngrid)%rtp                     ,basic_g(ngrid)%rv                      &
           ,micro_g(ngrid)%rcp                     ,micro_g(ngrid)%rrp                     &
           ,micro_g(ngrid)%rpp                     ,micro_g(ngrid)%rsp                     &
           ,micro_g(ngrid)%rap                     ,micro_g(ngrid)%rgp                     &
           ,micro_g(ngrid)%rhp                     ,micro_g(ngrid)%q6                      &
           ,micro_g(ngrid)%q7                      ,vctr5 ,vctr6                           )

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
subroutine drythrm(m1,m2,m3,ia,iz,ja,jz,thil,theta,rt,rv,level)
  !----------------------------------------------------------------------------------------!
  !    This routine calculates theta and rv for the case where no condensate is allowed.   !
  !----------------------------------------------------------------------------------------!

   implicit none
   !----- Arguments: ----------------------------------------------------------------------!
   integer                  , intent(in)    :: m1,m2,m3 ! Grid dimensions          [  ----]
   integer                  , intent(in)    :: ia,iz    ! Lateral boundaries       [  ----]
   integer                  , intent(in)    :: ja,jz    ! Lateral boundaries       [  ----]
   integer                  , intent(in)    :: level    ! # of H2O phases          [  ----]
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
         end do
         if (level == 1) then
            do k = 1,m1
               rv(k,i,j) = rt(k,i,j)
            end do
         end if
      end do
   end do
   return
end subroutine drythrm
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine satadjst(m1,m2,m3,ia,iz,ja,jz,pp,p,thil,theta,t,pi0,rtp,rv,rcp,rvls)
   !---------------------------------------------------------------------------------------!
   !    This routine diagnoses theta, rv, and rcp using a saturation adjustment for the    !
   ! case when water is in the liquid and vapour phases only.                              !
   !---------------------------------------------------------------------------------------!
   use rconstants, only: &
           cpi           & ! intent(in)
          ,p00           & ! intent(in)
          ,alvl          & ! intent(in)
          ,cp            & ! intent(in)
          ,cpor          ! ! intent(in)
   !----- External functions --------------------------------------------------------------!
   use therm_lib , only: &
           rslf    ! Sat. mixing ratio function
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
            p(k,i,j) = p00 * (cpi* exner) ** cpor
            !----- First guess for temperature and liquid mixing ratio --------------------!
            t(k,i,j)   = cpi * thil(k,i,j) * exner
            rcp(k,i,j) = max(0.,rtp(k,i,j) - rslf(p(k,i,j),t(k,i,j)))
            !----- Adjusting the accordingly to the saturation point ----------------------!
            call thil2tqliq(thil(k,i,j),exner,p(k,i,j),rtp(k,i,j),rcp(k,i,j),t(k,i,j)      &
                           ,rv(k,i,j),rvls(k,i,j))
            theta(k,i,j) = cp * t(k,i,j) / exner
         end do
      end do
   end do

   return
end subroutine satadjst
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine wetthrm3(m1,m2,m3,ia,iz,ja,jz,availcat,pi0,pp,thp,theta,rtp,rv,rcp,rrp,rpp,rsp  &
                   ,rap,rgp,rhp,q6,q7,rliq,rice)

   !---------------------------------------------------------------------------------------!
   !    This routine calculates theta and rv for "level 3 microphysics"
   ! given prognosed theta_il, cloud, rain, pristine ice, snow, aggregates,
   ! graupel, hail, q6, and q7.

   use rconstants, only: &
        cpi,      & ! intent(in)
        cp,       & ! intent(in)
        cpor,     & ! intent(in)
        p00,      & ! intent(in)
        ttripoli, & ! intent(in)
        alvl,     & ! intent(in)
        alvi,     & ! intent(in)
        cpi4,     & ! intent(in)
        htripolii ! ! intent(in)       
   !---- External function ----------------------------------------------------------------!
   use therm_lib , only: &
        thil2temp ! ! Theta_il => Temperature function

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
            pres         = p00 * (cpi * exner)**cpor
            !----- Finding the first guess ------------------------------------------------!
            til  = cpi * thp(k,i,j)   * exner
            temp = cpi * theta(k,i,j) * exner
            if (temp > ttripoli) then
               temp = 0.5 * (til + sqrt(til * (til + cpi4*(alvl*rliq(k)+alvi*rice(k)))))
            else
               temp = til * (1. + htripolii * (alvl*rliq(k)+alvi*rice(k)))
            endif
            !----- First guess for temperature --------------------------------------------!
            temp         = thil2temp(thp(k,i,j),exner,pres,rliq(k),rice(k),temp)
            theta(k,i,j) = cp * temp / exner
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
   use therm_lib, only: qtk
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
         call qtk(q6(k),tcoal,frcliq)
         rliq(k) = rliq(k) + rgp(k) * frcliq
         rice(k) = rice(k) + rgp(k) * (1. - frcliq)
      end do
   end if
   !---------------------------------------------------------------------------------------!

   !----- Hail is also mixed, do the same as graupel --------------------------------------!
   if (availcat(7)) then
      do k = 1,m1
         call qtk(q7(k),tcoal,frcliq)
         rliq(k) = rliq(k) + rhp(k) * frcliq
         rice(k) = rice(k) + rhp(k) * (1. - frcliq)
      end do
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine integ_liq_ice
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes a consistent set of temperature and condensated phases mix- !
! ing ratio for a given theta_il, Exner function, and total mixing ratio. This is very     !
! similar to the function thil2temp, except that now we don't know rliq and rice, and for  !
! this reason they also become functions of temperature, since they are defined as         !
! rtot-rsat(T,p), remembering that rtot and p are known. If the air is not saturated, we   !
! rather use the fact that theta_il = theta and skip the hassle. Otherwise, we use iter-   !
! ative methods. We will always try Newton's method, since it converges fast. The caveat   !
! is that Newton may fail, and it actually does fail very close to the triple point,       !
! because the saturation vapour pressure function has a "kink" at the triple point         !
! (continuous, but not differentiable). If that's the case, then we fall back to a modifi- !
! ed regula falsi (Illinois) method, which is a mix of secant and bisection and will       !
! converge.                                                                                !
!------------------------------------------------------------------------------------------!
subroutine thil2tqall(thil,exner,pres,rtot,rliq,rice,temp,rvap,rsat)
   use rconstants, only : alvl,alvi,allii,cp,cpi,t00,toodry,t3ple,ttripoli
   !----- External functions --------------------------------------------------------------!
   use therm_lib , only : &
          rslif           & ! Function to compute sat. mixing ratio.
         ,toler           & ! Tolerance 
         ,theta_iceliq    & ! Function to compute theta_il
         ,dthetail_dt     & ! d(Theta_il)/dT
         ,maxfpo          ! ! Maximum # of iterations before giving up.

   implicit none
   
   real, intent(in)    :: thil        ! Ice-liquid water potential temperature     [     K]
   real, intent(in)    :: exner       ! Exner function                             [J/kg/K]
   real, intent(in)    :: pres        ! Pressure                                   [    Pa]
   real, intent(in)    :: rtot        ! Total mixing ratio                         [ kg/kg]
   real, intent(out)   :: rliq        ! Liquid water mixing ratio                  [ kg/kg]
   real, intent(out)   :: rice        ! Ice mixing ratio                           [ kg/kg]
   real, intent(inout) :: temp        ! Temperature                                [     K]
   real, intent(out)   :: rvap        ! Water vapour mixing ratio                  [ kg/kg]
   real, intent(out)   :: rsat        ! Saturation water vapour mixing ratio       [ kg/kg]
   !----- Local variables -----------------------------------------------------------------!
   real                :: tempa,tempz ! Aux. vars for regula falsi iteration
   real                :: funa,funz   ! The functions for regula falsi
   real                :: funnow      ! Function at this iteration.
   real                :: delta       ! Aux. var in case we need regula falsi.
   real                :: deriv       ! Derivative of this function.
   integer             :: itn,itb,ii  ! Iteration counter
   logical             :: converged   ! Convergence handle
   logical             :: zside       ! Aux. Flag, for two purposes:
                                      ! 1. Found a good 2nd guess for regula falsi.
                                      ! 2. I retained the "zside" (T/F)

   !---------------------------------------------------------------------------------------!
   !      First check: try to find temperature assuming sub-saturation and check if this   !
   ! is the case. If it is, then there is no need to go through the iterative loop.        !
   !---------------------------------------------------------------------------------------!
   tempz  = cpi * thil * exner
   rsat   = max(toodry,rslif(pres,tempz))
   if (tempz >= t3ple) then
      rliq = max(0.,rtot-rsat)
      rice = 0.
   else
      rice = max(0.,rtot-rsat)
      rliq = 0.
   end if
   rvap = rtot-rliq-rice

   !---------------------------------------------------------------------------------------!
   !    If rtot < rsat, this is not saturated, we can leave the subroutine and bypass the  !
   ! iterative part.                                                                       !
   !---------------------------------------------------------------------------------------!
   if (rtot < rsat) then
      temp = tempz
      return
   end if

   !---------------------------------------------------------------------------------------!
   !   If not, then use the temperature the user gave as first guess and solve iterative-  !
   ! ly. We use the user instead of what we just found because if the air is saturated,    !
   ! then this can be too far off which may be bad for Newton's method.                    !
   !---------------------------------------------------------------------------------------!
   tempz = temp
   rsat   = max(toodry,rslif(pres,tempz))
   if (tempz >= t3ple) then
      rliq = max(0.,rtot-rsat)
      rice = 0.
   else
      rice = max(0.,rtot-rsat)
      rliq = 0.
   end if
   rvap = rtot-rliq-rice


   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   !write (unit=46,fmt='(a)') '-------------------------------------------------------------'
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!


   !---------------------------------------------------------------------------------------!
   !     Finding the function. We are seeking a temperature which is associated with the   !
   ! theta_il we provided. Thus, the function is simply the difference between the         !
   ! theta_il associated with our guess and the actual theta_il.                           !
   !---------------------------------------------------------------------------------------!
   funnow = theta_iceliq(exner,tempz,rliq,rice)
   !----- Updating the derivative. --------------------------------------------------------!
   deriv  = dthetail_dt(.false.,funnow,exner,pres,tempz,rliq,rice)
   funnow = funnow - thil

   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   !write (unit=46,fmt='(a,1x,i5,1x,6(a,1x,f11.4,1x),a,1x,es11.4,1x)')                      &
   !   'NEWTON: it=',0,'temp=',tempz-t00,'rsat=',1000.*rsat,'rliq=',1000.*rliq              &
   !  ,'rice=',1000.*rice,'rvap=',1000.*rvap,'fun=',funnow,'deriv=',deriv
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

   !---------------------------------------------------------------------------------------!
   !    Now we enter at the Newton's method iterative loop. We are always going to try     !
   ! this first, because it's fast, but if it turns out to be a dangerous choice or if it  !
   ! doesn't converge fast, we will fall back to regula falsi.                             !
   !    We start by initialising the flag and copying temp to tempz, the newest guess.     !
   !---------------------------------------------------------------------------------------!
   converged=.false.
   newloop: do itn=1,maxfpo/6
      !------------------------------------------------------------------------------------!
      !     Saving previous guess. We also save the function is in case we withdraw        !
      ! Newton's and switch to regula falsi.                                               !
      !------------------------------------------------------------------------------------!
      funa  = funnow
      tempa = tempz

      !----- Go to bisection if the derivative is too flat (too dangerous...) -------------!
      if (abs(deriv) < toler) exit newloop

      tempz = tempa - funnow / deriv

      !----- Finding the mixing ratios associated with this guess -------------------------!
      rsat  = max(toodry,rslif(pres,tempz))
      if (tempz >= t3ple) then
         rliq = max(0.,rtot-rsat)
         rice = 0.
      else
         rice = max(0.,rtot-rsat)
         rliq = 0.
      end if
      rvap = rtot-rliq-rice

      !----- Updating the function --------------------------------------------------------!
      funnow = theta_iceliq(exner,tempz,rliq,rice)
      !----- Updating the derivative. -----------------------------------------------------!
      deriv  = dthetail_dt(.false.,funnow,exner,pres,tempz,rliq,rice)
      funnow = funnow - thil

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=46,fmt='(a,1x,i5,1x,6(a,1x,f11.4,1x),a,1x,es11.4,1x)')                   &
      !   'NEWTON: it=',itn,'temp=',tempz-t00,'rsat=',1000.*rsat,'rliq=',1000.*rliq         &
      !  ,'rice=',1000.*rice,'rvap=',1000.*rvap,'fun=',funnow,'deriv=',deriv
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      
      converged = abs(tempa-tempz) < toler*tempz
      !------------------------------------------------------------------------------------!
      !   Convergence. The temperature will be the mid-point between tempa and tempz. Fix  !
      ! the mixing ratios and return.                                                      !
      !------------------------------------------------------------------------------------!
      if (converged) then
         temp = 0.5 * (tempa+tempz)
         rsat  = max(toodry,rslif(pres,temp))
         if (temp >= t3ple) then
            rliq = max(0.,rtot-rsat)
            rice = 0.
         else
            rice = max(0.,rtot-rsat)
            rliq = 0.
         end if
         rvap = rtot-rliq-rice
         exit newloop
      end if

   end do newloop
   !---------------------------------------------------------------------------------------!

   !----- For debugging only --------------------------------------------------------------!
   itb = itn+1

   if (.not. converged) then
      !------------------------------------------------------------------------------------!
      !    If I reach this point, then it means that Newton's method failed finding the    !
      ! equilibrium, so we are going to use the regula falsi instead. If Newton's method   !
      ! didn't  converge, we use tempa as one guess and now we seek a tempz with opposite  !
      ! sign.                                                                              !
      !------------------------------------------------------------------------------------!
      !----- Check funa and funnow have opposite signs. If so, we are ready to go ---------!
      if (funa*funnow < 0) then
         funz = funnow
         zside = .true.
      !------------------------------------------------------------------------------------!
      !    Looking for a guess. Extrapolate funa linearly, trying to get the -funa. We     !
      ! don't need it to be funa, just with the opposite sign. If that's not enough, we    !
      ! keep going further...                                                              !
      !------------------------------------------------------------------------------------!
      else
         if (abs(funnow-funa) < toler*tempa) then
            delta = 100.*toler*tempa
         else
            delta = max(abs(funa*(tempz-tempa)/(funnow-funa)),100.*toler*tempa)
         end if
         tempz = tempa + delta
         funz  = funa
         !----- Just to enter at least once. The 1st time tempz=tempa-2*delta -------------!
         zside = .false. 
         zgssloop: do itb=1,maxfpo
             tempz = tempa + real((-1)**itb * (itb+3)/2) * delta
             rsat   = max(toodry,rslif(pres,tempz))
             if (tempz >= t3ple) then
                rliq = max(0.,rtot-rsat)
                rice = 0.
             else
                rice = max(0.,rtot-rsat)
                rliq = 0.
             end if
             rvap = rtot-rliq-rice
             funz = theta_iceliq(exner,tempz,rliq,rice) - thil
             zside = funa*funz < 0
             if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside) then
            write (unit=*,fmt='(a)') ' No second guess for you...'
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempa=',tempa,'funa=',funa
            write (unit=*,fmt='(2(a,1x,es14.7))') 'tempz=',tempz,'func=',funz
            write (unit=*,fmt='(1(a,1x,es14.7))') 'delta=',delta
            call abort_run('Failed finding the second guess for regula falsi'              &
                          ,'thil2tqall','rthrm.f90')
         end if
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Now we loop until convergence is achieved. One important thing to notice is    !
      ! that Newton's method fail only when T is almost T3ple, which means that ice and    !
      ! liquid should be present, and we are trying to find the saturation point with all  !
      ! ice or all liquid. This will converge but the final answer will contain signifi-   !
      ! cant error. To reduce it we redistribute the condensates between ice and liquid    !
      ! conserving the total condensed mixing ratio.                                       !
      !------------------------------------------------------------------------------------!
      fpoloop: do itb=itn,maxfpo
         temp = (funz*tempa-funa*tempz)/(funz-funa)
         !temp = 0.5*(tempa+tempz)
         !----- Distributing vapour into the three phases ---------------------------------!
         rsat   = max(toodry,rslif(pres,temp))
         rvap   = min(rtot,rsat)
         if (temp >= t3ple) then
            rliq = max(0.,rtot-rsat)
            rice = 0.
         else
            rliq = 0.
            rice = max(0.,rtot-rsat)
         end if
         !----- Updating function ---------------------------------------------------------!
         funnow = theta_iceliq(exner,temp,rliq,rice) - thil

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write (unit=46,fmt='(a,1x,i5,1x,10(a,1x,f11.4,1x))')                              &
         !   'REGFAL: it=',itb,'temp=',temp-t00,'tempa=',tempa-t00,'tempz=',tempz-t00       &
         !  ,'rsat=',1000.*rsat,'rliq=',1000.*rliq,'rice=',1000.*rice,'rvap=',1000.*rvap    &
         !  ,'fun=',funnow,'funa=',funa,'funz=',funz
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         !---------------------------------------------------------------------------------!
         !    Checking for convergence or lucky guess. If it did, return, we found the     !
         ! solution. Otherwise, constrain the guesses.                                     !
         !---------------------------------------------------------------------------------!
         converged = abs(temp-tempa) < toler*temp .and. abs(temp-tempz) < toler*temp 
         if (funnow == 0. .or. converged) then
            converged = .true.
            exit fpoloop
         elseif (funnow*funa < 0.) then 
            tempz = temp
            funz  = funnow
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa=funa * 0.5
            !----- We just updated zside, setting zside to true. --------------------------!
            zside = .true.
         else
            tempa = temp
            funa  = funnow
            !----- If we are updating aside again, modify zside (Illinois method) ---------!
            if (.not.zside) funz = funz * 0.5
            !----- We just updated aside, setting zside to false --------------------------!
            zside = .false. 
         end if

      end do fpoloop

      !------------------------------------------------------------------------------------!
      !    Almost done... Usually when the method goes through regula falsi, it means that !
      ! the temperature is too close to the triple point, and often all three phases will  !
      ! coexist. The problem with the method is that it converges for temperature, but     !
      ! whenever regula falsi is called the function evaluation is usually far from zero.  !
      ! This can be improved by finding a better partition between ice and liquid given    !
      ! the temperature and saturation mixing ratio we just found. So just to round these  !
      ! edges, we will invert the ice-liquid potential temperature using the set of tem-   !
      ! perature and rsat, and fiding the liquid mixing ratio.                             !
      !------------------------------------------------------------------------------------!
      if (abs(temp-t3ple) < toler*temp) then
         rliq = min(rtot-rsat,max(0.,                                                      &
                    allii*(alvi*(rtot-rsat)+cp*max(temp,ttripoli)                          &
                         *log(cpi*exner*thil/temp))))
         rice = max(0.,rtot-rsat-rliq)
         funnow = theta_iceliq(exner,temp,rliq,rice) - thil
      end if

      itb=itb+1
   end if

   if (.not. converged) then
      write (unit=*,fmt='(60a1)')        ('-',ii=1,60)
      write (unit=*,fmt='(a)')           ' THIL2TQALL failed!'
      write (unit=*,fmt='(a)')           ' '
      write (unit=*,fmt='(a)')           ' -> Input: '
      write (unit=*,fmt='(a,1x,f12.5)')  '    THETA_IL [     K]:',thil
      write (unit=*,fmt='(a,1x,f12.5)')  '    EXNER    [J/kg/K]:',exner
      write (unit=*,fmt='(a,1x,f12.5)')  '    RTOT     [  g/kg]:',1000.*rtot
      write (unit=*,fmt='(a)')           ' '
      write (unit=*,fmt='(a)')           ' -> Output: '
      write (unit=*,fmt='(a,1x,i12)')    '    ITERATIONS       :',itb
      write (unit=*,fmt='(a,1x,f12.5)')  '    TEMP     [    °C]:',temp-t00
      write (unit=*,fmt='(a,1x,f12.5)')  '    RVAP     [  g/kg]:',1000.*rvap
      write (unit=*,fmt='(a,1x,f12.5)')  '    RLIQ     [  g/kg]:',1000.*rliq
      write (unit=*,fmt='(a,1x,f12.5)')  '    RICE     [  g/kg]:',1000.*rice
      write (unit=*,fmt='(a,1x,f12.5)')  '    TEMPA    [    °C]:',tempa-t00
      write (unit=*,fmt='(a,1x,f12.5)')  '    TEMPZ    [    °C]:',tempz-t00
      write (unit=*,fmt='(a,1x,es12.5)') '    FUNA     [     K]:',funnow
      write (unit=*,fmt='(a,1x,es12.5)') '    FUNZ     [     K]:',funnow
      write (unit=*,fmt='(a,1x,es12.5)') '    DERIV    [   ---]:',deriv
      write (unit=*,fmt='(a,1x,es12.5)') '    ERR_A    [   ---]:',abs(temp-tempa)/temp
      write (unit=*,fmt='(a,1x,es12.5)') '    ERR_Z    [   ---]:',abs(temp-tempz)/temp
      write (unit=*,fmt='(a)')           ' '
      write (unit=*,fmt='(60a1)')        ('-',ii=1,60)
      call abort_run('Failed finding equilibrium, I gave up!','thil2tqall','rthrm.f90')
   end if
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   !write (unit=46,fmt='(a,1x,i5,1x,6(a,1x,f11.4,1x))')                                     &
   !   'ANSWER: it=',itb,'funf=',funnow,'temp=',temp-t00                                    &
   !  ,'rsat=',1000.*rsat,'rliq=',1000.*rliq,'rice=',1000.*rice,'rvap=',1000.*rvap
   !write (unit=46,fmt='(a)') '-------------------------------------------------------------'
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   return
end subroutine thil2tqall
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes a consistent set of temperature and condensated phases mix- !
! ing ratio for a given theta_il, Exner function, and total mixing ratio. This is very     !
! similar to the function thil2temp, except that now we don't know rliq. Rliq becomes      !
! function of temperature, since it is defined as rtot-rsat(T,p), remembering that rtot    !
! and p are known. If the air is not saturated, we rather use the fact that theta_il =     !
! theta and skip the hassle. Otherwise, we use iterative methods. We will always try       !
! Newton's method, since it converges fast. Not always will Newton converge, and if that's !
! the case we use a modified regula falsi (Illinois) method. This method is a mix of sec-  !
! ant and bisection and will always converge.                                              !
!------------------------------------------------------------------------------------------!
subroutine thil2tqliq(thil,exner,pres,rtot,rliq,temp,rvap,rsat)
   use rconstants, only : alvl,cp,cpi,toodry,ttripoli

   !----- External functions --------------------------------------------------------------!
   use therm_lib , only : &
          rslf            & ! Function to compute sat. mixing ratio.
         ,theta_iceliq    & ! Function to compute theta_il for a given state
         ,dthetail_dt     & ! Function to compute d(theta_il)/dT
         ,toler           & ! Tolerance 
         ,maxfpo          ! ! Maximum # of iterations before giving up.
   implicit none
   
   real, intent(in)    :: thil        ! Ice-liquid water potential temperature     [     K]
   real, intent(in)    :: exner       ! Exner function                             [J/kg/K]
   real, intent(in)    :: pres        ! Pressure                                   [    Pa]
   real, intent(in)    :: rtot        ! Total mixing ratio                         [ kg/kg]
   real, intent(out)   :: rliq        ! Liquid water mixing ratio                  [ kg/kg]
   real, intent(inout) :: temp        ! Temperature                                [     K]
   real, intent(out)   :: rvap        ! Water vapour mixing ratio                  [ kg/kg]
   real, intent(out)   :: rsat        ! Saturation water vapour mixing ratio       [ kg/kg]
   !----- Local variables -----------------------------------------------------------------!
   real                :: tempa,tempz ! Aux. vars for regula falsi iteration
   real                :: funa,funz   ! The functions for regula falsi
   real                :: funnow      ! Function at this iteration.
   real                :: delta       ! Aux. var in case we need regula falsi.
   real                :: deriv       ! Derivative of this function.
   integer             :: itn,itb     ! Iteration counter
   logical             :: converged   ! Convergence handle
   logical             :: zside       ! Aux. Flag, for two purposes:
                                      ! 1. Found a good 2nd guess for regula falsi.
                                      ! 2. I retained the "zside" (T/F)

   !---------------------------------------------------------------------------------------!
   !      First check: try to find temperature assuming sub-saturation and check if this   !
   ! is the case. If it is, then there is no need to go through the iterative loop.        !
   !---------------------------------------------------------------------------------------!
   tempz = cpi * thil * exner
   rsat  = max(toodry,rslf(pres,tempz))
   rliq  = max(0.,rtot-rsat)
   rvap  = rtot-rliq

   !---------------------------------------------------------------------------------------!
   !    If rtot < rsat, this is not saturated, we can leave the subroutine and bypass the  !
   ! iterative part.                                                                       !
   !---------------------------------------------------------------------------------------!
   if (rtot < rsat) then
      temp = tempz
      return
   end if

   !---------------------------------------------------------------------------------------!
   !   If not, then use the temperature the user gave as first guess and solve iterative-  !
   ! ly. We use the user instead of what we just found because if the air is saturated,    !
   ! then this can be too far off which may be bad for Newton's method.                    !
   !---------------------------------------------------------------------------------------!
   tempz = temp
   rsat   = max(toodry,rslf(pres,tempz))
   rliq = max(0.,rtot-rsat)
   rvap = rtot-rliq


   !---------------------------------------------------------------------------------------!
   !     Finding the function and its derivative. We are seeking a temperature which is    !
   ! associated with the theta_il we provided. Thus, the function is simply the difference !
   ! between the theta_il associated with our guess and the actual theta_il.               !
   !     To find the derivative, we use the fact that rliq = rtot - rsat(T,p). When        !
   ! T < T(Tripoli), then the temperature at the denominator becomes constant so the       !
   ! derivative becomes different.                                                         !
   !---------------------------------------------------------------------------------------!
   funnow = theta_iceliq(exner,tempz,rliq,0.) ! Finding thil from our guess
   deriv  = dthetail_dt(.false.,funnow,exner,pres,tempz,rliq)
   funnow = funnow - thil ! Computing the function


   !---------------------------------------------------------------------------------------!
   !    Now we enter at the Newton's method iterative loop. We are always going to try     !
   ! this first, because it's fast, but if it turns out to be a dangerous choice or if it  !
   ! doesn't converge fast, we will fall back to regula falsi.                             !
   !    We start by initialising the flag and copying temp to tempz, the newest guess.     !
   !---------------------------------------------------------------------------------------!
   converged=.false.
   newloop: do itn=1,maxfpo/6
      !------------------------------------------------------------------------------------!
      !     Saving previous guess. We also save the function is in case we withdraw        !
      ! Newton's and switch to regula falsi.                                               !
      !------------------------------------------------------------------------------------!
      funa  = funnow
      tempa = tempz

      !----- Go to bisection if the derivative is too flat (too dangerous...) -------------!
      if (abs(deriv) < toler) exit newloop

      tempz = tempa - funnow / deriv

      !----- Finding the mixing ratios associated with this guess -------------------------!
      rsat  = max(toodry,rslf(pres,tempz))
      rliq = max(0.,rtot-rsat)
      rvap = rtot-rliq

      !----- Updating the function and its derivative -------------------------------------!
      funnow = theta_iceliq(exner,tempz,rliq,0.)
      deriv = dthetail_dt(.false.,funnow,exner,pres,tempz,rliq)
      funnow = funnow - thil

      converged = abs(tempa-tempz) < toler*tempz
      !------------------------------------------------------------------------------------!
      !   Convergence. The temperature will be the mid-point between tempa and tempz. Fix  !
      ! the mixing ratios and return.                                                      !
      !------------------------------------------------------------------------------------!
      if (converged) then
         temp = 0.5 * (tempa+tempz)
         rsat  = max(toodry,rslf(pres,temp))
         rliq = max(0.,rtot-rsat)
         rvap = rtot-rliq
         exit newloop
      end if

   end do newloop
   !---------------------------------------------------------------------------------------!


   if (.not. converged) then
      !------------------------------------------------------------------------------------!
      !    If I reach this point, then it means that Newton's method failed finding the    !
      ! equilibrium, so we are going to use the regula falsi instead. If Newton's method   !
      ! didn't  converge, we use tempa as one guess and now we seek a tempz with opposite  !
      ! sign.                                                                              !
      !------------------------------------------------------------------------------------!
      !----- Check funa and funnow have opposite signs. If so, we are ready to go ---------!
      if (funa*funnow < 0) then
         funz = funnow
         zside = .true.
      !------------------------------------------------------------------------------------!
      !    Looking for a guess. Extrapolate funa linearly, trying to get the -funa. We     !
      ! don't need it to be funa, just with the opposite sign. If that's not enough, we    !
      ! keep going further...                                                              !
      !------------------------------------------------------------------------------------!
      else
         if (abs(funnow-funa) < toler*tempa) then
            delta = 100.*toler*tempa
         else
            delta = max(abs(funa*(tempz-tempa)/(funnow-funa)),100.*toler*tempa)
         end if
         tempz = tempa + delta
         funz  = funa
         !----- Just to enter at least once. The 1st time tempz=tempa-2*delta -------------!
         zside = .false. 
         zgssloop: do itb=1,maxfpo
            tempz = tempz + real((-1)**itb * (itb+3)/2) * delta
            rsat  = max(toodry,rslf(pres,tempz))
            rliq  = max(0.,rtot-rsat)
            rvap  = rtot-rliq
            funz  = theta_iceliq(exner,tempz,rliq,0.) - thil
            zside = funa*funz < 0
            if (zside) exit zgssloop
         end do zgssloop
         if (.not. zside)                                                                  &
            call abort_run('Failed finding the second guess for regula falsi'              &
                          ,'thil2tqliq','rthrm.f90')
      end if
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !     Now we loop until convergence is achieved.                                     !
      !------------------------------------------------------------------------------------!
      fpoloop: do itb=itn,maxfpo
         temp = (funz*tempa-funa*tempz)/(funz-funa)
         !----- Distributing vapour into the three phases ---------------------------------!
         rsat   = max(toodry,rslf(pres,temp))
         rvap   = min(rtot,rsat)
         rliq   = max(0.,rtot-rsat)
         !----- Updating function ---------------------------------------------------------!
         funnow = theta_iceliq(exner,tempz,rliq,0.) - thil

         !---------------------------------------------------------------------------------!
         !    Checking for convergence or lucky guess. If it did, return, we found the     !
         ! solution. Otherwise, constrain the guesses.                                                    !
         !---------------------------------------------------------------------------------!
         converged = abs(temp-tempa)< toler*temp  .and. abs(temp-tempz) < toler*temp
         if (funnow == 0. .or. converged) then
            converged = .true.
            exit fpoloop
         elseif (funnow*funa < 0.) then 
            tempz = temp
            funz  = funnow
            !----- If we are updating zside again, modify aside (Illinois method) ---------!
            if (zside) funa=funa * 0.5
            !----- We just updated zside, setting zside to true. --------------------------!
            zside = .true.
         else
            tempa = temp
            funa  = funnow
            !----- If we are updating aside again, modify zside (Illinois method) ---------!
            if (.not.zside) funz = funz * 0.5
            !----- We just updated aside, setting zside to false --------------------------!
            zside = .false. 
         end if

      end do fpoloop

   end if

   if (.not. converged) call abort_run('Failed finding equilibrium, I gave up!'            &
                                      ,'thil2tqliq','rthrm.f90')
   return
end subroutine thil2tqliq
!==========================================================================================!
!==========================================================================================!
