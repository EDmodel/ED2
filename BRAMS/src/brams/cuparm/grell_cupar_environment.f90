!==========================================================================================!
! grell_cupar_environment.f90                                                              !
!                                                                                          !
!    This file contains subroutines that will calculate environment related stuff.         !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the values of some thermodynamic variables at the cloud level !
!------------------------------------------------------------------------------------------!
subroutine grell_thermo_cldlev(mkx,mgmzp,z_cup,exner,thil,t,qtot,qliq,qice,co2,exnersur    &
                              ,thilsur,tsur,qtotsur,qliqsur,qicesur,co2sur,exner_cup,p_cup &
                              ,t_cup,thil_cup,qtot_cup,qvap_cup,qliq_cup,qice_cup,qsat_cup &
                              ,co2_cup,rho_cup,theiv_cup,theivs_cup)

   use rconstants, only : p00,cpi,cpor,t3ple
   use therm_lib , only : & 
           thetaeiv       & ! Function that computes thetae_iv
          ,thetaeivs      & ! Function that computes sat. thetae_iv
          ,idealdens      & ! Function that computes density for ideal gas
          ,rslif          & ! Function that computes saturation mixing ratio
          ,theta_iceliq   ! ! Function that computes theta_il

   implicit none
   !------ Input variables ----------------------------------------------------------------!
   integer                  , intent(in)  :: mkx      ! # of Grell levels
   integer                  , intent(in)  :: mgmzp    ! # of levels for dimensions
   real   , dimension(mgmzp), intent(in)  :: z_cup    ! Height @ cloud levels      [     m]
   !------ Input variables, at model levels -----------------------------------------------!
   real   , dimension(mgmzp), intent(in)  :: exner    ! Exner function             [  J/kg]
   real   , dimension(mgmzp), intent(in)  :: thil     ! Theta_il                   [     K]
   real   , dimension(mgmzp), intent(in)  :: t        ! Temperature                [     K]
   real   , dimension(mgmzp), intent(in)  :: qtot     ! Total mixing ratio         [ kg/kg]
   real   , dimension(mgmzp), intent(in)  :: qliq     ! Liquid water mixing ratio  [ kg/kg]
   real   , dimension(mgmzp), intent(in)  :: qice     ! Ice mixing ratio           [ kg/kg]
   real   , dimension(mgmzp), intent(in)  :: co2      ! CO2 mixing ratio           [   ppm]
   !------ Input variables, at surface ----------------------------------------------------!
   real                     , intent(in)  :: exnersur ! Sfc. Exner function        [  J/kg]
   real                     , intent(in)  :: thilsur  ! Sfc. theta_il              [     K]
   real                     , intent(in)  :: tsur     ! Sfc. temperature           [     K]
   real                     , intent(in)  :: qtotsur  ! Sfc. total mixing ratio    [ kg/kg]
   real                     , intent(in)  :: qliqsur  ! Sfc. water mixing ratio    [ kg/kg]
   real                     , intent(in)  :: qicesur  ! Sfc. ice mixing ratio      [ kg/kg]
   real                     , intent(in)  :: co2sur   ! Sfc. CO2 mixing ratio      [   ppm]
   !------ Output variables, at cloud levels ----------------------------------------------!
   real   , dimension(mgmzp), intent(inout) :: exner_cup  ! Exner function         [J/kg/K]
   real   , dimension(mgmzp), intent(inout) :: p_cup      ! Pressure               [    Pa]
   real   , dimension(mgmzp), intent(inout) :: t_cup      ! Temperature            [     K]
   real   , dimension(mgmzp), intent(inout) :: thil_cup   ! Theta_il               [     K]
   real   , dimension(mgmzp), intent(inout) :: qtot_cup   ! Total mixing ratio     [ kg/kg]
   real   , dimension(mgmzp), intent(inout) :: qvap_cup   ! Vapour mixing ratio    [ kg/kg]
   real   , dimension(mgmzp), intent(inout) :: qliq_cup   ! Liq. water mix. ratio  [ kg/kg]
   real   , dimension(mgmzp), intent(inout) :: qice_cup   ! Ice mixing ratio       [ kg/kg]
   real   , dimension(mgmzp), intent(inout) :: qsat_cup   ! Sat. vapour mix. ratio [ kg/kg]
   real   , dimension(mgmzp), intent(inout) :: co2_cup    ! CO2 mixing ratio       [   ppm]
   real   , dimension(mgmzp), intent(inout) :: rho_cup    ! Air density            [ kg/m³]
   real   , dimension(mgmzp), intent(inout) :: theiv_cup  ! Thetae_iv              [     K]
   real   , dimension(mgmzp), intent(inout) :: theivs_cup ! Sat. thetae_iv         [     K]
   !------ Local variables ----------------------------------------------------------------!
   integer           :: k        ! Level counter
   real              :: qtotsat  ! Total mixing ratio if air was saturated.        [ kg/kg]
   real              :: dummy    ! Dummy variable for unwanted output. 
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Here I will interpolate Exner function, theta_il and total mixing ratio, because   !
   ! these are the conserved variables. Liquid and ice mixing ratios are interpolated here !
   ! but they are just the first guess, they may change later on the code.                 !
   !---------------------------------------------------------------------------------------!
   do k = 2,mkx
      exner_cup(k) = sqrt(exner(k-1) * exner(k)) !----- Log-interpolation -----------------!
      qtot_cup(k)  = sqrt(qtot(k-1)  * qtot(k) ) !----- Log-interpolation -----------------!
      co2_cup(k)   = sqrt(co2(k-1)   * co2(k)  ) !----- Log-interpolation -----------------!
      t_cup   (k)  = 0.5 * (t(k-1)   + t(k)    )
   end do
   !----- 1st. level, boundary condition only, just copying the surface variables. --------!
   exner_cup(1) = exnersur
   qtot_cup (1) = qtotsur
   co2_cup  (1) = co2sur
   t_cup    (1) = tsur

   !---------------------------------------------------------------------------------------!
   !    Now I'll compute the other variables based on the interpolated ones                !
   !---------------------------------------------------------------------------------------!
   do k = 1,mkx
      !------ Pressure, straightforward and it could be interpolated too ------------------!
      p_cup(k)   = p00*(cpi*exner_cup(k))**cpor

      !------ Finding liquid and ice mixing ratio -----------------------------------------!
      qsat_cup(k) = rslif(p_cup(k),t_cup(k))
      if (t_cup(k) >= t3ple) then
          qliq_cup(k) = max(0.,qtot_cup(k)-qsat_cup(k))
          qice_cup(k) = 0.
      else
          qliq_cup(k) = 0.
          qice_cup(k) = max(0.,qtot_cup(k)-qsat_cup(k))
      end if
      qvap_cup(k) = qtot_cup(k)-qliq_cup(k)-qice_cup(k)
      
      !------ Finding the new ice-liquid potential temperature ----------------------------!
      thil_cup(k)   = theta_iceliq(exner_cup(k),t_cup(k),qliq_cup(k),qice_cup(k))
      !------ Finding the air density -----------------------------------------------------!
      rho_cup(k)    = idealdens(p_cup(k),t_cup(k),qvap_cup(k),qtot_cup(k)) 
      !------ Finding the ice-vapour equivalent potential temperature ---------------------!
      theiv_cup(k)  = thetaeiv(thil_cup(k),p_cup(k),t_cup(k),qvap_cup(k),qtot_cup(k))
      !------ Finding the saturation ice-vapour equivalent potential temperature ----------!
      theivs_cup(k) = thetaeivs(thil_cup(k),t_cup(k),qsat_cup(k),qliq_cup(k),qice_cup(k))
   end do



   return
end subroutine grell_thermo_cldlev
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes the buoyancy acceleration at a given level, given the        !
! environment and draft properties. We use the definition of buoyancy based on the density !
! difference, which is equivalent to many other ways to compute buoyancy, like:            !
!                                                                                          !
!    Weissbluth, M. J.; Cotton, W. R., 1993: The representation of convection in mesoscale !
!         models. Part I: Scheme fabrication and calibration. J. Atmos. Sci, vol. 50(23),  !
!         3852-3872. (WC93)                                                                !
!                                                                                          !
! Their expression was derived based on the fact that:                                     !
!                 B ~= g (theta_v(cloud) - theta_v(envir.))/theta_v(envir.)                !
!                                                                                          !
! Even though we are using density, it is equivalent to theirs, because Tv is equivalent   !
! to density. As in (WC93), we are  neglecting the contribution of pressure difference     !
! between the draft and the environment, which may become more important at higher levels. !
! Our expression is just slightly different in formulation because WC93 approximated       !
! 1/(1+x) ~= 1-x, and we didn't do it here because there was no advantage to simplify in   !
! our case.                                                                                !
!---------------------------------------------------------------------------=--------------!
real function buoyancy_acc(rho_cup,rhoz_cld)

   use rconstants, only : g
   implicit none
   !----- Environment variables at cloud levels -------------------------------------------!
   real, intent(in) :: rho_cup   ! Environment density                             [ kg/m³]
   !----- Draft variables -----------------------------------------------------------------!
   real, intent(in) :: rhoz_cld  ! Draft density                                   [ kg/m³]
   !---------------------------------------------------------------------------------------!
   buoyancy_acc = g * (1. - rhoz_cld /rho_cup)

   return
end function buoyancy_acc
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine finds the estimates the PBL Height based on the turbulence and clouds. !
!------------------------------------------------------------------------------------------!
subroutine grell_find_pbl_height(mkx,mgmzp,z,tke,qliq,qice,pbllev)

   use grell_coms, only: pblhmax, rcpmin
   use rconstants, only: tkmin

   implicit none
   integer                  , intent(in)  :: mkx
   integer                  , intent(in)  :: mgmzp    ! Vertical dimension
   real   , dimension(mgmzp), intent(in)  :: z        ! Height at model level       [    m]
   real   , dimension(mgmzp), intent(in)  :: tke      ! Turb. Kinetic En. @ model   [ J/kg]
   real   , dimension(mgmzp), intent(in)  :: qliq     ! Liquid water mixing ratio   [kg/kg]
   real   , dimension(mgmzp), intent(in)  :: qice     ! Ice mixing ratio            [kg/kg]
   integer                  , intent(out) :: pbllev   ! Level of PBL top            [ ----]

   integer                                :: kpblmax  ! Maximum level for PBL.
   integer                                :: ktke_max ! Level of maximum cloud free TKE
   !---------------------------------------------------------------------------------------!
   ! First thing, find the maximum PBL height possible in terms of k levels                !
   !---------------------------------------------------------------------------------------!
   kpblmaxloop: do kpblmax=1,mkx
      if (z(kpblmax) >= pblhmax) exit kpblmaxloop
   end do kpblmaxloop

   !---------------------------------------------------------------------------------------!
   ! Now I look for the level of maximum TKE without clouds and below kzimax               !
   !---------------------------------------------------------------------------------------!
   ktke_max=maxloc(tke(1:kpblmax),dim=1,mask=qliq(1:kpblmax)+qice(1:kpblmax) < rcpmin)

   !---------------------------------------------------------------------------------------!
   !     Now I cycle between the level of maximum TKE and the maximum allowed height. as   !
   ! soon as I hit a level with little TKE and no cloud, that's going to be the PBL Height.!
   ! The "+1" is just because ktke_max is known to be the maximum.                         !
   !---------------------------------------------------------------------------------------!
   kpblloop: do pbllev=ktke_max+1,kpblmax
      if (tke(pbllev) <= 1.1*tkmin .and. qliq(pbllev)+qice(pbllev) < rcpmin) exit kpblloop
   end do kpblloop
   
   return
end subroutine grell_find_pbl_height
!==========================================================================================!
!==========================================================================================!
