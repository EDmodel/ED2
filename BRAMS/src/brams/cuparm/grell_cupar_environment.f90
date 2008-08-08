!==========================================================================================!
! grell_cupar_aux.f90                                                                      !
!                                                                                          !
!    This file contains subroutines that will calculate environment related stuff.         !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the values of some thermodynamic variables at the cloud level !
!------------------------------------------------------------------------------------------!
subroutine grell_thermo_cldlev(mkx,mgmzp,z_cup,t,q,p,tsur,qessur,qsur,hessur,hesur,psur    &
                              ,t_cup,qes_cup,q_cup,hes_cup,he_cup,p_cup,gamma_cup          )
   use rconstants, only : g,cp,alvl,aklv,lvordry,ep
   implicit none
   integer                  , intent(in)  :: mkx      ! # of Grell levels
   integer                  , intent(in)  :: mgmzp    ! # of levels for dimensions
   real   , dimension(mgmzp), intent(in)  :: z_cup    ! Height @ cloud levels        [    m]
   real   , dimension(mgmzp), intent(in)  :: t        ! Temperature @ model level    [    K]
   real   , dimension(mgmzp), intent(in)  :: q        ! Vapour mix. rat. @ model     [kg/kg]
   real   , dimension(mgmzp), intent(in)  :: p        ! Pressure @ model levels      [   Pa]
   real                     , intent(in)  :: tsur     ! Surface temperature          [    K]
   real                     , intent(in)  :: qessur   ! Surface Sat. mixing ratio    [kg/kg]
   real                     , intent(in)  :: qsur     ! Surface mixing ratio         [kg/kg]
   real                     , intent(in)  :: hessur   ! Surface sat. moist st. en.   [ J/kg]
   real                     , intent(in)  :: hesur    ! Surface moist static energy  [ J/kg]
   real                     , intent(in)  :: psur     ! Surface pressure             [   Pa]
   real   , dimension(mgmzp), intent(out) :: t_cup    ! Temperature @ cloud levels   [    K]
   real   , dimension(mgmzp), intent(out) :: qes_cup  ! Sat. mix. ratio @ cloud lev. [kg/kg]
   real   , dimension(mgmzp), intent(out) :: q_cup    ! Mixing ratio @ cloud levels  [kg/kg]
   real   , dimension(mgmzp), intent(out) :: hes_cup  ! Sat. MSE @ cloud levels      [ J/kg]
   real   , dimension(mgmzp), intent(out) :: he_cup   ! Moist static en. @ cloud lev.[ J/kg]
   real   , dimension(mgmzp), intent(out) :: p_cup    ! Pressure @ cloud levels      [   Pa]
   real   , dimension(mgmzp), intent(out) :: gamma_cup! Gamma, L/cp (dqes/dT)        [  ---]

   real   , external                      :: rslf     ! Function to find sat. mixing ratio
   
   integer                                :: k        ! Level counter
   !---------------------------------------------------------------------------------------!
   !    Here I will interpolate height, pressure, temperature and mixing ratio, and then   !
   ! recalculate the saturation mixing ratio, and moist static energy. By doing this we    !
   ! guarantee that thermodynamics will be consistent regarding phase changes.             !
   !---------------------------------------------------------------------------------------!
   do k = 2,mkx
      !----- Finding the temperature through simple interpolation -------------------------!
      t_cup  (k) = .5 *(t(k-1) + t (k))
      
      !----- Pressure is interpolating remembering that it has a log relationship with z --!
      p_cup  (k) = sqrt(p(k-1) * p (k))
      
      !----- Find mixing ratio to avoid super-saturation ----------------------------------!
      qes_cup(k) = max(1.e-8,rslf(p_cup(k),t_cup(k)))
      q_cup  (k) = min(.5 *(q(k-1)+q(k)),qes_cup(k))
      
      !----- Find the moist static energy at the cloud levels -----------------------------!
      he_cup (k)  = g * z_cup(k) + cp * t_cup(k) + alvl * q_cup(k)
      hes_cup(k)  = g * z_cup(k) + cp * t_cup(k) + alvl * qes_cup(k)

      !------------------------------------------------------------------------------------!
      ! Finding gamma. From Grell's paper, it is defined (after some algebra):             !
      !       _       _                   _          _                                     !
      !  Lv  |  d(qs)  |       Lv� qs    |  qs + eps  |    Lv� qs (qs + eps)               !
      ! ---- | ------- |  =  ----------  | ---------- | = -------------------              !
      !  Cp  |_  dT   _| p     Cp Rv T�  |_     eps  _|         Cp Ra T�                   !
      !------------------------------------------------------------------------------------!
      gamma_cup(k) = aklv*lvordry*qes_cup(k)*(qes_cup(k)+ep)/(t_cup(k)*t_cup(k) )

   end do
   
   !---------------------------------------------------------------------------------------!
   !    For the first level we will simply copy the surface values. This may be changed in !
   ! the future, to account for surface processes more effectively.                        !
   !---------------------------------------------------------------------------------------!
   t_cup    (1) = tsur
   qes_cup  (1) = qessur
   q_cup    (1) = qsur
   hes_cup  (1) = hessur
   he_cup   (1) = hesur
   p_cup    (1) = psur
   gamma_cup(1) = aklv*lvordry*qes_cup(1)*(qes_cup(1)+ep)/(t_cup(1)*t_cup(1) )
   
   return
end subroutine grell_thermo_cldlev
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes the values of some thermodynamic variables at the model      !
! level. This subroutine is currently being called to compute the variables modified by    !
! clouds.
!------------------------------------------------------------------------------------------!
subroutine grell_thermo_modlev(mkx,mgmzp,p,he,z,q,t,ql,qes,hes)
   use rconstants, only : cpi,g,cp,alvl,rgas
   implicit none
   integer                  , intent(in)    :: mkx   ! Number of levels
   integer                  , intent(in)    :: mgmzp ! Dimension in levels
   real   , dimension(mgmzp), intent(in)    :: p     ! Pressure
   real   , dimension(mgmzp), intent(in)    :: he    ! Moist static energy
   real   , dimension(mgmzp), intent(in)    :: z     ! Height
   
   real   , dimension(mgmzp), intent(inout) :: q     ! Mixing ratio - may be changed to
                                                     !   prevent supersaturation.

   real   , dimension(mgmzp), intent(out)   :: t     ! Temperature
   real   , dimension(mgmzp), intent(out)   :: ql    ! Condensed mixing ratio
   real   , dimension(mgmzp), intent(out)   :: qes   ! Saturation mixing ratio
   real   , dimension(mgmzp), intent(out)   :: hes   ! Saturation moist static energy
   
   real   , external                        :: rslf  ! Saturation mixing ratio function

   integer                                  :: k     ! Level counter
   real                                     :: qtemp ! Temporary mixing ratio.
   do k = 1,mkx
      !------------------------------------------------------------------------------------!
      !     Pressure is retrieved assuming that the density is almost the same as before.  !
      ! This should be done until we have a set of pressure, temperature and vapour mixing !
      ! ratio that is in equilibrium or not saturated (super-saturated). If it happens to  !
      ! be supersaturated, then I remove the excess and re-compute the pressure and sat.   !
      ! mixing ratio.                                                                      !
      !------------------------------------------------------------------------------------!
      qtemp   = 9999.
      qes(k)  = q(k)      !---- First guess for qes is q...
      do while ( qtemp-qes(k) > 1.e-6)
         !---- The first time this will be q(k), if not saturated it won't cycle ----------!
         qtemp  = qes(k) 
         t(k)   = cpi * (he(k) - g *z(k) - alvl * qtemp)
         qes(k) = max(1.e-8,rslf(p(k),t(k)))
      end do
      
      !----- Fixing the mixing ratio in case it was super-saturated -----------------------!
      if (qes(k) < q(k)) then
         ql(k) = q(k)-qes(k)
         q(k)  = qes(k)
      else
         ql(k) = 0.
      end if

      !----- Find the moist static energy at the cloud levels -----------------------------!
      hes(k)  = g * z(k) + cp * t(k) + alvl * qes(k)

   end do
   
   return
end subroutine grell_thermo_modlev
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine finds the estimates the PBL Height based on the turbulence and clouds. !
!------------------------------------------------------------------------------------------!
subroutine grell_find_pbl_height(mkx,mgmzp,z,tkeg,rcpg,pbllev)

   use grell_coms, only: pblhmax, rcpmin
   use rconstants, only: tkmin

   implicit none
   integer                  , intent(in)  :: mkx
   integer                  , intent(in)  :: mgmzp    ! Vertical dimension
   real   , dimension(mgmzp), intent(in)  :: z        ! Height at model level       [    m]
   real   , dimension(mgmzp), intent(in)  :: tkeg     ! Turb. Kinetic En. @ model   [ J/kg]
   real   , dimension(mgmzp), intent(in)  :: rcpg     ! Cloud droplet mixing        [kg/kg]
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
   ktke_max=maxloc(tkeg(1:kpblmax),dim=1,mask=rcpg(1:kpblmax) < rcpmin)

   !---------------------------------------------------------------------------------------!
   !     Now I cycle between the level of maximum TKE and the maximum allowed height. as   !
   ! soon as I hit a level with little TKE and no cloud, that's going to be the PBL Height.!
   ! The "+1" is just because ktke_max is known to be the maximum.                         !
   !---------------------------------------------------------------------------------------!
   kpblloop: do pbllev=ktke_max+1,kpblmax
      if (tkeg(pbllev) <= 1.1*tkmin .and. rcpg(pbllev) < rcpmin) exit kpblloop
   end do kpblloop
   
   return
end subroutine grell_find_pbl_height
!==========================================================================================!
!==========================================================================================!
