!==========================================================================================!
! grell_cupar_updraft.f90                                                                  !
!                                                                                          !
!    This file contains subroutines that will calculate updraft related stuff.             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine finds the level of origin of updrafts. Depending on the profile, this  !
! may not be the definite level, it may change during the search for the level of free     !
! convection.                                                                              !
!------------------------------------------------------------------------------------------!
subroutine grell_updraft_origin(mkx,mgmzp,iupmethod,kpbl,kbmax,z,tke,qice,qliq,theiv_cup   &
                               ,ierr,k22)

   implicit none
   integer               , intent(in)     :: mkx,mgmzp ! Grid dimensions
   integer               , intent(in)     :: iupmethod ! Method to find the level
   integer               , intent(in)     :: kpbl      ! Level of PBL top if pre-computed.
   integer               , intent(in)     :: kbmax     ! Maximum allowed level for updraft
   !----- Variables at the model levels ---------------------------------------------------!
   real, dimension(mgmzp), intent(in)     :: z         ! Height                    [     m]
   real, dimension(mgmzp), intent(in)     :: tke       ! Turbulent Kinetic Energy  [  J/kg]
   real, dimension(mgmzp), intent(in)     :: qliq      ! Liquid water mixing ratio [ kg/kg]
   real, dimension(mgmzp), intent(in)     :: qice      ! Ice mixing ratio          [ kg/kg]
   !----- Variables at the cloud levels ---------------------------------------------------!
   real, dimension(mgmzp), intent(in)     :: theiv_cup ! Thetae_iv                 [     K]
   !----- Output and sort of output variables ---------------------------------------------!
   integer               , intent(inout)  :: ierr      ! Error flag
   integer               , intent(out)    :: k22       ! Updraft origin level
   !----- Constant. Avoding using too levels too close to the surface ---------------------!
   integer               , parameter      :: kstart=2  ! Minimum level
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Four possibilities, depending on the user's preference:                            !
   !    a. The user wants to use the PBL height (iupmethod=2), and the user is running     !
   !       Nakanishi/Niino turbulence (positive pblidx, previously computed);              !
   !    b. The user wants to use the PBL height (iupmethod=2), other turbulence was used   !
   !    (pblidx is always zero in this case). Find PBL here.                               !
   !    c. The user wants to use the maximum thetae_iv as the first guess, for maximum     !
   !       instability.                                                                    !
   !    d. The user wants to use the most turbulent level as the updraft origin.           !
   !---------------------------------------------------------------------------------------!
   select case (iupmethod)
   case (1) ! Maximum Thetae_iv
      k22 = (kstart-1) + maxloc(theiv_cup(kstart:kbmax),dim=1)
   case (2) ! PBL top  
      if (kpbl /= 0) then
         k22 = kpbl
      else
         call grell_find_pbl_height(mkx,mgmzp,z,tke,qliq,qice,k22)
      end if
   case (3) ! Most turbulent
      k22 = (kstart-1) + maxloc(tke(kstart:kbmax),dim=1)
   end select
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !   If the level of updraft origin is too high, then cumulus parameterization should    !
   ! not be called here.                                                                   !
   !---------------------------------------------------------------------------------------!
   if (k22 >= kbmax) ierr = 2

   return
end subroutine grell_updraft_origin
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine finds the level of free convection. Note that this subroutine is       !
! recursive, and this is the case because if we don't like the distance between the cloud  !
! base and the level that updrafts origin, we may push the latter up and try again.        !
!------------------------------------------------------------------------------------------!
recursive subroutine grell_find_cloud_lfc(mkx,mgmzp,kbmax,cap_max,wnorm_max,wwind,sigw     &
                                         ,exner_cup,p_cup,theiv_cup,thil_cup,t_cup         &
                                         ,qtot_cup,qvap_cup,qliq_cup,qice_cup,qsat_cup     &
                                         ,rho_cup,dzd_cld,mentru_rate,theivu_cld,thilu_cld &
                                         ,tu_cld,qtotu_cld,qvapu_cld,qliqu_cld,qiceu_cld   &
                                         ,qsatu_cld,rhou_cld,dbyu,k22,ierr,kbcon,wbuoymin)
   use rconstants, only : epi,rgas
   use therm_lib , only : idealdens
   implicit none

   integer, intent(in)                  :: mkx        ! # of vertical layers
   integer, intent(in)                  :: mgmzp      ! Vertical dimension
   integer, intent(in)                  :: kbmax      ! Top level allowed for LFC
   real   , intent(in)                  :: cap_max    ! Depth of capping inversion  [   Pa]
   real   , intent(in)                  :: wnorm_max  ! Maximum normalised velocity in
   !----- Input environment variables -----------------------------------------------------!
   real, dimension(mgmzp), intent(in)   :: wwind      ! Vertical velocity          [   m/s]
   real, dimension(mgmzp), intent(in)   :: sigw       ! wwind standard deviation   [   m/s]
   real, dimension(mgmzp), intent(in)   :: exner_cup  ! Exner f. @ cloud level     [J/kg/K]
   real, dimension(mgmzp), intent(in)   :: p_cup      ! Pressure @ cloud level     [    Pa]
   real, dimension(mgmzp), intent(in)   :: theiv_cup  ! Thetae_iv                  [     K]
   real, dimension(mgmzp), intent(in)   :: thil_cup   ! Theta_il                   [     K]
   real, dimension(mgmzp), intent(in)   :: t_cup      ! Temperature                [     K]
   real, dimension(mgmzp), intent(in)   :: qtot_cup   ! Total mixing ratio         [ kg/kg]
   real, dimension(mgmzp), intent(in)   :: qvap_cup   ! Vapour mixing ratio        [ kg/kg]
   real, dimension(mgmzp), intent(in)   :: qliq_cup   ! Liquid mixing ratio        [ kg/kg]
   real, dimension(mgmzp), intent(in)   :: qice_cup   ! Ice mixing ratio           [ kg/kg]
   real, dimension(mgmzp), intent(in)   :: qsat_cup   ! Sat. vapour mixing ratio   [ kg/kg]
   real, dimension(mgmzp), intent(in)   :: rho_cup    ! Density                    [ kg/m³]
   real, dimension(mgmzp), intent(in)   :: dzd_cld    ! Top-bottom cloud thickness [     m]
   real, dimension(mgmzp), intent(in)   :: mentru_rate! Entrainment rate           [   1/m]
   !----- Updraft variables. These will are scratch now, they will be adjusted later on ---!
   real, dimension(mgmzp), intent(inout):: theivu_cld ! Updraft theta_il           [     K]
   real, dimension(mgmzp), intent(inout):: thilu_cld  ! Updraft theta_il           [     K]
   real, dimension(mgmzp), intent(inout):: tu_cld     ! Updraft temperature        [     K]
   real, dimension(mgmzp), intent(inout):: qtotu_cld  ! Updraft total mixing ratio [ kg/kg]
   real, dimension(mgmzp), intent(inout):: qvapu_cld  ! Updraft vapour mix. ratio  [ kg/kg]
   real, dimension(mgmzp), intent(inout):: qliqu_cld  ! Updraft liquid mix.  ratio [ kg/kg]
   real, dimension(mgmzp), intent(inout):: qiceu_cld  ! Updraft ice    mix.  ratio [ kg/kg]
   real, dimension(mgmzp), intent(inout):: qsatu_cld  ! Updraft satur. mix. ratio  [ kg/kg]
   real, dimension(mgmzp), intent(inout):: rhou_cld   ! Updraft density            [ kg/m³]
   real, dimension(mgmzp), intent(inout):: dbyu       ! Buoyancy acceleration      [  m/s²]
   !----- These variables may or may not be assigned here so use inout --------------------!
   integer               , intent(inout):: k22        ! Level of origin of updrafts
   integer               , intent(inout):: ierr       ! Error flag
   !----- Output variable -----------------------------------------------------------------!
   integer               , intent(out)  :: kbcon      ! Level of free convection
   real                  , intent(out)  :: wbuoymin   ! Min. buoyanct velocity     [   m/s]
   !----- External functions --------------------------------------------------------------!
   real   , external                    :: buoyancy_acc ! Buoyancy acceleration funtion.
   !----- Local variables -----------------------------------------------------------------!
   integer                              :: k       ! Counter
   real                                 :: pcdiff  ! Pres. diff. b/w k22 and kbcon [    Pa]
   logical                              :: pushup  ! Push k22 upwards              [   T/F]
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Start with the first guess: the Level of free convection (kbcon) is at the level   !
   ! in which updrafts originate.                                                          !
   !---------------------------------------------------------------------------------------!
   pushup = .false.
   kbcon  = k22
   
   !---------------------------------------------------------------------------------------!
   !    Initialise draft-related variables. Some levels will be overwritten later on.      !
   !---------------------------------------------------------------------------------------!
   theivu_cld = theiv_cup
   thilu_cld  = thil_cup
   tu_cld     = t_cup
   qtotu_cld  = qtot_cup
   qvapu_cld  = qvap_cup
   qliqu_cld  = qliq_cup
   qiceu_cld  = qice_cup
   qsatu_cld  = qsat_cup
   rhou_cld   = rho_cup
   dbyu       = 0.  
   
   !---------------------------------------------------------------------------------------!
   !   First step: finding the level in which the air lifted from the level that updrafts  !
   ! originate would become buoyant.                                                       !
   !---------------------------------------------------------------------------------------!
   kbconloop: do
      kbcon=kbcon+1
      if (kbcon > kbmax + 2) then 
         !------ Gave up... Cloud would be too high to be a cumulus. ----------------------!
         ierr = 3
         return
      end if
      
      !----- No entrainment below the LFC, so theiv, thil, and qtot shouldn't change ------!
      theivu_cld(kbcon) = theiv_cup(k22)
      thilu_cld (kbcon) = thil_cup (k22)
      qtotu_cld (kbcon) = qtot_cup (k22)
      qliqu_cld (kbcon) = qliqu_cld(kbcon-1) !----- For faster convergence ----------------!
      qiceu_cld (kbcon) = qiceu_cld(kbcon-1) !----- For faster convergence ----------------!
      !------ Finding a consistent set of temperature and mixing ratios -------------------!
      call thil2tqall(thilu_cld(kbcon),exner_cup(kbcon),p_cup(kbcon),qtotu_cld(kbcon)      &
                     ,qliqu_cld(kbcon),qiceu_cld(kbcon),tu_cld(kbcon),qvapu_cld(kbcon)     &
                     ,qsatu_cld(kbcon))
      !------ Finding the draft density, assuming pu_cld = p_cup... -----------------------!
      rhou_cld(kbcon) = idealdens(p_cup(kbcon),tu_cld(kbcon),qvapu_cld(kbcon)              &
                                 ,qtotu_cld(kbcon))

      !------ Finding buoyancy ------------------------------------------------------------!
      dbyu(kbcon) = buoyancy_acc(rho_cup(kbcon),rhou_cld(kbcon))
      !------ First guess for buoyancy. LFC is the first one to have it positive ----------!
      if (dbyu(kbcon) > 0. ) exit kbconloop

   end do kbconloop

   !---------------------------------------------------------------------------------------!
   !     The LFC is really somewhere betweeen kbcon-1 and kbcon. If we assume it is kbcon  !
   ! we are going to make it less likely for convection to happen, whereas assuming it to  !
   ! be kbcon-1 may cause convection where it shouldn't really happen.                     !
   !---------------------------------------------------------------------------------------!
   kbcon = kbcon - 1

   !---------------------------------------------------------------------------------------!
   !   Finding the minimum velocity an updraft would need to have to reach the tentative   !
   ! LFC. This may not be needed here depending on the test the user asked for, but it     !
   ! will be needed when computing the fractional area covered by clouds.                  !
   !---------------------------------------------------------------------------------------!
   wbuoymin=0.
   do k=kbcon-1,k22,-1
      !------------------------------------------------------------------------------------!
      !     The only risk for this sqrt have negative value is if positive buoyancy at     !
      ! kbcon. But in this case, the LFC is actually closer to kbcon-1 rather than kbcon,  !
      ! so keep it zero again. Again, the updraft experiences no entrainment nor           !
      ! detrainment below kbcon, which simplifies the equation below. However, we must     !
      ! account the friction, which is defined based on the constant entrainment rate,     !
      ! following Zhang and Fritsch (1986) parametrisation).                               !
      !------------------------------------------------------------------------------------!
      wbuoymin=sqrt(max(0., wbuoymin*wbuoymin - (dbyu(k+1)+dbyu(k))*dzd_cld(k)))
   end do
   !---------------------------------------------------------------------------------------!
   !    Now I use either cap_max or wnorm_max to decide whether convection can happen with !
   ! the k22/kbcon combination we got at this point.                                       !
   !---------------------------------------------------------------------------------------!
   if (cap_max > 0.) then
      !------------------------------------------------------------------------------------!
      !    So the LFC is kind of far from the level in which updrafts originate. If there  !
      ! is an inversion capping entirely within the layer, then I may not develop any      !
      ! convection, so I try an alternative approach. If the "gap" between k22 and kbcon   !
      ! is larger than the maximum depth of the inversion capping, I will try to push k22  !
      ! upwards. If kbcon is just the next level, take it.                                 !
      !------------------------------------------------------------------------------------!
      pcdiff = p_cup(k22) - p_cup(kbcon) !----- Pressure decreases with height... ---------!
      pushup = (pcdiff > cap_max) .and. (kbcon-k22 > 1)
   else
      !------------------------------------------------------------------------------------!
      !    We now compare this buoyant velocity with the velocity at k22. If this is too   !
      ! unlikely to happen, move k22 up and give one more try, otherwise, we found the     !
      ! combination.                                                                       !
      !------------------------------------------------------------------------------------!
      pushup = wbuoymin > wwind(k22) + wnorm_max*sigw(k22)
   end if

   if (pushup) then
      k22 = k22 + 1
      call grell_find_cloud_lfc(mkx,mgmzp,kbmax,cap_max,wnorm_max,wwind,sigw,exner_cup     &
                              ,p_cup,theiv_cup,thil_cup,t_cup,qtot_cup,qvap_cup,qliq_cup   &
                              ,qice_cup,qsat_cup,rho_cup,dzd_cld,mentru_rate,theivu_cld    &
                              ,thilu_cld,tu_cld,qtotu_cld,qvapu_cld,qliqu_cld,qiceu_cld    &
                              ,qsatu_cld,rhou_cld,dbyu,k22,ierr,kbcon,wbuoymin)

   end if

   return
end subroutine grell_find_cloud_lfc
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine should be used only for "0" and "x" thermodynamics. This is using the  !
! k22 and kbcon already found, and computing the buoyancy at this lowest part.             !
!------------------------------------------------------------------------------------------!
subroutine grell_buoy_below_lfc(mkx,mgmzp,k22,kbcon,exner_cup,p_cup,theiv_cup,thil_cup     &
                               ,t_cup,qtot_cup,qvap_cup,qliq_cup,qice_cup,qsat_cup,rho_cup &
                               ,theivu_cld,thilu_cld,tu_cld,qtotu_cld,qvapu_cld,qliqu_cld  &
                               ,qiceu_cld,qsatu_cld,rhou_cld,dbyu)
   use rconstants, only : epi,rgas
   use therm_lib , only : idealdens
   implicit none

   integer               , intent(in)   :: mkx        ! # of vertical layers
   integer               , intent(in)   :: mgmzp      ! Vertical dimension
   integer               , intent(in)   :: k22        ! Level of origin of updrafts
   integer               , intent(in)   :: kbcon      ! Level of free convection
   !----- Input environment variables -----------------------------------------------------!
   real, dimension(mgmzp), intent(in)   :: exner_cup  ! Exner f. @ cloud level     [J/kg/K]
   real, dimension(mgmzp), intent(in)   :: p_cup      ! Pressure @ cloud level     [    Pa]
   real, dimension(mgmzp), intent(in)   :: theiv_cup  ! Thetae_iv                  [     K]
   real, dimension(mgmzp), intent(in)   :: thil_cup   ! Theta_il                   [     K]
   real, dimension(mgmzp), intent(in)   :: t_cup      ! Temperature                [     K]
   real, dimension(mgmzp), intent(in)   :: qtot_cup   ! Total mixing ratio         [ kg/kg]
   real, dimension(mgmzp), intent(in)   :: qvap_cup   ! Vapour mixing ratio        [ kg/kg]
   real, dimension(mgmzp), intent(in)   :: qliq_cup   ! Liquid mixing ratio        [ kg/kg]
   real, dimension(mgmzp), intent(in)   :: qice_cup   ! Ice mixing ratio           [ kg/kg]
   real, dimension(mgmzp), intent(in)   :: qsat_cup   ! Sat. vapour mixing ratio   [ kg/kg]
   real, dimension(mgmzp), intent(in)   :: rho_cup    ! Density                    [ kg/m³]
   !----- Updraft variables. --------------------------------------------------------------!
   real, dimension(mgmzp), intent(inout):: theivu_cld ! Updraft theta_il           [     K]
   real, dimension(mgmzp), intent(inout):: thilu_cld  ! Updraft theta_il           [     K]
   real, dimension(mgmzp), intent(inout):: tu_cld     ! Updraft temperature        [     K]
   real, dimension(mgmzp), intent(inout):: qtotu_cld  ! Updraft total mixing ratio [ kg/kg]
   real, dimension(mgmzp), intent(inout):: qvapu_cld  ! Updraft vapour mix. ratio  [ kg/kg]
   real, dimension(mgmzp), intent(inout):: qliqu_cld  ! Updraft liquid mix.  ratio [ kg/kg]
   real, dimension(mgmzp), intent(inout):: qiceu_cld  ! Updraft ice    mix.  ratio [ kg/kg]
   real, dimension(mgmzp), intent(inout):: qsatu_cld  ! Updraft satur. mix. ratio  [ kg/kg]
   real, dimension(mgmzp), intent(inout):: rhou_cld   ! Updraft density            [ kg/m³]
   real, dimension(mgmzp), intent(inout):: dbyu       ! Buoyancy acceleration      [  m/s²]
   !----- External functions --------------------------------------------------------------!
   real   , external                    :: buoyancy_acc ! Buoyancy acceleration funtion.
   !----- Local variables -----------------------------------------------------------------!
   integer                              :: k       ! Counter
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   ! 1. Initialise draft-related variables. Some levels will be overwritten later on.      !
   !---------------------------------------------------------------------------------------!
   theivu_cld = theiv_cup
   thilu_cld  = thil_cup
   tu_cld     = t_cup
   qtotu_cld  = qtot_cup
   qvapu_cld  = qvap_cup
   qliqu_cld  = qliq_cup
   qiceu_cld  = qice_cup
   qsatu_cld  = qsat_cup
   rhou_cld   = rho_cup
   dbyu       = 0.  
   
   !---------------------------------------------------------------------------------------!
   ! 2. Between the updraft origin and the level of free convection we assume no entrain-  !
   !    ment, so theivu_cld, thilu_cld and qtotu_cld are conserved. Then we find the other !
   !    variables.                                                                         !
   !---------------------------------------------------------------------------------------!
   do k=k22,kbcon
      theivu_cld(k) = theiv_cup(k22)
      thilu_cld (k) = thil_cup (k22)
      qtotu_cld (k) = qtot_cup (k22)
      !------ Finding a consistent set of temperature and mixing ratios -------------------!
      call thil2tqall(thilu_cld(k),exner_cup(k),p_cup(k),qtotu_cld(k),qliqu_cld(k)         &
                     ,qiceu_cld(k),tu_cld(k),qvapu_cld(k),qsatu_cld(k))
      !------ Finding the draft density, assuming pu_cld = p_cup... -----------------------!
      rhou_cld(k) = idealdens(p_cup(k),tu_cld(k),qvapu_cld(k),qtotu_cld(k))

      !------ Finding buoyancy ------------------------------------------------------------!
      dbyu(k) = buoyancy_acc(rho_cup(k),rhou_cld(k))
   end do

   return
end subroutine grell_buoy_below_lfc
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine finds the level of convective cloud top. The method as is has two    !
! approximations:                                                                          !
! 1. The cloud top is roughly the equilibrium level (the cloud top is often higher because !
!    the air parcel still has momentum at EL);                                             !
! 2. We can determine the equilibrium level from thetae_iv only. The most accurate way     !
!    would be using the buoyancy, however, we need to know the  cloud top before proceed-  !
!    ing to the moisture balance. Nonetheless, the buoyancy will be correctly computed for !
!    cloud work function after we know the updraft water phase partition.                  !
!------------------------------------------------------------------------------------------!
subroutine grell_find_cloud_top(mkx,mgmzp,kbcon,cdu,mentru_rate,theiv_cup,theivs_cup       &
                               ,dzu_cld,theivu_cld,dbyu,ierr,ktop)

   implicit none
   !----- Input variables -----------------------------------------------------------------!
   integer, intent(in)                   :: mkx,mgmzp   ! Grid dimensions
   integer, intent(in)                   :: kbcon       ! Level of free convection
   real, dimension(mgmzp), intent(in)    :: cdu         ! Detrainment function      [  1/m]
   real, dimension(mgmzp), intent(in)    :: mentru_rate ! Entrainment rate          [  1/m]

   real, dimension(mgmzp), intent(in)    :: theiv_cup   ! Thetae_iv @ cloud levels  [    K]
   real, dimension(mgmzp), intent(in)    :: theivs_cup  ! Thetaes_iv @ cloud levels [    K]
   real, dimension(mgmzp), intent(in)    :: dzu_cld     ! Layer thickness           [    m]
   !----- Transit variables, these will be reset above the cloud top ----------------------!
   real, dimension(mgmzp), intent(inout) :: theivu_cld  ! Thetae_iv                 [    K]
   real, dimension(mgmzp), intent(inout) :: dbyu        ! Buoyancy acceleration     [ m/s²]
   !----- Output variables, error and the cloud top ---------------------------------------!
   integer               , intent(inout) :: ierr        ! Error flag
   integer               , intent(out)   :: ktop        ! Cloud top level
   !----- Auxiliary variables -------------------------------------------------------------!
   logical                               :: foundit     ! Top check Boolean flag
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !    Initialize foundit                                                                 !
   !---------------------------------------------------------------------------------------!
   foundit=.false.
   
   !---------------------------------------------------------------------------------------!
   !   Loop over levels looking for the first level above the cloud top that has negative  !
   ! buoyancy.                                                                             !
   !---------------------------------------------------------------------------------------!
   ktoploop: do ktop=kbcon+1,mkx-2
      foundit=theivu_cld(ktop) < theivs_cup(ktop)
      if (foundit) exit ktoploop
   end do ktoploop
 
   !---------------------------------------------------------------------------------------!
   !    The cloud top truthfully should be the level in which it loses its positive verti- !
   ! cal velocity, but as a first approximation we can leave it at the equilibrium level.  !
   !---------------------------------------------------------------------------------------!
   ktop = ktop - 1

   if (foundit) then
      !------------------------------------------------------------------------------------!
      !    Fixing the updraft properties aloft to be the same as the environment.          !
      !------------------------------------------------------------------------------------!
      theivu_cld(ktop+1:mkx) = theiv_cup(ktop+1:mkx)
      dbyu(ktop+1:mkx)       = 0.
   else
      !------------------------------------------------------------------------------------!
      !    I won't allow those clouds that extend all the way up to the mesosphere, they   !
      ! are way too scary...                                                               !
      !------------------------------------------------------------------------------------!
      ktop=0
      dbyu=0.
      theivu_cld = theiv_cup
      ierr=5
   end if

   return
end subroutine grell_find_cloud_top
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the incloud moist static energy                               !
!------------------------------------------------------------------------------------------!
subroutine grell_theiv_updraft(mkx,mgmzp,k22,kbcon,cdu,mentru_rate,theiv,theiv_cup,dzu_cld &
                              ,theivu_cld)
   implicit none
   integer               , intent(in)    :: mkx, mgmzp  ! Grid dimesnsions
   integer               , intent(in)    :: k22         ! Level in which updrafts begin
   integer               , intent(in)    :: kbcon       ! Level of free convection

   real, dimension(mgmzp), intent(in)    :: cdu         ! Updraft detrainment function;
   real, dimension(mgmzp), intent(in)    :: mentru_rate ! Updraft entrainment rate
   real, dimension(mgmzp), intent(in)    :: theiv       ! Thetae_iv @ model levels;
   real, dimension(mgmzp), intent(in)    :: theiv_cup   ! Thetae_iv @ cloud levels;
   real, dimension(mgmzp), intent(in)    :: dzu_cld     ! Delta-z for updrafts;
   real, dimension(mgmzp), intent(inout) :: theivu_cld  ! Cloud thetae_iv

   integer                             :: k           ! Counter
   !---------------------------------------------------------------------------------------!
   !    Below the LFC thil is already set up, move on to the entrainment layer. Inside the !
   ! cloud, the updraft will no longer conserve thetae_iv because there is entrainment and !
   ! detrainment (the phase change only doesn't change thetae_iv).                         !
   !---------------------------------------------------------------------------------------!
   do k=kbcon+1,mkx-1
      theivu_cld(k) = (theivu_cld(k-1)*(1.-.5*cdu(k)*dzu_cld(k))                           &
                      + mentru_rate(k)*dzu_cld(k)*theiv(k-1))                              &
                      / (1.+mentru_rate(k)*dzu_cld(k) - .5*cdu(k)*dzu_cld(k))
   end do

   return
end subroutine grell_theiv_updraft
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the normalized mass flux associated with updrafts             !
!------------------------------------------------------------------------------------------!
subroutine grell_nms_updraft(mkx,mgmzp,k22,kbcon,ktop,mentru_rate,cdu,dzu_cld,etau_cld)
   implicit none

   integer               , intent(in)    :: mkx, mgmzp  ! Grid dimesnsions
   integer               , intent(in)    :: k22         ! Level in which updrafts begin
   integer               , intent(in)    :: kbcon       ! Level of free convection
   integer               , intent(in)    :: ktop        ! Cloud top
                                        
   real, dimension(mgmzp), intent(in)    :: mentru_rate ! Updraft entrainment rate
   real, dimension(mgmzp), intent(in)    :: cdu         ! Updraft detrainment function;
   real, dimension(mgmzp), intent(in)    :: dzu_cld     ! Delta-z for clouds           
   real, dimension(mgmzp), intent(inout) :: etau_cld    ! Normalized updraft flux

   integer                             :: k           ! Counter

   !---------------------------------------------------------------------------------------!
   ! 1. Below the updraft origin there is no upward mass flux, set it to zero.             !
   !---------------------------------------------------------------------------------------!
   etau_cld(1:(k22-1)) = 0.
   
   
   
   !---------------------------------------------------------------------------------------!
   ! 2. There is no entrainment/detrainment between the updraft origin and the level of    !
   !    free convection,so I assume that the normalized mass flux is one.                  !
   !---------------------------------------------------------------------------------------!
   etau_cld(k22:kbcon) = 1.



   !---------------------------------------------------------------------------------------!
   ! 3. Between the LFC and cloud top, need to consider entrainment and detrainment        !
   !    contributions, loop through levels.                                                !
   !---------------------------------------------------------------------------------------!
   do k=kbcon+1,ktop
      etau_cld(k)=etau_cld(k-1)*(1.+(mentru_rate(k)-cdu(k))*dzu_cld(k))
   end do


   !---------------------------------------------------------------------------------------!
   ! 1. Above the cloud top, no mass flux, set it to zero.                                 !
   !---------------------------------------------------------------------------------------!
   etau_cld(ktop+1:mkx) = 0.
   
   return
end subroutine grell_nms_updraft
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes most thermodynamic variables associated with the updraft,   !
! in particular those affected by phase change.                                            !
!------------------------------------------------------------------------------------------!
subroutine grell_most_thermo_updraft(comp_down,mkx,mgmzp,kbcon,ktop,cdu,mentru_rate,qtot   &
                                    ,p_cup,exner_cup,thil_cup,t_cup,qtot_cup,qvap_cup      &
                                    ,qliq_cup,qice_cup,qsat_cup,rho_cup,theivu_cld         &
                                    ,etau_cld,dzu_cld,thilu_cld,tu_cld,qtotu_cld,qvapu_cld &
                                    ,qliqu_cld,qiceu_cld,qsatu_cld,rhou_cld,dbyu,pwu_cld   &
                                    ,pwavu)
   use rconstants, only : epi,rgas, t00, toodry
   use therm_lib , only : thetaeiv2thil, idealdens, toler, maxfpo
   implicit none
   logical               , intent(in)    :: comp_down  ! Flag for downdraft/precipitation
   integer               , intent(in)    :: mkx        ! Levels 
   integer               , intent(in)    :: mgmzp      ! Levels
   integer               , intent(in)    :: kbcon      ! Level of free convection
   integer               , intent(in)    :: ktop       ! Cloud top

   real, dimension(mgmzp), intent(in)    :: cdu         ! Detrainment function     [   1/m]
   real, dimension(mgmzp), intent(in)    :: mentru_rate ! Entrainment function     [   1/m]
   !----- Variables at model levels -------------------------------------------------------!
   real, dimension(mgmzp), intent(in)    :: qtot        ! Total mixing ratio       [ kg/kg]
   !----- Variables at cloud levels -------------------------------------------------------!
   real, dimension(mgmzp), intent(in)    :: p_cup       ! Pressure @ cloud levels  [   1/m]
   real, dimension(mgmzp), intent(in)    :: exner_cup   ! Exner fctn. @ cloud lev. [J/kg/K]
   real, dimension(mgmzp), intent(in)    :: thil_cup    ! Theta_il                 [     K]
   real, dimension(mgmzp), intent(in)    :: t_cup       ! Temperature              [     K]
   real, dimension(mgmzp), intent(in)    :: qtot_cup    ! Total mixing ratio       [ kg/kg]
   real, dimension(mgmzp), intent(in)    :: qvap_cup    ! Vapour mixing ratio      [ kg/kg]
   real, dimension(mgmzp), intent(in)    :: qliq_cup    ! Liquid mixing ratio      [ kg/kg]
   real, dimension(mgmzp), intent(in)    :: qice_cup    ! Ice mixing ratio         [ kg/kg]
   real, dimension(mgmzp), intent(in)    :: qsat_cup    ! Sat. mixing ratio        [ kg/kg]
   real, dimension(mgmzp), intent(in)    :: rho_cup     ! Density                  [ kg/m³]
   !----- Input variables at updraft ------------------------------------------------------!
   real, dimension(mgmzp), intent(in)    :: theivu_cld  ! Thetae_iv @ updraft      [     K]
   real, dimension(mgmzp), intent(in)    :: etau_cld    ! Normalized mass flux     [   ---]
   real, dimension(mgmzp), intent(in)    :: dzu_cld     ! Layer thickness          [     m]
   !----- Transit variables, which will be changed between kbcon and ktop -----------------!
   real, dimension(mgmzp), intent(inout) :: thilu_cld   ! Theta_il                 [     K]
   real, dimension(mgmzp), intent(inout) :: tu_cld      ! Temperature              [     K]
   real, dimension(mgmzp), intent(inout) :: qtotu_cld   ! Total mixing ratio       [ kg/kg]
   real, dimension(mgmzp), intent(inout) :: qvapu_cld   ! Vapour mixing ratio      [ kg/kg]
   real, dimension(mgmzp), intent(inout) :: qliqu_cld   ! Liquid mixing ratio      [ kg/kg]
   real, dimension(mgmzp), intent(inout) :: qiceu_cld   ! Ice mixing ratio         [ kg/kg]
   real, dimension(mgmzp), intent(inout) :: qsatu_cld   ! Sat. mixing ratio        [ kg/kg]
   real, dimension(mgmzp), intent(inout) :: rhou_cld    ! Density                  [ kg/m³]
   real, dimension(mgmzp), intent(inout) :: dbyu        ! Buoyancy acceleration    [  m/s²]
   !----- Output variables ----------------------------------------------------------------!
   real, dimension(mgmzp), intent(inout) :: pwu_cld     ! Fall-out rain            [ kg/kg]
   real                  , intent(out)   :: pwavu       ! Total normalized integrated cond.
   !----- Local variables -----------------------------------------------------------------!
   integer                :: k              ! Counter                              [  ----]
   integer                :: it             ! Iteration counter                    [  ----]
   real                   :: c0             ! Rainfall conversion parameter. As pointed out
                                            !   in Grell (1993), this could be a function 
                                            !   of wind shear and cloud spectral size.
   logical                :: converged      ! Flag to test convergence             [   T|F]
   logical                :: bisection      ! Flag to use bisection                [   T|F]
   real                   :: qeverything    ! Sum of all mixing ratios             [ kg/kg]
   real                   :: qtotua, qtotuz ! Aux. vars for bisection iteration    [ kg/kg]
   real                   :: qtotuc, qtotup ! Aux. vars for secant iteration       [ kg/kg]
   real                   :: funa, funz     ! Function evaluation for bisection    [ kg/kg]
   real                   :: func, funp     ! Function evaluation for secant       [ kg/kg]
   real                   :: funnow         ! Function at this iteration           [ kg/kg]
   real                   :: denomin        ! Denominator, just to clean the eqn.  [  ----]
   real                   :: denomini       ! 1./denominator                       [  ----]
   real, dimension(mgmzp) :: leftu_cld      ! Total condensation that left cloud   [ kg/kg]
   real                   :: tubis          ! Scratch var. for temperature         [     K]
   real                   :: delta          ! Aux. var. for bisection 2nd guess    [ kg/kg]

   !----- External functions --------------------------------------------------------------!
   real, external                     :: buoyancy_acc   ! Buoyancy acceleration funtion.
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! 1. First of all, I'll check whether this is a precipitating cloud or not. Currently,  !
   !    precipitating clouds are the ones that contain downdrafts. This does not need to   !
   !    be the requirement, it was just for convenience. The conversion rate from cloud to !
   !    rain should be a function of the cloud size and wind shear, but it is only a step  !
   !    function, 0 if the radius is small, or a non-zero constant otherwise. Perhaps this !
   !    should be done in a different way in the future.                                   !
   !---------------------------------------------------------------------------------------!
   if (comp_down) then
      c0 = 0.002
   else 
      c0 = 0.
   end if
   
   !---------------------------------------------------------------------------------------!
   ! 2.  Initialise integrated condensation                                                !
   !---------------------------------------------------------------------------------------!
   pwavu     = 0.
   pwu_cld   = 0.
   leftu_cld = 0.
   !---------------------------------------------------------------------------------------!
   ! 3. Between the surface and the level of free convection, everything is already set up !
   !    (it was done at the grell_find_cloud_lfc subroutine) but between the LFC and the   !
   !    cloud top I will need a first guess for most variables, so I will initialise with  !
   !    the environment. Above the cloud top, no mass flux, nothing happens in terms of    !
   !    updraft, so I will set up the variables with environment values.                   !
   !---------------------------------------------------------------------------------------!
   thilu_cld((kbcon+1):mkx) = thil_cup((kbcon+1):mkx)
   tu_cld   ((kbcon+1):mkx) = t_cup   ((kbcon+1):mkx)
   qvapu_cld((kbcon+1):mkx) = qvap_cup((kbcon+1):mkx)
   qliqu_cld((kbcon+1):mkx) = qliq_cup((kbcon+1):mkx)
   qiceu_cld((kbcon+1):mkx) = qice_cup((kbcon+1):mkx)
   qtotu_cld((kbcon+1):mkx) = qtot_cup((kbcon+1):mkx)
   rhou_cld ((kbcon+1):mkx) = rho_cup ((kbcon+1):mkx)
   dbyu     ((kbcon+1):mkx) = 0.



   !---------------------------------------------------------------------------------------!
   ! 4. Between the LFC and the cloud top, we need to consider entrainment and detrainment !
   !    contributions. Also, as the air parcel moves upward, condensation will happen, and !
   !    a fraction of this condensation will fall out as rainfall. Here we will compute    !
   !    all these variables, some of them in an iterative way.                             !
   !---------------------------------------------------------------------------------------!
   do k=kbcon,ktop
      !------------------------------------------------------------------------------------!
      !    The total mixing ratio is what came from level immediately beneath us plus      !
      ! entrainment and detrainment, plus the water and ice that is about to leave the     !
      ! updraft.                                                                           !
      !------------------------------------------------------------------------------------!
      denomin     = 1.+(mentru_rate(k)-0.5*cdu(k))*dzu_cld(k)
      denomini    = 1./denomin
      qeverything = ( qtotu_cld(k-1)*(1.-.5*cdu(k)*dzu_cld(k))                             &
                     + qtot(k-1)*mentru_rate(k)*dzu_cld(k) - 0.5*leftu_cld(k-1) )          &
                    * denomini

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=38,fmt='(a)') '-----------------------------------------------------------'
      !write(unit=38,fmt='(a,1x,i5,1x,3(a,1x,f12.4,1x))')                                   &
      !   'Input values. k= ',k,'qeverything=',1000.*qeverything                            &
      !  ,'theivu_cld=',theivu_cld(k),'p_cup=',0.01*p_cup(k)
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!


      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      !    The solution of Grell (1993) equation's (A.12) is now done iteratively rather   !
      ! than using the original method. By doing this we don't need to use the approxima-  !
      ! tion for qsat, although we still assume the downdraft is saturated. We will use    !
      ! the zeroin method, which is just a combination of secant and bisection to find the !
      ! new qtotd_cld. In zero-in method, secant is the standard because it usually        !
      ! converges fast. If secant turns out to be a bad choice (because the secant is too  !
      ! flat so the new guess is too far) we use bisection instead. Likewise, the triple   !
      ! point is an obstacle to secant which may create a "bouncing" effect, never con-    !
      ! verging, so bisection becomes the standard if it doesn't converge after a few      !
      ! iterations.                                                                        !
      !    To begin with, we need three guesses, two in which the function has opposite    !
      ! signs for bisection, and a third one which will be the secant's second "chrono-    !
      ! logical" guess.                                                                    !
      !------------------------------------------------------------------------------------!
      
      !------------------------------------------------------------------------------------!
      ! a. Initialising the convergence and bisection flags.                               !
      !------------------------------------------------------------------------------------!
      converged   = .false.
      bisection   = .false.
      
      !------------------------------------------------------------------------------------!
      ! b. 1st guess outside the loop. For the 1st guess we assume no fall-out condensate. !
      !    It may turn out to be the case because the environment may be bringing too dry  !
      !    air, bringing qtotua under saturation. If that's the case, we already have a    !
      !    solution, so we don't need to iterate (funa will be zero, which means that this !
      !    is a root). Previous tests have shown that even if we shift it from this guess  !
      !    it will eventually return to qtotua. 
      !------------------------------------------------------------------------------------!
      qtotua = qeverything
      !----- Finding the equilibrium state ------------------------------------------------!
      thilu_cld(k) = thetaeiv2thil(theivu_cld(k),p_cup(k),qtotua)
      call thil2tqall(thilu_cld(k),exner_cup(k),p_cup(k),qtotua,qliqu_cld(k)               &
                     ,qiceu_cld(k),tu_cld(k),qvapu_cld(k),qsatu_cld(k))
      leftu_cld(k) = c0 * (qliqu_cld(k) + qiceu_cld(k))*dzu_cld(k)
      funa         = qtotua - qeverything + 0.5 * leftu_cld(k) * denomini

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=38,fmt='(2(a,1x,i5,1x),a,1x,l1,1x,8(a,1x,f10.4,1x),a,1x,es12.5)')        &
      !    'k=',k,'it=',-1,'bisection=',.false.,'qall=',1000.*qeverything                   &
      !   ,'qtot=',1000.*qtotua,'left=',1000.*leftu_cld(k),'qsat=',1000.*qsatu_cld(k)       &
      !   ,'qvap=',1000.*qvapu_cld(k),'qliq=',1000.*qliqu_cld(k)                            &
      !   ,'qice=',1000.*qiceu_cld(k),'temp=',tu_cld(k)-t00,'funa=',1000.*funa
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
    
      !----- We will iterate only if funa is not zero -------------------------------------!
      converged = funa == 0.
      iterif: if (.not. converged) then

         !---------------------------------------------------------------------------------!
         ! c. 3rd guess, the second for the secant method, which will be computed on top   !
         !    of the first one. qtotda will be also qtotdp.                                !
         !---------------------------------------------------------------------------------!
         qtotup = qtotua
         funp   = funa
         !------ Finding the current guess ------------------------------------------------!
         qtotuc = max(toodry,qeverything - 0.5 * leftu_cld(k) * denomini)
         thilu_cld(k) = thetaeiv2thil(theivu_cld(k),p_cup(k),qtotuc)
         call thil2tqall(thilu_cld(k),exner_cup(k),p_cup(k),qtotuc,qliqu_cld(k)            &
                        ,qiceu_cld(k),tu_cld(k),qvapu_cld(k),qsatu_cld(k))
         leftu_cld(k) = c0 * (qliqu_cld(k) + qiceu_cld(k)) * dzu_cld(k)
         func         = qtotuc - qeverything + 0.5 * leftu_cld(k) * denomini

         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write (unit=38,fmt='(2(a,1x,i5,1x),a,1x,l1,1x,8(a,1x,f10.4,1x),a,1x,es12.5)')     &
         !    'k=',k,'it=',0,'bisection=',.false.,'qall=',1000.*qeverything                 &
         !   ,'qtot=',1000.*qtotuc,'left=',1000.*leftu_cld(k),'qsat=',1000.*qsatu_cld(k)    &
         !   ,'qvap=',1000.*qvapu_cld(k),'qliq=',1000.*qliqu_cld(k)                         &
         !   ,'qice=',1000.*qiceu_cld(k),'temp=',tu_cld(k)-t00,'func=',1000.*func
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

         
         !---------------------------------------------------------------------------------!
         ! d. 3rd guess. This will seek a function evaluation with opposite sign for bi-   !
         !    section. We have already two guesses, so first we check whether they have    !
         !    opposite signs. If not, then we use the secant between funa and func to      !
         !    extrapolate to a new guess that would lead to -funa should the function be   !
         !    linear. Usually this will lead to the opposite sign at the first trial, if   !
         !    not, extrapolate further...                                                  !
         !---------------------------------------------------------------------------------!
         browsegss: if (funa * func < 0.) then
            qtotuz = qtotuc
            funz   = func
         else 
            if (abs(func-funa) < toler*qtotua) then
               delta = 100.*toler*qtotua
            else
               delta = max(abs(funa*(qtotuc-qtotua)/(func-funa)),100.*toler*qtotua)
            end if
            qtotuz    = qtotua + delta
            tubis     = tu_cld(k) !---- Using a scratch to avoid sending td_cld too far ---!
            funz      = funa
            bisection = .false.
         !----- Just to enter at least once. The 1st time qtotdz=qtotda-2*delta -----------!
            zgssloop: do it=1,maxfpo
               qtotuz = max(toodry,qtotua + real((-1)**it * (it+3)/2) * delta)

               !----- Finding this equilibrium state --------------------------------------!
               thilu_cld(k) = thetaeiv2thil(theivu_cld(k),p_cup(k),qtotuz)
               call thil2tqall(thilu_cld(k),exner_cup(k),p_cup(k),qtotuz                   &
                              ,qliqu_cld(k),qiceu_cld(k),tubis,qvapu_cld(k)                &
                              ,qsatu_cld(k))
               leftu_cld(k) = c0 * (qliqu_cld(k) + qiceu_cld(k)) * dzu_cld(k)
               funz         = qtotuz - qeverything + 0.5 * leftu_cld(k) * denomini

               bisection = funa*funz < 0.
               if (bisection) exit zgssloop
            end do zgssloop
            if (.not. bisection) then
               write (unit=*,fmt='(a)') ' No second guess for you...'
               write (unit=*,fmt='(2(a,1x,es14.7))') 'qtota=',qtotua,'funa=',funa
               write (unit=*,fmt='(2(a,1x,es14.7))') 'qtotc=',qtotuc,'func=',func
               write (unit=*,fmt='(2(a,1x,es14.7))') 'qtotz=',qtotuz,'funz=',funz
               write (unit=*,fmt='(2(a,1x,es14.7))') 'delta=',delta, 'tubis=',tubis
               call abort_run('Failed finding the second guess for bisection'              &
                             ,'grell_most_thermo_updraft','grell_cupar_updraft.f90'    )
            end if
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=38,fmt='(a,1x,i5,1x,a,1x,l1,1x,8(a,1x,f10.4,1x),a,1x,es12.5)')     &
            !    'k=',k,'it=     ½ bisection=',.true.,'qall=',1000.*qeverything             &
            !   ,'qtot=',1000.*qtotuc,'left=',1000.*leftu_cld(k),'qsat=',1000.*qsatu_cld(k) &
            !   ,'qvap=',1000.*qvapu_cld(k),'qliq=',1000.*qliqu_cld(k)                      &
            !   ,'qice=',1000.*qiceu_cld(k),'temp=',tubis-t00,'funz=',1000.*funz
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         end if browsegss
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         ! e. Now I will enter the loop. After updating, I will check the new function     !
         !    evaluation and update the (A Z) set accordingly, even if I am not using bi-  !
         !    section. Since (C P) are "current" and "previous", they will be always up-   !
         !    dated.                                                                       !
         !---------------------------------------------------------------------------------!
         bisection=.false. 
         itloop: do it=1,maxfpo
            !------------------------------------------------------------------------------!
            ! e1. Deciding whether to go with bisection or not. I should go with bisection !
            !     if the secant is dangerously small (derivative too flat, which causes    !
            !     divergence). Also if it didn't converge fast with secant, fall back to   !
            !     bisection.                                                               !
            !------------------------------------------------------------------------------!
            bisection= it > maxfpo/4 .or. abs(func-funp) < toler * qtotuc

            !------------------------------------------------------------------------------!
            ! e2. Setting the new guess. Still not sure with which method I should go, so  !
            !     establish the new guess using secant. If the guess is outside the range  !
            !     defined by the A Z pair, use bisection this time.                        !
            !------------------------------------------------------------------------------!
            if(.not. bisection) then
               qtotu_cld(k) = max(toodry,( func*qtotup - qtotuc*funp ) / (func-funp))
               !----- Checking whether this new guess represents an improvement -----------!
               bisection    = abs(qtotu_cld(k)-qtotua) > abs(qtotuz-qtotua) .or.           &
                              abs(qtotu_cld(k)-qtotuz) > abs(qtotuz-qtotua)
            end if
            if (bisection) qtotu_cld(k) = 0.5 * (qtotua + qtotuz)
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            ! e3. Finding the new function evaluation.                                     !
            !------------------------------------------------------------------------------!
            thilu_cld(k) = thetaeiv2thil(theivu_cld(k),p_cup(k),qtotu_cld(k))
            call thil2tqall(thilu_cld(k),exner_cup(k),p_cup(k),qtotu_cld(k),qliqu_cld(k)   &
                           ,qiceu_cld(k),tu_cld(k),qvapu_cld(k),qsatu_cld(k))
            leftu_cld(k) = c0 * (qliqu_cld(k) + qiceu_cld(k)) * dzu_cld(k)
            funnow       = qtotu_cld(k) - qeverything + 0.5 * leftu_cld(k) * denomini

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=38,fmt='(2(a,1x,i5,1x),a,1x,l1,1x,8(a,1x,f10.4,1x),a,1x,es12.5)')  &
            !    'k=',k,'it=',it,'bisection=',bisection,'qall=',1000.*qeverything           &
            !   ,'qtot=',1000.*qtotu_cld(k),'left=',1000.*leftu_cld(k)                      &
            !   ,'qsat=',1000.*qsatu_cld(k),'qvap=',1000.*qvapu_cld(k)                      &
            !   ,'qliq=',1000.*qliqu_cld(k),'qice=',1000.*qiceu_cld(k)                      &
            !   ,'temp=',tu_cld(k)-t00,'funnow=',1000.*funnow
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

            !------------------------------------------------------------------------------!
            ! e4. Testing for convergence, depending on the method. We may be lucky and    !
            !     hit the zero-error, depending on 
            !------------------------------------------------------------------------------!
            if (funnow == 0.) then
               converged = .true.
            elseif (bisection) then 
               converged = abs(qtotu_cld(k)-qtotua) < toler*qtotu_cld(k)
            else
               converged = abs(qtotu_cld(k)-qtotuc) < toler*qtotu_cld(k)
            end if
            !----- Found a good set, leaving... -------------------------------------------!
            if (converged) exit itloop
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            ! e5. Checking which side from the A Z pair I can update, based on funnow.     !
            !------------------------------------------------------------------------------!
            if (funnow*funa < 0.) then
               funz   = funnow
               qtotuz = qtotu_cld(k)
            else
               funa   = funnow
               qtotua = qtotu_cld(k)
            end if

            !------------------------------------------------------------------------------!
            ! e6. Updating the Previous-Current pair for the next secant attempt.          !
            !------------------------------------------------------------------------------!
            qtotup = qtotuc
            funp   = func
            qtotuc = qtotu_cld(k)
            func   = funnow


         end do itloop
      end if iterif

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=38,fmt='(a)') '-----------------------------------------------------------'
      !write(unit=38,fmt='(a)') ' '
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      if (.not. converged) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' Fall-out water finding didn''t converge!!!'
         write (unit=*,fmt='(a,1x,i5,1x,a)') ' I gave up, after',maxfpo,'iterations...'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Input values.'
         write (unit=*,fmt='(a,1x,i5)') 'Level k= ',k
         write (unit=*,fmt='(a,1x,f12.4)' ) 'qeverything       [g/kg] =',1000.*qeverything
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Upd. Thetae_iv    [   K] =',theivu_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'p_cup             [ hPa] =',0.01*p_cup(k)
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome (updraft values).'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Total  mix. ratio [g/kg] =',1000.*qtotu_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Vapour mix. ratio [g/kg] =',1000.*qvapu_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Liquid mix. ratio [g/kg] =',1000.*qliqu_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Ice    mix. ratio [g/kg] =',1000.*qiceu_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Loss   mix. ratio [g/kg] =',1000.*leftu_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Theta_il          [   K] =',thilu_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Temperature       [  °C] =',tu_cld(k)-t00
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Function          [g/kg] =',1000.*funnow
         write (unit=*,fmt='(a,1x,es12.4)') 'Error (secant)    [ ---] ='                   &
                                            ,abs(qtotu_cld(k)-qtotuc)/qtotu_cld(k) 
         write (unit=*,fmt='(a,1x,es12.4)') 'Error (bisection) [ ---] ='                   &
                                            ,abs(qtotu_cld(k)-qtotua)/qtotu_cld(k) 
         write (unit=*,fmt='(a)') '-------------------------------------------------------' 

         call abort_run('Couldn''t find the fall-out water, the zeroin method diverged...' &
                       ,'grell_most_thermo_updraft','grell_cupar_updraft.f90')
      end if

      !------------------------------------------------------------------------------------!
      !    Liquid water that leaves the cloud as rain, in kg[liq. water]/kg[air], remember-!
      ! ing that c0 is in m^-1.                                                            !
      !------------------------------------------------------------------------------------!
      pwu_cld(k) = leftu_cld(k)*etau_cld(k)
      !----- Integrate condensation -------------------------------------------------------!
      pwavu = pwavu + pwu_cld(k)
      !------ Finding density, assuming pu_cld(k) ~= p_cup(k)... --------------------------!
      rhou_cld(k) = idealdens(p_cup(k),tu_cld(k),qvapu_cld(k),qtotu_cld(k))
      !------ Finding buoyancy ------------------------------------------------------------!
      dbyu(k) = buoyancy_acc(rho_cup(k),rhou_cld(k))

   end do

   return
end subroutine grell_most_thermo_updraft
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the cloud work function associated with updrafts. Since the   !
! buoyancy acceleration was computed before, we can use it instead of the equation more    !
! often used. From Grell (1993) paper, we are using the cloud work definition from equa-   !
! tion A.38 rather than equation A.40, remembering that both are equivalent.               !
!------------------------------------------------------------------------------------------!
subroutine grell_cldwork_updraft(mkx,mgmzp,kbcon,ktop,dbyu,dzu_cld,etau_cld,aau)

   implicit none

   integer               , intent(in)  :: mkx, mgmzp ! Grid dimesnsions
   integer               , intent(in)  :: kbcon      ! Level of free convection
   integer               , intent(in)  :: ktop       ! Cloud top

   real, dimension(mgmzp), intent(in)  :: dbyu       ! Buoyancy acceleration       [  m/s²]
   real, dimension(mgmzp), intent(in)  :: dzu_cld    ! Delta-z for updrafts        [     m]
   real, dimension(mgmzp), intent(in)  :: etau_cld   ! Normalized mass flux        [   ---]
   real                  , intent(out) :: aau        ! Cloud work                  [  J/kg]

   integer                             :: k          ! Counter

   !----- Initialize cloud work to zero. --------------------------------------------------!
   aau = 0.
   
   !---------------------------------------------------------------------------------------!
   !    The cloud work is a measure of efficiency of kinetic energy generation inside the  !
   ! cloud, and therefore it is directly proportional to the mass flux and buoyancy. The   !
   ! final value should represent the cloud function for the entire cloud thus the         !
   ! integral between the LFC and cloud top.                                               !
   !---------------------------------------------------------------------------------------!
   do k=kbcon,ktop-1
      aau = aau + etau_cld(k)*dbyu(k-1) *dzu_cld(k)
   end do

   !----- Include ktop only if buoyancy is positive there ---------------------------------!
   if (dbyu(ktop-1) > 0) aau = aau + etau_cld(ktop)*dbyu(ktop-1)*dzu_cld(ktop)
   aau = max(0.,aau)

   return
end subroutine grell_cldwork_updraft
!==========================================================================================!
!==========================================================================================!
