!==========================================================================================!
! grell_cupar_updraft.f90                                                                  !
!                                                                                          !
!    This file contains subroutines that will calculate updraft related stuff.             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine finds the level of origin of updrafts. Depending on the profile, this  !
! may not be the definite level, it may change during the search for a cloud base.         !
!------------------------------------------------------------------------------------------!
subroutine grell_updraft_origin(mkx,mgmzp,iupmethod,kpbl,kbmax,z,tkeg,rcpg,he_cup,ierr,k22)

   implicit none
   integer               , intent(in)     :: mkx,mgmzp ! Grid dimensions
   integer               , intent(in)     :: iupmethod ! Method to find the level
   integer               , intent(in)     :: kpbl      ! Level of PBL top if pre-computed.
   integer               , intent(in)     :: kbmax     ! Maximum allowed level for updraft
   real, dimension(mgmzp), intent(in)     :: z         ! Height
   real, dimension(mgmzp), intent(in)     :: tkeg      ! Turbulent Kinetic Energy
   real, dimension(mgmzp), intent(in)     :: rcpg      ! Cloud mixing ratio
   real, dimension(mgmzp), intent(in)     :: he_cup    ! Moist static energy @ cloud levels
   integer               , intent(inout)  :: ierr      ! Error flag
   integer               , intent(out)    :: k22       ! Updraft origin level
   
   integer               , parameter      :: kstart=3  ! Minimum level
   !---------------------------------------------------------------------------------------!
   !    Three possibilities exist:                                                         !
   !    a. The user wants to use the PBL height (iupmethod=2), and the user is running     !
   !       Nakanishi/Niino turbulence (positive pblidx, previously computed);              !
   !    b. The user wants to use the PBL height (iupmethod=2), other turbulence was used   !
   !    (pblidx is always zero in this case). Find PBL here.                               !
   !    c. The user wants to use the maximum moist static energy as the first guess.       !
   !---------------------------------------------------------------------------------------!
   if (iupmethod == 2 .and. kpbl /= 0) then
      k22 = kpbl
   elseif (iupmethod == 2 .and. kpbl == 0) then
      call grell_find_pbl_height(mkx,mgmzp,z,tkeg,rcpg,k22)
   else
      k22 = (kstart-1) + maxloc(he_cup(kstart:kbmax),dim=1)
   end if
   
   !---------------------------------------------------------------------------------------!
   !   If the level of updraft origin is too high, then cumulus parameterization should    !
   ! not be called.                                                                        !
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
recursive subroutine grell_find_cloud_lfc(mkx,mgmzp,kbmax,cap_max,he_cup,hes_cup,p_cup    &
                                          ,k22,ierr,kbcon)

   implicit none
   integer, intent(in)                :: mkx     ! # of vertical layers
   integer, intent(in)                :: mgmzp   ! Vertical dimension
   integer, intent(in)                :: kbmax   ! Top level allowed for cloud base
   real   , intent(in)                :: cap_max ! Depth of capping inversion       [   Pa]

   real, dimension(mgmzp), intent(in) :: he_cup  ! Moist static en. @ cloud level   [ J/kg]
   real, dimension(mgmzp), intent(in) :: hes_cup ! Sat. MSE @ cloud level           [ J/kg]
   real, dimension(mgmzp), intent(in) :: p_cup   ! Pressure @ cloud level           [   Pa]

   !----- These variables may not be assigned here so use inout ---------------------------!
   integer, intent(inout)             :: k22     ! Level of origin of updrafts
   integer, intent(inout)             :: ierr    ! Error flag
   
   integer, intent(out)               :: kbcon   ! Cloud base level 

   real                               :: pcdiff  ! Pres. diff. between k22 and kbcon 
   
   !---------------------------------------------------------------------------------------!
   !    Start with the first guess: the cloud base (kbcon) is at the level in  which       !
   ! updrafts originate.                                                                   !
   !---------------------------------------------------------------------------------------!
   kbcon = k22
   
   !---------------------------------------------------------------------------------------!
   !   First step: finding the level in which the air lifted from the level that updrafts  !
   ! originate would saturate.                                                             !
   !---------------------------------------------------------------------------------------!
   do while (he_cup(k22) < hes_cup(kbcon))
      kbcon=kbcon+1
      if (kbcon > kbmax + 2) then ! Gave up... Cloud would be too high to be a cumulus.
         ierr = 3
         return
      end if
   end do
   !----- If kbcon is the next above level, I will take it and happily leave --------------!
   if (kbcon-k22 == 1) then
      return
   else
      !------------------------------------------------------------------------------------!
      !    So the cloud base is kind of far from the level in which updrafts originate.    !
      ! If there is an inversion capping, then I may not develop any convection, so I try  !
      ! an alternative approach. If the "gap" between k22 and kbcon is larger than the     !
      ! maximum depth of the inversion capping, I will try to push k22 upwards.            !
      !------------------------------------------------------------------------------------!
      pcdiff = p_cup(k22) - p_cup(kbcon) ! Pressure decreases with height...
      if (pcdiff > cap_max) then
         k22 = k22 + 1
         call grell_find_cloud_lfc(mkx,mgmzp,kbmax,cap_max,he_cup,hes_cup,p_cup,k22        &
                                   ,ierr,kbcon)
      end if
   end if
   return
end subroutine grell_find_cloud_lfc
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine finds the level of convective cloud top.                               !
!------------------------------------------------------------------------------------------!
subroutine grell_find_cloud_top(mkx,mgmzp,kbcon,he_cup,heu_cld,dbyu,ierr,ktop)

   implicit none
   integer, intent(in)                   :: mkx,mgmzp ! Grid dimensions
   integer, intent(in)                   :: kbcon     ! Cloud base level
   real, dimension(mgmzp), intent(in)    :: he_cup    ! Moist static energy @ cloud levels


   real, dimension(mgmzp), intent(inout) :: heu_cld   ! Updraft MSE
   real, dimension(mgmzp), intent(inout) :: dbyu      ! Updraft buoyancy term
   integer               , intent(inout) :: ierr      ! Error flag

   integer               , intent(out)   :: ktop      ! Cloud top level
   
   logical                               :: foundit   ! Boolean variable for top check

   !---------------------------------------------------------------------------------------!
   !    Initialize foundit                                                                 !
   !---------------------------------------------------------------------------------------!
   foundit=.false.
   
   !---------------------------------------------------------------------------------------!
   !   Loop over levels looking for the first level above the cloud top that has negative  !
   ! buoyancy.                                                                             !
   !---------------------------------------------------------------------------------------!
   ktoploop: do ktop=kbcon+1,mkx-2
      foundit=dbyu(ktop) < 0.
      if (foundit) exit ktoploop
   end do ktoploop

   if (foundit) then
      !------------------------------------------------------------------------------------!
      !    If I find the negative buoyancy, then the cloud top is the previous one. This   !
      ! is because the cloud top should be below a negative buoyancy level. Then I flush   !
      ! dbyu above the cloud top to zero.                                                  !
      !------------------------------------------------------------------------------------!
      ktop=ktop-1
      dbyu(ktop+1:mkx) = 0.
      heu_cld(ktop+1:mkx)=he_cup(ktop+1:mkx)
   else
      !------------------------------------------------------------------------------------!
      !    I won't allow those clouds that extend all the way up to the mesosphere, they   !
      ! are way too scary...                                                               !
      !------------------------------------------------------------------------------------!
      ktop=0
      dbyu=0.
      heu_cld=he_cup
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
subroutine grell_he_updraft(mkx,mgmzp,k22,kbcon,cdu,mentru_rate,he,he_cup,hes_cup,dzu_cld  &
                           ,heu_cld,dbyu)
   implicit none
   integer               , intent(in)  :: mkx, mgmzp  ! Grid dimesnsions
   integer               , intent(in)  :: k22         ! Level in which updrafts begin
   integer               , intent(in)  :: kbcon       ! Cloud base

   real, dimension(mgmzp), intent(in)  :: mentru_rate ! Updraft entrainment rate
   real, dimension(mgmzp), intent(in)  :: cdu         ! Updraft detrainment function;
   real, dimension(mgmzp), intent(in)  :: he          ! Moist static energy @ model levels;
   real, dimension(mgmzp), intent(in)  :: he_cup      ! Moist Static Energy @ cloud levels;
   real, dimension(mgmzp), intent(in)  :: hes_cup     ! Saturation MSE @ cloud levels;
   real, dimension(mgmzp), intent(in)  :: dzu_cld     ! Delta-z for updrafts;
   real, dimension(mgmzp), intent(out) :: heu_cld     ! Cloud moist static energy;
   real, dimension(mgmzp), intent(out) :: dbyu        ! Buoyancy term

   integer                             :: k           ! Counter
   real                                :: he22        ! Just a shortcut for he_cup(k22)
   !---------------------------------------------------------------------------------------!
   ! 1. Below the updraft origin, the cloud moist static energy is simply the same thing   !
   !    as the environment. Buoyancy factor is set to zero.                                !
   !---------------------------------------------------------------------------------------!
   do k=1,k22-1
      heu_cld(k) = he_cup(k)
      dbyu(k)    = 0.
   end do
   
   !---------------------------------------------------------------------------------------!
   ! 2. Between the updraft origin and the cloud base, the cloud moist static energy is    !
   !   constant because it is a dry adiabatic process. Keep the same moist static energy.  !
   !---------------------------------------------------------------------------------------!
   do k=k22,kbcon-1
      heu_cld(k) = he_cup(k22)
      dbyu(k)    = 0.
   end do

   !---------------------------------------------------------------------------------------!
   ! 3. At the cloud base, the moist static energy is still he_cup, but need to find the   !
   !    buoyancy level.                                                                    !
   !---------------------------------------------------------------------------------------!
   heu_cld(kbcon) = he_cup(k22)
   dbyu(kbcon)    = heu_cld(kbcon)-hes_cup(kbcon)

   !---------------------------------------------------------------------------------------!
   ! 4. Inside the cloud, the updraft is no longer adiabatic because there is entrainment  !
   !    and detrainment (the phase change and precipitation has not started yet).          !
   !---------------------------------------------------------------------------------------!
   do k=kbcon+1,mkx-1
      heu_cld(k) = (heu_cld(k-1)*(1.-.5*cdu(k)*dzu_cld(k))                                 &
                 + mentru_rate(k)*dzu_cld(k)*he(k-1))                                      &
                 / (1+mentru_rate(k)*dzu_cld(k) - .5*cdu(k)*dzu_cld(k))
      dbyu(k)    = heu_cld(k)-hes_cup(k)
   end do
   
   return
end subroutine grell_he_updraft
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the normalized mass flux associated with updrafts             !
!------------------------------------------------------------------------------------------!
subroutine grell_nms_updraft(mkx,mgmzp,k22,kbcon,ktop,mentru_rate,cdu,dzu_cld,etau_cld)
   implicit none

   integer               , intent(in)  :: mkx, mgmzp  ! Grid dimesnsions
   integer               , intent(in)  :: k22         ! Level in which updrafts begin
   integer               , intent(in)  :: kbcon       ! Cloud base
   integer               , intent(in)  :: ktop        ! Cloud top
   
   real, dimension(mgmzp), intent(in)  :: mentru_rate ! Updraft entrainment rate
   real, dimension(mgmzp), intent(in)  :: cdu         ! Updraft detrainment function;
   real, dimension(mgmzp), intent(in)  :: dzu_cld     ! Delta-z for clouds
   real, dimension(mgmzp), intent(out) :: etau_cld    ! Normalized updraft flux

   integer                             :: k           ! Counter

   !---------------------------------------------------------------------------------------!
   ! 1. Below the updraft origin there is no upward mass flux, set it to zero.             !
   !---------------------------------------------------------------------------------------!
   etau_cld(1:(k22-1)) = 0.
   
   
   
   !---------------------------------------------------------------------------------------!
   ! 2. Between the updraft origin and the cloud base, there is no entrainment/detrainment,!
   !    so I assume that the normalized mass flux is one.                                  !
   !---------------------------------------------------------------------------------------!
   etau_cld(k22:kbcon) = 1.



   !---------------------------------------------------------------------------------------!
   ! 3. Between the cloud base and cloud top, need to consider entrainment and detrainment !
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
!   This subroutine computes the normalized mass flux associated with updrafts             !
!------------------------------------------------------------------------------------------!
subroutine grell_moist_updraft(comp_down,mkx,mgmzp,k22,kbcon,ktop,radius,q,q_cup,qes_cup   &
                              ,gamma_cup,mentru_rate,cdu,dbyu,dzu_cld,etau_cld,qtu_cld     &
                              ,qlu_cld,qu_cld,pwu_cld,pwavu)
   use rconstants, only : alvli
   implicit none
   logical               , intent(in)  :: comp_down   ! Flag for downdraft/precipitation 
   integer               , intent(in)  :: mkx, mgmzp  ! Grid dimesnsions
   integer               , intent(in)  :: k22         ! Level in which updrafts begin
   integer               , intent(in)  :: kbcon       ! Cloud base
   integer               , intent(in)  :: ktop        ! Cloud top

   real                  , intent(in)  :: radius      ! Cloud radius (for precipitation)
   real, dimension(mgmzp), intent(in)  :: q           ! Mixing ratio @ model levels
   real, dimension(mgmzp), intent(in)  :: q_cup       ! Mixing ratio @ cloud levels
   real, dimension(mgmzp), intent(in)  :: qes_cup     ! Sat. mixing ratio @ cloud levels
   real, dimension(mgmzp), intent(in)  :: gamma_cup   ! Gamma factor
   real, dimension(mgmzp), intent(in)  :: mentru_rate ! Updraft entrainment rate
   real, dimension(mgmzp), intent(in)  :: cdu         ! Updraft detrainment function;
   real, dimension(mgmzp), intent(in)  :: dbyu        ! Buoyancy term
   real, dimension(mgmzp), intent(in)  :: dzu_cld     ! Delta-z for updrafts
   real, dimension(mgmzp), intent(in)  :: etau_cld    ! Normalized mass flux 
   real, dimension(mgmzp), intent(out) :: qtu_cld     ! Updraft total mixing ratio
   real, dimension(mgmzp), intent(out) :: qlu_cld     ! Updraft water mixing ratio 
   real, dimension(mgmzp), intent(out) :: qu_cld      ! Updraft vapour mixing ratio 
   real, dimension(mgmzp), intent(out) :: pwu_cld     ! Level-dependent condensation that
                                                      !   will fall out.
   real                  , intent(out) :: pwavu       ! Total normalized integrated cond.

   integer                             :: k           ! Counter
   real                                :: c0          ! Conversion rate (cloud -> rain)
   real                                :: qeverything ! Sum of all mixing ratios

   !---------------------------------------------------------------------------------------!
   !     First of all, I figure out whether this is a precipitating cloud or not. Current- !
   ! ly precipitating clouds are the ones that contain downdrafts. This does not need to   !
   ! be the requirement, it was just for convenience. The conversion rate from cloud to    !
   ! rain should be a function of the cloud size and wind shear, but it is only a step     !
   ! function, 0 if the radius is small, or a non-zero constant otherwise. Maybe this will !
   ! be updated in the future.                                                             !
   !---------------------------------------------------------------------------------------!
   if (comp_down) then
      c0 = 0.002
   else 
      c0 = 0.
   end if
   
   !---------------------------------------------------------------------------------------!
   !   Initialize integrated condensation                                                  !
   !---------------------------------------------------------------------------------------!
   pwavu = 0.

   !---------------------------------------------------------------------------------------!
   ! 1. Between the surface and the updraft origin, nothing happens in terms of updraft, so!
   !    I make condensation to be zero and the total mixing ratio to be the environment    !
   !---------------------------------------------------------------------------------------!
   qu_cld(1:(k22-1))  = q_cup(1:(k22-1))
   qtu_cld(1:(k22-1)) = q_cup(1:(k22-1))
   qlu_cld(1:(k22-1)) = 0.
   pwu_cld(1:(k22-1)) = 0.



   !---------------------------------------------------------------------------------------!
   ! 2. Between the level in which updrafts begin and the cloud base, no condensation      !
   !    happens. Mixing ratio is a conserved property should entrainment and detrainment   !
   !    not happen. Since it is under the saturated layer, no condensation happens.        !
   !---------------------------------------------------------------------------------------!
   qu_cld(k22:(kbcon-1))  = q_cup(k22)
   qtu_cld(k22:(kbcon-1)) = q_cup(k22)
   qlu_cld(k22:(kbcon-1)) = 0.
   pwu_cld(k22:(kbcon-1)) = 0.
   


   !---------------------------------------------------------------------------------------!
   ! 3. Between the cloud base and cloud top, need to consider entrainment and detrainment !
   !    contributions. Also, as the air parcel moves upward, condensation will happen, and !
   !    a fraction of this condensation will fall out as rainfall. Here we will compute    !
   !    all these variables.                                                               !
   !---------------------------------------------------------------------------------------!
   do k=kbcon,ktop
      !------------------------------------------------------------------------------------!
      !    The total mixing ratio is what came from level immediately beneath us plus      !
      ! entrainment and detrainment.                                                       !
      !------------------------------------------------------------------------------------!
      qeverything = ( qtu_cld(k-1)*(1.-.5*cdu(k)*dzu_cld(k))                               &
                     + q(k-1)*mentru_rate(k)*dzu_cld(k)     )                              &
                  / (1.+(mentru_rate(k)-0.5*cdu(k))*dzu_cld(k))
      
      !------------------------------------------------------------------------------------!
      !    This is the updraft saturation mixing ratio for this level                      !
      !------------------------------------------------------------------------------------!
      qu_cld(k) = qes_cup(k) + alvli * gamma_cup(k)*dbyu(k) / (1. + gamma_cup(k))
      
      
      
      !------------------------------------------------------------------------------------!
      !    Liquid water that remains in the cloud after rainout.                           !
      !------------------------------------------------------------------------------------!
      qlu_cld(k) = max(0.,(qeverything-qu_cld(k))/(1. + c0*dzu_cld(k)))
      
      
      !------------------------------------------------------------------------------------!
      !    Liquid water that leaves the cloud as rain, in kg[liq. water]/kg[air], remember-!
      ! ing that c0 is in m^-1.                                                            !
      !------------------------------------------------------------------------------------!
      pwu_cld(k) = c0 * dzu_cld(k)*qlu_cld(k)*etau_cld(k)
      !----- Integrate condensation -------------------------------------------------------!
      pwavu = pwavu + pwu_cld(k)
      
      !------------------------------------------------------------------------------------!
      !    Setting total updraft mixing ratio at this level                                !
      !------------------------------------------------------------------------------------!
      qtu_cld(k) = qu_cld(k) + qlu_cld(k)
      
   end do


   !---------------------------------------------------------------------------------------!
   ! 4. Above the cloud top, no mass flux, nothing happens in terms of updraft, setting    !
   !    environment values for qtu_cld and no condensation otherwise                       !
   !---------------------------------------------------------------------------------------!
   qu_cld((ktop+1):mkx) = q_cup((ktop+1):mkx)
   qlu_cld((ktop+1):mkx) = 0.
   pwu_cld((ktop+1):mkx) = 0.
   qtu_cld((ktop+1):mkx) = q_cup((ktop+1):mkx)
   
   return
end subroutine grell_moist_updraft
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the cloud work function associated with updrafts              !
!------------------------------------------------------------------------------------------!
subroutine grell_cldwork_updraft(mkx,mgmzp,kbcon,ktop,t_cup,gamma_cup,dbyu,dzu_cld         &
                                ,etau_cld,aau)
   use rconstants, only : gocp
   implicit none

   integer               , intent(in)  :: mkx, mgmzp  ! Grid dimesnsions
   integer               , intent(in)  :: kbcon       ! Cloud base
   integer               , intent(in)  :: ktop        ! Cloud top

   real, dimension(mgmzp), intent(in)  :: t_cup       ! Height @ cloud levels;
   real, dimension(mgmzp), intent(in)  :: gamma_cup   ! Gamma factor
   real, dimension(mgmzp), intent(in)  :: dbyu        ! Buoyancy term
   real, dimension(mgmzp), intent(in)  :: dzu_cld     ! Delta-z for updrafts
   real, dimension(mgmzp), intent(in)  :: etau_cld    ! Normalized mass flux 
   real                  , intent(out) :: aau         ! Total normalized integrated cond.

   integer                             :: k           ! Counter

   !----- Initialize cloud work to zero. --------------------------------------------------!
   aau = 0.
   
   !---------------------------------------------------------------------------------------!
   !    The cloud work is a measure of efficiency of kinetic energy generation inside the  !
   ! cloud, and therefore it is directly proportional to the mass flux and buoyancy. The   !
   ! final value should represent the cloud function for the entire cloud thus the         !
   ! integral between the cloud base and cloud top.                                        !
   !---------------------------------------------------------------------------------------!
   do k=kbcon,ktop-1
      aau = aau + gocp *etau_cld(k)*dbyu(k-1) *dzu_cld(k) / (t_cup(k)*(1 + gamma_cup(k)))
   end do

   !----- Include ktop only if buoyancy is positive there ---------------------------------!
   if (dbyu(ktop-1) > 0) then
      aau = aau + gocp*etau_cld(ktop)*dbyu(ktop-1)*dzu_cld(k)                              &
                / (t_cup(ktop)*(1 + gamma_cup(ktop)))
   end if
   aau = max(0.,aau)

   return
end subroutine grell_cldwork_updraft
!==========================================================================================!
!==========================================================================================!
