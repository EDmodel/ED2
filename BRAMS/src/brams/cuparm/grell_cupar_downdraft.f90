!==========================================================================================!
! grell_cupar_downdraft.f90                                                                !
!                                                                                          !
!    This file contains subroutines that will calculate downdraft related stuff.           !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the level in which the downdraft originates.                  !
!------------------------------------------------------------------------------------------!
subroutine grell_find_downdraft_origin(mkx,mgmzp,k22,ktop,relheight_down,zcutdown,z_cup    &
                                      ,hes_cup,dzd_cld,ierr,kdet,jmin)

   implicit none
   
   integer, intent(in)                   :: mkx, mgmzp     ! Grid variables
   integer, intent(in)                   :: k22            ! Updrafts originating level
   integer, intent(in)                   :: ktop           ! Cloud top level
   
   !------ The next two variables define the uppermost possible level for downdrafts ------!
   real                  , intent(in)    :: relheight_down ! Max. relative to the cloud top
   real                  , intent(in)    :: zcutdown       ! Maximum absolute height;
   real, dimension(mgmzp), intent(in)    :: z_cup          ! Height @ cloud levels.
   real, dimension(mgmzp), intent(in)    :: hes_cup        ! Environ. Sat. Moist Static En.
   real, dimension(mgmzp), intent(in)    :: dzd_cld        ! Delta-z for downdrafts.

   integer               , intent(out)   :: jmin           ! Downdraft originating level
   integer               , intent(inout) :: kdet           ! Top of dndraft detrainm. layer
   integer               , intent(inout) :: ierr           ! Error flag

   real, dimension(mgmzp)                :: hed_tmp        ! Downdraft moist static energy
   real                                  :: ssbuoy         ! Integrated buoyancy term
   real                                  :: zktop          ! Top Height allowed for dndrafts.
   integer                               :: kzdown         ! Top level allowed for dndrafts
   integer                               :: k              ! Level counter

   !---------------------------------------------------------------------------------------!
   !   Initializing hed_tmp. This is temporary because it neglects entrainment and         !
   ! detrainment, so it won't be the definite moist static energy associated with          !
   ! downdrafts.                                                                           !
   !---------------------------------------------------------------------------------------!
   hed_tmp = 0.



   !----- Finding the top layer in which downdrafts are allowed to originate. -------------!
   zktop = min(relheight_down*z_cup(ktop),zcutdown)
   
   
   
   !---------------------------------------------------------------------------------------!
   !    The above definition should be sufficient to not allow the level to be above ktop  !
   !  anyway, but here we constrain it to be at least two levels apart.                    !
   !---------------------------------------------------------------------------------------!
   kzdownloop: do kzdown=1,ktop-2
      if (z_cup(kzdown) > zktop) exit kzdownloop
   end do kzdownloop



   !---------------------------------------------------------------------------------------!
   !    Between k22 and the downdraft upper bound, look for the minimum hes0_cup. This     !
   ! will be the 1st. guess for the level that downdrafts originate (aka jmin).            !
   !---------------------------------------------------------------------------------------!
   jmin=(k22-1) + minloc(hes_cup(k22:kzdown),dim=1)



   !---------------------------------------------------------------------------------------!
   !    This is going to be done iteractively until a suitable combination of jmin that    !
   ! will be associated with a layer of negative buoyancy, do I will try to find such      !
   ! combination. In case I fail, this cloud won't exist.                                  !
   !---------------------------------------------------------------------------------------!
   jminloop: do

      !------------------------------------------------------------------------------------!
      !    First thing to check: Is this jmin too low? If so, I won't allow this cloud to  !
      ! happen.                                                                            !
      !------------------------------------------------------------------------------------!
      if (jmin <= 3) then
         ierr = 4
         return
      end if

      !------------------------------------------------------------------------------------!
      !    If jmin is at the same level or below the level in which downdrafts should      !
      ! detrain all their mass, move the detrainment level to the one immediately below    !
      ! jmin.                                                                              !
      !------------------------------------------------------------------------------------!
      if (jmin <= kdet) kdet = jmin - 1
            
      !----- The downdraft starts with the saturation moist static energy -----------------!
      hed_tmp(jmin) = hes_cup(jmin)
      
      !----- Initialize the integrated buoyancy term --------------------------------------!
      ssbuoy = 0

      !----- Loop over layers beneath jmin, find the integrated buoyancy factor -----------!
      do k=jmin-1,1,-1
      
         !---------------------------------------------------------------------------------!
         !    Downdrafts are moist adiabatic at this point. Since they detrain all mass in !
         ! a single level, kdet, between jmin and kdet the downdraft moist static energy   !
         ! is conserved. Therefore, it should be the same as at jmin.                      !
         !---------------------------------------------------------------------------------!
         hed_tmp(k)=hes_cup(jmin)

         !----- Finding the integrated buoyancy term --------------------------------------!
         ssbuoy = ssbuoy + dzd_cld(k)*(hed_tmp(k)-hed_tmp(k))
         
         !---------------------------------------------------------------------------------!
         !   Check the sign of the integrated buoyancy. If it is positive, then the        !
         ! current jmin is not suitable, so I will try the level below and repeat the      !
         ! procedure.                                                                      !
         !---------------------------------------------------------------------------------!
         if (ssbuoy > 0) then
            jmin = jmin -1
            cycle jminloop
         end if
      end do
      
      !------------------------------------------------------------------------------------!
      ! If I reached this point, it means that I found a decent jmin, stick with it and    !
      ! leave the subroutine                                                               !
      !------------------------------------------------------------------------------------!
      exit jminloop

   end do jminloop

   return
end subroutine grell_find_downdraft_origin
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the normalized mass flux associated with downdrafts           !
!------------------------------------------------------------------------------------------!
subroutine grell_nms_downdraft(mkx,mgmzp,kdet,jmin,part_detr,pmass_left,mentrd_rate,cdd    &
                              ,z_cup,dzd_cld,etad_cld)
   implicit none

   integer               , intent(in)    :: mkx, mgmzp  ! Grid dimesnsions
   integer               , intent(in)    :: kdet        ! Level in which downdrafts detrain
   integer               , intent(in)    :: jmin        ! Level in which downdrafts begin

   logical               , intent(in)    :: part_detr   ! Perform partial detrainment?
   real                  , intent(in)    :: pmass_left  ! Percentage of mass left when 
                                                        !    hitting the ground.
   real, dimension(mgmzp), intent(in)    :: mentrd_rate ! Updraft entrainment rate
   real, dimension(mgmzp), intent(inout) :: cdd         ! Updraft detrainment function;
   real, dimension(mgmzp), intent(in)    :: z_cup       ! Height @ cloud levels;
   real, dimension(mgmzp), intent(in)    :: dzd_cld     ! Delta-z for downdrafts
   real, dimension(mgmzp), intent(out)   :: etad_cld    ! Normalized updraft flux

   integer                               :: k           ! Counter
   real                                  :: pmass_kept  ! 1.-pmass_left
   
   pmass_kept = 1.-pmass_left
   
   !---------------------------------------------------------------------------------------!
   ! 1. Above the level in which downdrafts begin, set the mass flux to zero.              !
   !---------------------------------------------------------------------------------------!
   etad_cld((jmin+1):mkx) = 0.
   
   
   
   !---------------------------------------------------------------------------------------!
   ! 2. At the level in which downdrafts begin, set the mass flux to one.                  !
   !---------------------------------------------------------------------------------------!
   etad_cld(jmin) = 1.
   
   
   !---------------------------------------------------------------------------------------!
   ! 3. Between the level in which downdrafts originate and the level of detrainment,      !
   !    include the entrainment effect.                                                    !
   !---------------------------------------------------------------------------------------!
   do k=jmin-1,kdet+1,-1
      etad_cld(k)=etad_cld(k+1)*(1+mentrd_rate(k)*dzd_cld(k))
   end do
   
   
   
   !---------------------------------------------------------------------------------------!
   ! 4. Between the detrainment level and the surface, consider the detrainment effect.    !
   !    Unless the input is asking to recompute the detrainment, it will be zero           !
   !---------------------------------------------------------------------------------------!
   do k=kdet,1,-1
      if (part_detr) then
         cdd(k) = mentrd_rate(k)+(1. - (pmass_kept*z_cup(k) + pmass_left*z_cup(kdet))/     &
                                    (pmass_kept*z_cup(k+1) + pmass_left*z_cup(kdet)) )/    &
                                    dzd_cld(k)
      end if
      etad_cld(k) = etad_cld(k+1) * (1. + (mentrd_rate(k) - cdd(k))*dzd_cld(k))
   end do


   return
end subroutine grell_nms_downdraft
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the moist static energy and buoyancy effect associated with   !
! downdrafts.                                                                              !
!------------------------------------------------------------------------------------------!
subroutine grell_he_downdraft(mkx,mgmzp,jmin,forbid_pos_buoy,cdd,mentrd_rate,he            &
                             ,hes_cup,dzd_cld,ierr,hed_cld,dbyd)
  implicit none

   integer               , intent(in)    :: mkx, mgmzp      ! Grid dimesnsions;
   integer               , intent(in)    :: jmin            ! Downdrafts begin here;
   integer               , intent(inout) :: ierr            ! Error flag;
   logical               , intent(in)    :: forbid_pos_buoy ! Positive column buoy. check;

   real, dimension(mgmzp), intent(in)    :: mentrd_rate     ! Downdraft entrainment rate;
   real, dimension(mgmzp), intent(in)    :: cdd             ! Downdraft detrain. function;
   real, dimension(mgmzp), intent(in)    :: he              ! MSE @ model levels;
   real, dimension(mgmzp), intent(in)    :: hes_cup         ! Saturation MSE @ cloud levels;
   real, dimension(mgmzp), intent(in)    :: dzd_cld         ! Delta-z for downdrafts;
   real, dimension(mgmzp), intent(out)   :: hed_cld         ! Downdraft moist static energy;
   real, dimension(mgmzp), intent(out)   :: dbyd            ! Downdraft buoyancy term

   integer                               :: k               ! Counter
   real                                  :: totbu           ! Total buoyancy
   


   !---------------------------------------------------------------------------------------!
   ! 1. Between jmin and the top, the moist static energy is set to be the same as the     !
   !    environment saturation moist static energy. This makes sense at jmin, as the       !
   !    downdrafts originate there so they must have the same energy as the level. Above   !
   !    jmin is just a way not to leave it empty.                                          !
   !---------------------------------------------------------------------------------------!
   hed_cld(jmin:mkx) = hes_cup(jmin:mkx)
   dbyd(jmin:mkx)     = 0.
   totbu              = 0.
   


   !---------------------------------------------------------------------------------------!
   ! 2. Below jmin, the downdraft should have the same moist static energy as the level    !
   !    right above if no entrainment or entrainment ever happened. Since they may happen  !
   !    we need to take them into account.                                                 !
   !---------------------------------------------------------------------------------------!
   do k = jmin-1,1,-1
      hed_cld(k) = (hed_cld(k+1)*(1.-.5*cdd(k)*dzd_cld(k))                                 &
                 + mentrd_rate(k)*dzd_cld(k)*he(k))                                        &
                 / (1.+mentrd_rate(k)*dzd_cld(k) - 0.5*cdd(k)*dzd_cld(k))
      dbyd(k) = hed_cld(k) - hes_cup(k)
      totbu   = totbu + dbyd(k)*dzd_cld(k)
   end do
   
   !---------------------------------------------------------------------------------------!
   !    Now the buoyancy contains the effect of entrainment and detrainment, and the layer !
   ! should contain integrated negative buoyancy. If it does not (and this is checked only !
   ! at the forced computation).                                                           !
   !---------------------------------------------------------------------------------------!
   if (totbu > 0 .and. forbid_pos_buoy) ierr = 8

   return
end subroutine grell_he_downdraft
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the moisture profile, as well as the evaporation rate       !
! associated with each level.                                                              !
!------------------------------------------------------------------------------------------!
subroutine grell_moist_downdraft(mkx,mgmzp,jmin,q,q_cup,qes_cup,gamma_cup                  &
                                ,mentrd_rate,cdd,dbyd,dzd_cld,etad_cld,qtd_cld,qd_cld      &
                                ,qld_cld,pwd_cld,pwevd)
   use rconstants, only: alvli
   implicit none
   
   integer               , intent(in)  :: mkx, mgmzp  ! Grid dimesnsions
   integer               , intent(in)  :: jmin        ! Level in which downdrafts begin

   real, dimension(mgmzp), intent(in)  :: q           ! Mixing ratio @ model levels
   real, dimension(mgmzp), intent(in)  :: q_cup       ! Mixing ratio @ cloud levels
   real, dimension(mgmzp), intent(in)  :: qes_cup     ! Sat. mixing ratio @ cloud levels
   real, dimension(mgmzp), intent(in)  :: gamma_cup   ! Gamma factor
   real, dimension(mgmzp), intent(in)  :: mentrd_rate ! Downdraft entrainment rate
   real, dimension(mgmzp), intent(in)  :: cdd         ! Downdraft detrainment function;
   real, dimension(mgmzp), intent(in)  :: dbyd        ! Buoyancy term for downdrafts
   real, dimension(mgmzp), intent(in)  :: dzd_cld     ! Delta-z for downdrafts
   real, dimension(mgmzp), intent(in)  :: etad_cld    ! Normalized mass flux 
   real, dimension(mgmzp), intent(out) :: qtd_cld     ! Downdraft total mixing ratio
   real, dimension(mgmzp), intent(out) :: qld_cld     ! Downdraft condensate mixing ratio 
   real, dimension(mgmzp), intent(out) :: qd_cld      ! Downdraft mixing ratio 
   real, dimension(mgmzp), intent(out) :: pwd_cld     ! Level-dependent evaporation
   real                  , intent(out) :: pwevd       ! Total normalized integrated evap.

   integer                             :: k           ! Counter
   real                                :: qeverything ! Total mixing ratio.
   !---------------------------------------------------------------------------------------!
   !    We will compute the evaporation in this subroutine. By definition, evaporation     !
   ! happens only if the downdraft mixing ratio is less than the saturation (then there    !
   ! will be "pressure" over condensates to return to vapour phase).                       !
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! 1. Initialize the output variables. All column variables should be equal to the       !
   !    environment above the level of origin of downdrafts. They won't be used anywhere,  !
   !    but they need to be initialized anyway.                                            !
   !---------------------------------------------------------------------------------------!
   qtd_cld((jmin+1):mkx)  = q_cup((jmin+1):mkx)
   qld_cld((jmin+1):mkx)  = 0. !---- No condensated water... ------------------------------!
   qd_cld((jmin+1):mkx)   = q_cup((jmin+1):mkx)
   pwd_cld((jmin+1):mkx)  = 0. !---- No evaporation... ------------------------------------!
   
   !---------------------------------------------------------------------------------------!
   ! 2. At the level in which downdrafts begin, the downdraft total mixing ratio must      !
   !    equals the environment saturation mixing ratio. This is done in two steps. First,  !
   !    I find the total mixing ratio, considering all effects of entrainment and          !
   !    detrainment. If this leads to a value that is less than the saturation mixing      !
   !    ratio, I evaporate rain drops until the downdraft reaches equilibrium (saturation).!
   !    Since there is no previous level at jmin, I assume the air starts with the         !
   !    environment mixing ratio and no entrainment or air from aloft.                     !
   !---------------------------------------------------------------------------------------!
   qtd_cld(jmin)  = qes_cup(jmin)
   qld_cld(jmin)  = 0.
   qd_cld(jmin)   = qes_cup(jmin)
   pwd_cld(jmin)  = max(0.,(q_cup(jmin)-qtd_cld(jmin))*etad_cld(jmin))

   !----- This is the first level in which evaporation can happen, so I just copy it. -----!
   pwevd          = pwd_cld(jmin)

   !---------------------------------------------------------------------------------------!
   ! 3. Now loop through the levels beneath jmin, using the same idea as above, but now    !
   !    considering entrainment/detrainment. In principle q should be conserved should     !
   !    entrainment and detrainment not happen.                                            !
   !---------------------------------------------------------------------------------------!
   do k=jmin-1,1,-1
      qeverything = (qd_cld(k+1)*(1-0.5*cdd(k)*dzd_cld(k))                                 &
                  + mentrd_rate(k)*dzd_cld(k)*q(k))                                        &
                  / (1+(mentrd_rate(k)-0.5*cdd(k))*dzd_cld(k))
      qd_cld(k)   = qes_cup(k) + alvli * gamma_cup(k) * dbyd(k) / (1.+gamma_cup(k))
      pwd_cld(k)  = min(0.,(qeverything - qd_cld(k))*etad_cld(k))
      qtd_cld(k)  = qeverything-pwd_cld(k)
      qld_cld(k)  = max(0.,qtd_cld(k)-qd_cld(k))
      pwevd       = pwevd +pwd_cld(k)
   end do
   
   return
end subroutine grell_moist_downdraft
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!  This subroutine computes the cloud work function associated with downdrafts             !
!------------------------------------------------------------------------------------------!
subroutine grell_cldwork_downdraft(mkx,mgmzp,jmin,t_cup,gamma_cup,dbyd,dzd_cld,etad_cld,aad)
   use rconstants, only : gocp
   implicit none

   integer               , intent(in)  :: mkx, mgmzp  ! Grid dimesnsions
   integer               , intent(in)  :: jmin        ! Downdraft origin level

   real, dimension(mgmzp), intent(in)  :: t_cup       ! Temperature @ cloud levels  [    K]
   real, dimension(mgmzp), intent(in)  :: gamma_cup   ! Gamma                       [  ---]
   real, dimension(mgmzp), intent(in)  :: dbyd        ! Buoyancy term               [ J/kg]
   real, dimension(mgmzp), intent(in)  :: dzd_cld     ! Delta-z for downdrafts      [    m]
   real, dimension(mgmzp), intent(in)  :: etad_cld    ! Normalized mass flux        [  ---]
   real                  , intent(out) :: aad         ! Downdraft work function     [ J/kg]

   integer                             :: k           ! Counter
   real                                :: buoyancy    ! Buoyancy 
   real, external                      :: virtt       ! Virtual temperature fctn    [    K]

   !----- Initialize cloud work to zero. --------------------------------------------------!
   aad = 0.
   
   !---------------------------------------------------------------------------------------!
   !    The cloud work is a measure of efficiency of kinetic energy generation inside the  !
   ! cloud, and therefore it is directly proportional to the mass flux and buoyancy. The   !
   ! final value should represent the cloud function for the entire downdraft layer thus   !
   ! the integral between the surface and the level in which downdrafts originate cloud    !
   ! base and cloud top. Note that dbyd is negative so this will give a positive function. !
   !---------------------------------------------------------------------------------------!
   do k=jmin-1,1,-1
      buoyancy = gocp * dbyd(k)/(t_cup(k)*(1+gamma_cup(k)))
      aad = aad - etad_cld(k) * buoyancy *dzd_cld(k)
   end do
   if (aad < 0) aad = 0.

   return
end subroutine grell_cldwork_downdraft
!==========================================================================================!
!==========================================================================================!
