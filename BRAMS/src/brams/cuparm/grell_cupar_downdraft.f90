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
subroutine grell_find_downdraft_origin(mkx,mgmzp,klou,ktop,relheight_down,zcutdown,z_cup   &
                                      ,theivs_cup,dzd_cld,ierr,kdet,klod)

   implicit none
   
   integer, intent(in)                   :: mkx, mgmzp     ! Grid variables
   integer, intent(in)                   :: klou           ! Level of origin of updrafts
   integer, intent(in)                   :: ktop           ! Cloud top level
   
   !------ The next two variables define the uppermost possible level for downdrafts ------!
   real                  , intent(in)    :: relheight_down ! Max. relative to the cloud top
   real                  , intent(in)    :: zcutdown       ! Maximum absolute height;
   real, dimension(mgmzp), intent(in)    :: z_cup          ! Height @ cloud levels.
   real, dimension(mgmzp), intent(in)    :: theivs_cup     ! Sat. thetae_iv        
   real, dimension(mgmzp), intent(in)    :: dzd_cld        ! Delta-z for downdrafts.

   integer               , intent(out)   :: klod           ! Level of origin of downdrafts
   integer               , intent(inout) :: kdet           ! Top of dndraft detrainm. layer
   integer               , intent(inout) :: ierr           ! Error flag

   real, dimension(mgmzp)                :: theivd_tmp     ! Temporary downdraft thetae_iv
   real                                  :: ssbuoy         ! Integrated buoyancy term
   real                                  :: zktop          ! Top Height allowed for dndfts.
   integer                               :: kzdown         ! Top level allowed for dndrafts
   integer                               :: k              ! Level counter

   !---------------------------------------------------------------------------------------!
   !   Initializing theivd_tmp. This is temporary because it neglects entrainment and      !
   ! detrainment, so it won't be the definite ice-vapour equivalent potential temperature  !
   ! associated with downdrafts.                                                           !
   !---------------------------------------------------------------------------------------!
   theivd_tmp = 0.



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
   !    If klou is above kzdown, then this cloud doesn't make sense.  Assign error flag    !
   ! and leave the subroutine.                                                             !
   !---------------------------------------------------------------------------------------!
   if (kzdown < klou) then
      ierr = 12
      return
   end if


   !---------------------------------------------------------------------------------------!
   !    Between klou and the downdraft upper bound, look for the minimum theivs_cup.  This !
   ! will be the 1st. guess for the level that downdrafts originate (klod).                !
   !---------------------------------------------------------------------------------------!
   klod=(klou-1) + minloc(theivs_cup(klou:kzdown),dim=1)



   !---------------------------------------------------------------------------------------!
   !    This is going to be done iteractively until a suitable klod that will be associat- !
   ! ed with a layer of negative buoyancy, do I will try to find such  combination. In     !
   ! I can't dind one, this cloud won't exist.                                             !
   !---------------------------------------------------------------------------------------!
   klodloop: do

      !------------------------------------------------------------------------------------!
      !    First thing to check: Is this klod too low? If so, I won't allow this cloud to  !
      ! happen.                                                                            !
      !------------------------------------------------------------------------------------!
      if (klod <= 3) then
         ierr = 4
         return
      end if

      !------------------------------------------------------------------------------------!
      !    If klod is at the same level or below the level in which downdrafts should      !
      ! start detraining all their mass, move the detrainment level to the one immediately !
      ! below klod.                                                                        !
      !------------------------------------------------------------------------------------!
      if (klod <= kdet) kdet = klod - 1
            
      !----- The downdraft starts with the saturation moist static energy -----------------!
      theivd_tmp(klod) = theivs_cup(klod)
      
      !----- Initialize the integrated buoyancy term --------------------------------------!
      ssbuoy = 0

      !----- Loop over layers beneath klod, find the integrated buoyancy factor -----------!
      do k=klod-1,1,-1
      
         !---------------------------------------------------------------------------------!
         !    Downdrafts are moist adiabatic at this point. Since they detrain all mass    !
         ! below kdet, between klod and kdet the downdraft moist static energy is conserv- !
         ! ed. Therefore, it should be the same as at klod.                                !
         !---------------------------------------------------------------------------------!
         theivd_tmp(k)=theivs_cup(klod)

         !----- Finding the integrated buoyancy term --------------------------------------!
         ssbuoy = ssbuoy + dzd_cld(k)*(theivd_tmp(k)-theivs_cup(k))
         
         !---------------------------------------------------------------------------------!
         !   Check the sign of the integrated buoyancy. If it is positive, then the        !
         ! current klod is not suitable, so I will try the level below and repeat the      !
         ! procedure.                                                                      !
         !---------------------------------------------------------------------------------!
         if (ssbuoy > 0) then
            klod = klod -1
            cycle klodloop
         end if
      end do
      
      !------------------------------------------------------------------------------------!
      ! If I reached this point, it means that I found a decent klod, stick with it and    !
      ! leave the subroutine                                                               !
      !------------------------------------------------------------------------------------!
      exit klodloop

   end do klodloop

   return
end subroutine grell_find_downdraft_origin
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the normalized mass flux associated with downdrafts           !
!------------------------------------------------------------------------------------------!
subroutine grell_nms_downdraft(mkx,mgmzp,kdet,klod,mentrd_rate,cdd,z_cup,dzd_cld,etad_cld)
   implicit none

   integer               , intent(in)    :: mkx, mgmzp  ! Grid dimesnsions
   integer               , intent(in)    :: kdet        ! Level in which downdrafts detrain
   integer               , intent(in)    :: klod        ! Level in which downdrafts begin

   real, dimension(mgmzp), intent(in)    :: mentrd_rate ! Dndraft entrainment rate
   real, dimension(mgmzp), intent(in)    :: cdd         ! Dndraft detrainment function;
   real, dimension(mgmzp), intent(in)    :: z_cup       ! Height @ cloud levels;
   real, dimension(mgmzp), intent(in)    :: dzd_cld     ! Delta-z for downdrafts
   real, dimension(mgmzp), intent(inout) :: etad_cld    ! Normalized updraft flux

   integer                               :: k           ! Counter
   
   !---------------------------------------------------------------------------------------!
   ! 1. Above the level in which downdrafts begin, set the mass flux to zero.              !
   !---------------------------------------------------------------------------------------!
   etad_cld((klod+1):mkx) = 0.
   
   
   
   !---------------------------------------------------------------------------------------!
   ! 2. At the level in which downdrafts begin, set the mass flux to one.                  !
   !---------------------------------------------------------------------------------------!
   etad_cld(klod) = 1.
   
   
   !---------------------------------------------------------------------------------------!
   ! 3. Between the level in which downdrafts originate and the level of detrainment,      !
   !    include the entrainment effect.                                                    !
   !---------------------------------------------------------------------------------------!
   do k=klod-1,kdet+1,-1
      etad_cld(k)=etad_cld(k+1)*(1+mentrd_rate(k)*dzd_cld(k))
   end do
   
   
   
   !---------------------------------------------------------------------------------------!
   ! 4. Between the detrainment level and the surface, consider the detrainment effect.    !
   !---------------------------------------------------------------------------------------!
   do k=kdet,1,-1
      etad_cld(k) = etad_cld(k+1) * (1. + (mentrd_rate(k) - cdd(k))*dzd_cld(k))
   end do


   return
end subroutine grell_nms_downdraft
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine computes the downdraft ice-vapour equivalent potential temperature     !
! and buoyancy effect associated with downdrafts.                                          !
!------------------------------------------------------------------------------------------!
subroutine grell_theiv_downdraft(mkx,mgmzp,klod,cdd,mentrd_rate,theiv,theiv_cup,theivs_cup &
                                ,dzd_cld,ierr,theivd_cld)
  implicit none

   integer               , intent(in)    :: mkx, mgmzp  ! Grid dimesnsions;
   integer               , intent(in)    :: klod        ! Downdrafts begin here;
   integer               , intent(inout) :: ierr        ! Error flag;

   real, dimension(mgmzp), intent(in)    :: mentrd_rate ! Entrainment rate;
   real, dimension(mgmzp), intent(in)    :: cdd         ! Detrainment function;
   real, dimension(mgmzp), intent(in)    :: theiv       ! Thetae_iv @ model levels;
   real, dimension(mgmzp), intent(in)    :: theiv_cup   ! Thetae_iv @ cloud levels;
   real, dimension(mgmzp), intent(in)    :: theivs_cup  ! Sat. thetae_iv @ cloud levels;
   real, dimension(mgmzp), intent(in)    :: dzd_cld     ! Delta-z for downdrafts;
   real, dimension(mgmzp), intent(inout) :: theivd_cld  ! Downdraft thetae_iv;

   integer                               :: k           ! Counter
   


   !---------------------------------------------------------------------------------------!
   ! 1. Between klod and the top, the moist static energy is set to be the same as the     !
   !    environment saturation moist static energy. This makes sense at klod, as the       !
   !    downdrafts originate there so they must have the same energy as the level. Above   !
   !    klod is just a way not to leave it empty.                                          !
   !---------------------------------------------------------------------------------------!
   theivd_cld(klod:mkx) = theiv_cup(klod:mkx)
   


   !---------------------------------------------------------------------------------------!
   ! 2. Below klod, the downdraft would have the same moist static energy as the level     !
   !    right above if no entrainment or entrainment ever happened. Since they may happen  !
   !    we need to take them into account.                                                 !
   !---------------------------------------------------------------------------------------!
   do k = klod-1,1,-1
      theivd_cld(k) = (theivd_cld(k+1)*(1.-0.5*cdd(k)*dzd_cld(k))                          &
                    + mentrd_rate(k)*dzd_cld(k)*theiv(k))                                  &
                    / (1.+(mentrd_rate(k)- 0.5*cdd(k))*dzd_cld(k))
   end do
   
   return
end subroutine grell_theiv_downdraft
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes the moisture profile, as well as the evaporation rate       !
! associated with each level.                                                              !
!------------------------------------------------------------------------------------------!
subroutine grell_most_thermo_downdraft(mkx,mgmzp,klod,qtot,co2,mentrd_rate,cdd,p_cup       &
                                      ,exner_cup,thil_cup,t_cup,qtot_cup,qvap_cup,qliq_cup &
                                      ,qice_cup,qsat_cup,co2_cup,rho_cup,pwav,theivd_cld   &
                                      ,etad_cld,dzd_cld,thild_cld,td_cld,qtotd_cld         &
                                      ,qvapd_cld,qliqd_cld,qiced_cld,qsatd_cld,co2d_cld    &
                                      ,rhod_cld,dbyd,pwd_cld,pwev,ierr)
   use rconstants, only: epi,rdry,t00,cpi,toodry
   use therm_lib , only: thetaeiv2thil,toler,maxfpo,idealdens
   implicit none
   
   integer               , intent(in)    :: mkx, mgmzp  ! Grid dimesnsions
   integer               , intent(in)    :: klod        ! Level in which downdrafts begin

   !----- Variables at model levels -------------------------------------------------------!
   real, dimension(mgmzp), intent(in)    :: qtot        ! Total mixing ratio       [ kg/kg]
   real, dimension(mgmzp), intent(in)    :: co2         ! CO2 mixing ratio         [   ppm]
   !----- Downdraft mass exchange rates ---------------------------------------------------!
   real, dimension(mgmzp), intent(in)    :: mentrd_rate ! Entrainment rate;        [   1/m]
   real, dimension(mgmzp), intent(in)    :: cdd         ! Detrainment function;    [   1/m]
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
   real, dimension(mgmzp), intent(in)    :: co2_cup     ! CO2 mixing ratio         [   ppm]
   real, dimension(mgmzp), intent(in)    :: rho_cup     ! Density                  [ kg/m³]
   real                  , intent(in)    :: pwav        ! Total available water    [ kg/kg]
   !----- Input variables at downdraft ----------------------------------------------------!
   real, dimension(mgmzp), intent(in)    :: theivd_cld  ! Thetae_iv                [     K]
   real, dimension(mgmzp), intent(in)    :: etad_cld    ! Normalized mass flux     [   ---]
   real, dimension(mgmzp), intent(in)    :: dzd_cld     ! Layer thickness          [     m]
   !----- Transit variables, which will be changed between the surface and klod -----------!
   real, dimension(mgmzp), intent(inout) :: thild_cld   ! Theta_il                 [     K]
   real, dimension(mgmzp), intent(inout) :: td_cld      ! Temperature              [     K]
   real, dimension(mgmzp), intent(inout) :: qtotd_cld   ! Total mixing ratio       [ kg/kg]
   real, dimension(mgmzp), intent(inout) :: qvapd_cld   ! Vapour mixing ratio      [ kg/kg]
   real, dimension(mgmzp), intent(inout) :: qliqd_cld   ! Liquid mixing ratio      [ kg/kg]
   real, dimension(mgmzp), intent(inout) :: qiced_cld   ! Ice mixing ratio         [ kg/kg]
   real, dimension(mgmzp), intent(inout) :: qsatd_cld   ! Sat. mixing ratio        [ kg/kg]
   real, dimension(mgmzp), intent(inout) :: co2d_cld    ! CO2 mixing ratio         [   ppm]
   real, dimension(mgmzp), intent(inout) :: rhod_cld    ! Density                  [ kg/m³]
   real, dimension(mgmzp), intent(inout) :: dbyd        ! Buoyancy acceleration    [  m/s²]
   !----- Output variables ----------------------------------------------------------------!
   real, dimension(mgmzp), intent(inout) :: pwd_cld     ! Normal. evap. flux       [ kg/kg]
   real                  , intent(out)   :: pwev        ! Total evaporation flux   [ kg/kg]
   !----- Transit variable ----------------------------------------------------------------!
   integer               , intent(inout) :: ierr        ! Error flag
   !----- Local variables -----------------------------------------------------------------!
   integer                :: k              ! Counter                              [  ----]
   integer                :: it             ! Iteration counter                    [  ----]
   logical                :: converged      ! Flag to test convergence             [   T|F]
   logical                :: bisection      ! Flag to use bisection                [   T|F]
   real                   :: qtotd_0_evap   ! Mixing ratio without evaporation     [ kg/kg]
   real                   :: qtotda, qtotdz ! Aux. vars for bisection iteration    [ kg/kg]
   real                   :: qtotdc, qtotdp ! Aux. vars for secant iteration       [ kg/kg]
   real                   :: funa, funz     ! Function evaluation for bisection    [ kg/kg]
   real                   :: func, funp     ! Function evaluation for secant       [ kg/kg]
   real                   :: funnow         ! Function at this iteration           [ kg/kg]
   real                   :: bytot          ! Integrated buoyancy                  [  J/kg]
   real                   :: denomin        ! Denominator, just to clean the eqn.  [  ----]
   real                   :: denomini       ! 1./denominator                       [  ----]
   real, dimension(mgmzp) :: evapd_cld      ! Evaporated mix. ratio                [ kg/kg]
   real                   :: tdbis          ! Scratch var. for temperature         [     K]
   real                   :: delta          ! Aux. var. for bisection 2nd guess    [ kg/kg]
   !----- External functions --------------------------------------------------------------!
   real, external         :: buoyancy_acc   ! Buoyancy acceleration funtion.
   !----- Constants -----------------------------------------------------------------------!
   real, parameter        :: toowet=0.030   ! Max. mix. ratio for initial guesses  [ kg/kg]
   !---------------------------------------------------------------------------------------!




   !---------------------------------------------------------------------------------------!
   !    We will compute the evaporation in this subroutine. By definition, evaporation     !
   ! happens only if the downdraft mixing ratio is less than the saturation, because in    !
   ! this case the liquid/ice phase will not be in equilibrium.                            !
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   ! 1. Initialize the output variables. All column variables should be equal to the       !
   !    environment above the level of origin of downdrafts. They won't be used anywhere,  !
   !    but they need to be initialized anyway. Within the layer in which downdrafts       !
   !    occur, we need a first guess for most variables, so we just copy the entire        !
   !    vectors.                                                                           !
   !---------------------------------------------------------------------------------------!
   thild_cld(1:mkx)  = thil_cup(1:mkx)
   td_cld   (1:mkx)  = t_cup   (1:mkx)
   qtotd_cld(1:mkx)  = qtot_cup(1:mkx)
   qvapd_cld(1:mkx)  = qvap_cup(1:mkx)
   qliqd_cld(1:mkx)  = qliq_cup(1:mkx)
   qiced_cld(1:mkx)  = qice_cup(1:mkx)
   co2d_cld (1:mkx)  = co2_cup (1:mkx)
   rhod_cld (1:mkx)  = rho_cup (1:mkx)
   evapd_cld(1:mkx)  = 0.
   dbyd     (1:mkx)  = 0.
   pwd_cld  (1:mkx)  = 0. !---- No evaporation above klod. --------------------------------!
   pwev              = 0.
   bytot             = 0. !---- This will have the column integrated buoyancy. ------------!


   !---------------------------------------------------------------------------------------!
   ! 2. Now loop through the levels beneath klod, using the same idea as above, but now    !
   !    considering entrainment/detrainment. In principle q should be conserved should     !
   !    entrainment and detrainment not happen.                                            !
   !---------------------------------------------------------------------------------------!
   do k=klod,1,-1
      !------------------------------------------------------------------------------------!
      !    This is the mixing ratio the downdraft would have did evaporation not happen.   !
      !------------------------------------------------------------------------------------!
      denomin      = 1.+(mentrd_rate(k)-0.5*cdd(k))*dzd_cld(k)
      denomini     = 1./denomin
      qtotd_0_evap = (qtotd_cld(k+1)*(1.-0.5*cdd(k)*dzd_cld(k))                            &
                   + mentrd_rate(k)*dzd_cld(k)*qtot(k) -0.5*evapd_cld(k+1))                &
                   * denomini

      !----- CO2 is assumed to be an inert gas (i.e., no sources or sinks). ---------------!
      co2d_cld(k) = (co2d_cld(k+1)*(1.-0.5*cdd(k)*dzd_cld(k))                              &
                  + mentrd_rate(k)*dzd_cld(k)*co2(k) )                                     &
                  * denomini

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=28,fmt='(a)') '-----------------------------------------------------------'
      !write(unit=28,fmt='(a,1x,i5,1x,3(a,1x,f12.4,1x))')                                   &
      !   'Input values. k= ',k,'qtotd_0_evap=',1000.*qtotd_0_evap                          &
      !  ,'theivu_cld=',theivd_cld(k),'p_cup=',0.01*p_cup(k)
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!



      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      !    The solution of Grell (1993) equation's (A.13) is now done iteratively rather   !
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
      ! b. 1st guess outside the loop. For the 1st guess we assume no evaporation state.   !
      !    It may turn out to be the case because the environment entrainment may bring    !
      !    saturated air. If that's the case we have the solution and we don't need to     !
      !    iterate. (funa will be zero so that's the root we are looking for). Previous    !
      !    tests have shown that even if we shift it from this state for the first guess,  !
      !    it will eventually return to qtotda.                                            !
      !------------------------------------------------------------------------------------!
      qtotda       = qtotd_0_evap
      !----- Finding the equilibrium state ------------------------------------------------!
      thild_cld(k) = thetaeiv2thil(theivd_cld(k),p_cup(k),qtotda)
      call thil2tqall(thild_cld(k),exner_cup(k),p_cup(k),qtotda,qliqd_cld(k)               &
                     ,qiced_cld(k),td_cld(k),qvapd_cld(k),qsatd_cld(k))
      evapd_cld(k) = min(0.,2. * (qtotd_0_evap - qsatd_cld(k)) * denomin)
      funa         = qtotda - qtotd_0_evap + 0.5 * evapd_cld(k) * denomini

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=28,fmt='(2(a,1x,i5,1x),a,1x,l1,1x,8(a,1x,f10.4,1x),a,1x,es12.5)')        &
      !     'k=',k,'it=',-1,'bisection=',.false.,'q000=',1000.*qtotd_0_evap                 &
      !    ,'qtot=',1000.*qtotda,'evap=',1000.*evapd_cld(k),'qsat=',1000.*qsatd_cld(k)      &
      !    ,'qvap=',1000.*qvapd_cld(k),'qliq=',1000.*qliqd_cld(k)                           &
      !    ,'qice=',1000.*qiced_cld(k),'temp=',td_cld(k)-t00,'funa=',funa
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      !----- We will iterate only if funa is not zero. ------------------------------------!
      converged = funa == 0.
      iterif: if (.not. converged) then

         !---------------------------------------------------------------------------------!
         ! c. 2nd guess, the second for the secant method, which will be computed on top   !
         !    of the first one. qtotda will be also qtotdp.However, because the first      !
         !    guess started assuming that no evaporation have occurred, the error is like- !
         !    ly to be large, which may overestimate evapd_cld and put the second guess    !
         !    way off. Should this happen, we impose a maximum qtotd of 30 g/kg, which is  !
         !    significantly above normal values, but not too large to the point that could !
         !    cause the model to crash.                                                    !
         !---------------------------------------------------------------------------------!
         qtotdp = qtotda
         funp   = funa
         !---------------------------------------------------------------------------------!
         !     Finding the current guess. As a 2nd guess, qtotdc can be way off, in which  !
         ! case we put a lid of 30. g/kg, which is a faily large number and should never   !
         ! be below the maximum under normal conditions.                                   !
         !---------------------------------------------------------------------------------!
         qtotdc = min(max(toodry,qtotd_0_evap - 0.5 * evapd_cld(k) * denomini),toowet)
         thild_cld(k) = thetaeiv2thil(theivd_cld(k),p_cup(k),qtotdc)
         call thil2tqall(thild_cld(k),exner_cup(k),p_cup(k),qtotdc,qliqd_cld(k)            &
                        ,qiced_cld(k),td_cld(k),qvapd_cld(k),qsatd_cld(k))
         evapd_cld(k) = min(0.,2. * (qtotd_0_evap - qsatd_cld(k)) * denomin)
         func         = qtotdc - qtotd_0_evap + 0.5 * evapd_cld(k) * denomini
         !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
         !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
         !write (unit=28,fmt='(2(a,1x,i5,1x),a,1x,l1,1x,8(a,1x,f10.4,1x),a,1x,es12.5)')     &
         !     'k=',k,'it=',0,'bisection=',.false.,'q000=',1000.*qtotd_0_evap               &
         !    ,'qtot=',1000.*qtotdc,'evap=',1000.*evapd_cld(k),'qsat=',1000.*qsatd_cld(k)   &
         !    ,'qvap=',1000.*qvapd_cld(k),'qliq=',1000.*qliqd_cld(k)                        &
         !    ,'qice=',1000.*qiced_cld(k),'temp=',td_cld(k)-t00,'func=',func
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
            qtotdz = qtotdc
            funz   = func
         else 
            if (abs(func-funa) < toler*qtotda) then
               delta = 100.*toler*qtotda
            else
               delta = max(abs(funa*(qtotdc-qtotda)/(func-funa)),100.*toler*qtotda)
            end if
            qtotdz    = qtotda + delta
            tdbis     = td_cld(k) !---- Using a scratch to avoid sending td_cld too far ---!
            funz      = funa
            bisection = .false.
         !----- Just to enter at least once. The 1st time qtotdz=qtotda-2*delta -----------!
            zgssloop: do it=1,maxfpo
               qtotdz = min(max(toodry,qtotda + real((-1)**it * (it+3)/2) * delta),toowet)
               !----- Finding this equilibrium state --------------------------------------!
               thild_cld(k) = thetaeiv2thil(theivd_cld(k),p_cup(k),qtotdz)
               call thil2tqall(thild_cld(k),exner_cup(k),p_cup(k),qtotdz                   &
                              ,qliqd_cld(k),qiced_cld(k),tdbis,qvapd_cld(k)                &
                              ,qsatd_cld(k))
               evapd_cld(k) = min(0., 2. * (qtotd_0_evap - qsatd_cld(k)) * denomin )
               funz         = qtotdz - qtotd_0_evap + 0.5 * evapd_cld(k) * denomini

               bisection = funa*funz < 0.
               if (bisection) exit zgssloop
            end do zgssloop
            if (.not. bisection) then
               write (unit=*,fmt='(a)') ' No second guess for you...'
               write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'qtota=',qtotda,'funa=',funa
               write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'qtotc=',qtotdc,'func=',func
               write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'qtotz=',qtotdz,'funz=',funz
               write (unit=*,fmt='(2(a,1x,es14.7,1x))') 'delta=',delta, 'tdbis=',tdbis
               call abort_run('Failed finding the second guess for bisection'              &
                             ,'grell_most_thermo_downdraft','grell_cupar_downdraft.f90'    )
            end if
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=28,fmt='(a,1x,i5,1x,a,1x,l1,1x,8(a,1x,f10.4,1x),a,1x,es12.5)')     &
            !     'k=',k,'it=     ½ bisection=',.true.,'q000=',1000.*qtotd_0_evap           &
            !    ,'qtot=',1000.*qtotdc,'evap=',1000.*evapd_cld(k)                           &
            !    ,'qsat=',1000.*qsatd_cld(k),'qvap=',1000.*qvapd_cld(k)                     &
            !    ,'qliq=',1000.*qliqd_cld(k),'qice=',1000.*qiced_cld(k)                     &
            !    ,'temp=',tdbis-t00,'funz=',funz
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
            bisection= it > maxfpo/4 .or. abs(func-funp) < toler * qtotdc

            !------------------------------------------------------------------------------!
            ! e2. Setting the new guess. Still not sure with which method I should go, so  !
            !     establish the new guess using secant. If the guess is outside the range  !
            !     defined by the A Z pair, use bisection this time.                        !
            !------------------------------------------------------------------------------!
            if(.not. bisection) then
               qtotd_cld(k) = max(toodry,( func*qtotdp - qtotdc*funp ) / (func-funp))
               !----- Checking whether this new guess represents an improvement -----------!
               bisection    = abs(qtotd_cld(k)-qtotda) > abs(qtotdz-qtotda) .or.           &
                              abs(qtotd_cld(k)-qtotdz) > abs(qtotdz-qtotda)
            end if
            if (bisection) qtotd_cld(k) = 0.5 * (qtotda + qtotdz)
            !------------------------------------------------------------------------------!

            !------------------------------------------------------------------------------!
            ! e3. Finding the new function evaluation.                                     !
            !------------------------------------------------------------------------------!
            thild_cld(k) = thetaeiv2thil(theivd_cld(k),p_cup(k),qtotd_cld(k))
            call thil2tqall(thild_cld(k),exner_cup(k),p_cup(k),qtotd_cld(k),qliqd_cld(k)   &
                           ,qiced_cld(k),td_cld(k),qvapd_cld(k),qsatd_cld(k))
            evapd_cld(k) = min(0., 2. * (qtotd_0_evap - qsatd_cld(k)) * denomin )
            funnow       = qtotd_cld(k) - qtotd_0_evap + 0.5 * evapd_cld(k) * denomini

            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
            !write (unit=28,fmt='(2(a,1x,i5,1x),a,1x,l1,1x,8(a,1x,f10.4,1x),a,1x,es12.5)')  &
            !    'k=',k,'it=',it,'bisection=',.false.,'q000=',1000.*qtotd_0_evap            &
            !   ,'qtot=',1000.*qtotd_cld(k),'evap=',1000.*evapd_cld(k)                      &
            !   ,'qsat=',1000.*qsatd_cld(k),'qvap=',1000.*qvapd_cld(k)                      &
            !   ,'qliq=',1000.*qliqd_cld(k),'qice=',1000.*qiced_cld(k)                      &
            !   ,'temp=',td_cld(k)-t00,'funnow=',funnow
            !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
            !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

            !------------------------------------------------------------------------------!
            ! e4. Testing for convergence, depending on the method.                        !
            !------------------------------------------------------------------------------!
            if (funnow == 0.) then
               converged = .true.
            elseif (bisection) then 
               converged = abs(qtotd_cld(k)-qtotda) < toler*qtotd_cld(k)
            else
               converged = abs(qtotd_cld(k)-qtotdc) < toler*qtotd_cld(k)
            end if
            !----- Found a good set, leaving... -------------------------------------------!
            if (converged) exit itloop
            !------------------------------------------------------------------------------!


            !------------------------------------------------------------------------------!
            ! e5. Checking which side from the A Z pair I can update, based on funnow.     !
            !------------------------------------------------------------------------------!
            if (funnow*funa < 0.) then
               funz   = funnow
               qtotdz = qtotd_cld(k)
            else
               funa   = funnow
               qtotda = qtotd_cld(k)
            end if

            !------------------------------------------------------------------------------!
            ! e6. Updating the Previous-Current pair for the next secant attempt.          !
            !------------------------------------------------------------------------------!
            qtotdp = qtotdc
            funp   = func
            qtotdc = qtotd_cld(k)
            func   = funnow

         end do itloop
      end if iterif
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write(unit=28,fmt='(a)') '-----------------------------------------------------------'
      !write(unit=28,fmt='(a)') ' '
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

      if (.not. converged) then
         write (unit=*,fmt='(a)') '-------------------------------------------------------'
         write (unit=*,fmt='(a)') ' Fall-out water finding didn''t converge!!!'
         write (unit=*,fmt='(a,1x,i5,1x,a)') ' I gave up, after',maxfpo,'iterations...'
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Input values.'
         write (unit=*,fmt='(a,1x,i5)') 'Level k= ',k
         write (unit=*,fmt='(a,1x,f12.4)' ) 'qtotd_0_evap      [g/kg] =',1000.*qtotd_0_evap
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Dnd. Thetae_iv    [   K] =',theivd_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'p_cup             [ hPa] =',0.01*p_cup(k)
         write (unit=*,fmt='(a)') ' '
         write (unit=*,fmt='(a)') ' Last iteration outcome (downdraft values).'
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Total  mix. ratio [g/kg] =',1000.*qtotd_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Vapour mix. ratio [g/kg] =',1000.*qvapd_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Liquid mix. ratio [g/kg] =',1000.*qliqd_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Ice    mix. ratio [g/kg] =',1000.*qiced_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Evap.  mix. ratio [g/kg] =',1000.*evapd_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Theta_il          [   K] =',thild_cld(k)
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Temperature       [  °C] =',td_cld(k)-t00
         write (unit=*,fmt='(a,1x,f12.4)' ) 'Function          [g/kg] =',1000.*funnow
         write (unit=*,fmt='(a,1x,es12.4)') 'Error (secant)    [ ---] ='                   &
                                            ,abs(qtotd_cld(k)-qtotdc)/qtotd_cld(k) 
         write (unit=*,fmt='(a,1x,es12.4)') 'Error (bisection) [ ---] ='                   &
                                            ,abs(qtotd_cld(k)-qtotda)/qtotd_cld(k) 
         write (unit=*,fmt='(a)') '-------------------------------------------------------' 

         call abort_run('Couldn''t find the evaporation, the zeroin method diverged...'    &
                       ,'grell_most_thermo_downdraft','grell_cupar_downdraft.f90')
      end if
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!

      !------------------------------------------------------------------------------------!
      !    Evaporation is supplied by the rain that fell out from the updrafts and some    !
      ! small amount of liquid/ice that may exist in the downdraft due to entrainment of   !
      ! saturated air.                                                                     !
      !------------------------------------------------------------------------------------!
      pwd_cld(k) = etad_cld(k) * evapd_cld(k)
      pwev       = pwev        + pwd_cld(k)

      !------ Finding density, assuming pd_cld(k) ~= p_cup(k)... --------------------------!
      rhod_cld(k) = idealdens(p_cup(k),td_cld(k),qvapd_cld(k),qtotd_cld(k))

      !------ Finding buoyancy ------------------------------------------------------------!
      dbyd(k) = buoyancy_acc(rho_cup(k),rhod_cld(k))

      bytot   = bytot + dbyd(k) * dzd_cld(k)

   end do

   !---------------------------------------------------------------------------------------!
   ! 3. Now we should check whether the downdraft buoyancy makes sense. Downdrafts, as we  !
   !    can imagine from their names, must go down... This will be possible only if there  !
   !    is negative acceleration, provided by buoyancy. If that's not the case this cloud  !
   !    makes no sense so its existence will be denied...                                  !
   !---------------------------------------------------------------------------------------!
   if (bytot > 0.) ierr = 8

   return
end subroutine grell_most_thermo_downdraft
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!  This subroutine computes the cloud work function associated with downdrafts             !
!------------------------------------------------------------------------------------------!
subroutine grell_cldwork_downdraft(mkx,mgmzp,klod,dbyd,dzd_cld,etad_cld,aad)
   use rconstants, only : gocp
   implicit none

   integer               , intent(in)  :: mkx, mgmzp  ! Grid dimesnsions
   integer               , intent(in)  :: klod        ! Downdraft origin level

   real, dimension(mgmzp), intent(in)  :: dbyd        ! Buoyancy term               [ J/kg]
   real, dimension(mgmzp), intent(in)  :: dzd_cld     ! Delta-z for downdrafts      [    m]
   real, dimension(mgmzp), intent(in)  :: etad_cld    ! Normalized mass flux        [  ---]
   real                  , intent(out) :: aad         ! Downdraft work function     [ J/kg]

   integer                             :: k           ! Counter

   !----- Initialize cloud work to zero. --------------------------------------------------!
   aad = 0.
   
   !---------------------------------------------------------------------------------------!
   !    The cloud work is a measure of efficiency of kinetic energy generation inside the  !
   ! cloud, and therefore it is directly proportional to the mass flux and buoyancy. The   !
   ! final value should represent the cloud function for the entire downdraft layer thus   !
   ! the integral between the surface and the level in which downdrafts originate cloud    !
   ! base and cloud top. Note that dbyd is negative so this will give a positive function. !
   !---------------------------------------------------------------------------------------!
   do k=klod-1,1,-1
      aad = aad - etad_cld(k) * dbyd(k) *dzd_cld(k)
   end do
   aad = max(0.,aad)

   return
end subroutine grell_cldwork_downdraft
!==========================================================================================!
!==========================================================================================!
