!==========================================================================================!
! grell_cupar_aux.f90                                                                      !
!                                                                                          !
!    This file contains subroutines that will handle the copying from BRAMS to Grell and   !
! vice-versa. This means that all initialization and output routines are here.             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This tiny subroutire initialises Grell's grid. This grid goes skips the bottom bound- !
! ary condition and in case of adaptive coordinate, the levels below surface               !
!------------------------------------------------------------------------------------------!
subroutine initial_grid_grell(m1,deltax,deltay,zt,zm,flpw,rtgt,confrq,pblidx)
   use mem_scratch_grell, only: &
           mkx                  & ! intent(out) - Number of Grell levels.
          ,kgoff                & ! intent(out) - BRAMS offset related to Grell
          ,kpbl                 & ! intent(out) - Level of PBL top
          ,ktpse                & ! intent(out) - Maximum cloud top allowed
          ,lpw                  & ! intent(out) - Lowest thermodynamic point
          ,tscal_kf             & ! intent(out) - Kain-Fritsch(1990) time scale
          ,z                    & ! intent(out) - Height
          ,z_cup                & ! intent(out) - Height at cloud levels
          ,dzu_cld              & ! intent(out) - Delta z for updraft calculations
          ,dzd_cld              ! ! intent(out) - Delta z for downdraft calculations
   use grell_coms       , only: &
           zmaxtpse             ! ! intent(in)  - Maximum height allowed for cloud top.
   implicit none
   !------ I/O variables ------------------------------------------------------------------!
   integer               , intent(in)  :: m1        ! Number of vertical levels
   real                  , intent(in)  :: deltax    ! Zonal resolution
   real                  , intent(in)  :: deltay    ! Meridional resolution
   real   , dimension(m1), intent(in)  :: zt        ! Height at thermodynamic levels
   real   , dimension(m1), intent(in)  :: zm        ! Height at momentum levels
   real                  , intent(in)  :: rtgt      ! Correction for terrain-following
   real                  , intent(in)  :: flpw      ! Level of lowest point above ground
   real                  , intent(in)  :: confrq    ! How often is this cloud called?
   integer               , intent(in)  :: pblidx    ! Level of PBL top (Nakanishi/Niino)

   !------ Local variables ----------------------------------------------------------------!
   integer                             :: k

  !------ Initializing grid variables -----------------------------------------------------!
   lpw   = nint(flpw) ! This is always 2 for terrain-following, but varies for adaptive.
   kgoff = lpw-1      ! K-offset between BRAMS and Grell's grids.
   mkx   = m1-kgoff   ! Number of Grell levels.
   
   !------ Grid-dependent Kain-Fritsch time scale, based on Betts suggestion --------------!
   !tscal_kf = 0.02 * sqrt(deltax*deltay)
   !------ Using the frequency of call as the time scale ----------------------------------!
   !tscal_kf = confrq
   !------ Using the original: 3000.s -----------------------------------------------------!
   tscal_kf  = 3000.

   !------ PBL top. If running Nakanishi/Niino only, otherwise it needs to be found. ------!
   if (pblidx /= 0) then
      kpbl = pblidx - kgoff
   else
      kpbl = 0
   end if
   !------ Height, used to compute wind average (unecessary for pressure means) -----------!
   do k=1,mkx
      z(k)   = (zt(k+kgoff)-zm(kgoff))*rtgt
   end do
   !----- Zcup 1 is always 0, height uses ground as reference -----------------------------!
   z_cup(1) = 0.
   do k=2,mkx
      z_cup(k)= 0.5 * (z(k-1) + z(k))
   end do

   !---------------------------------------------------------------------------------------!
   !    Computing dz associated with downdrafts and updrafts. Downdraft properties are     !
   ! usually computed from top to bottom, whereas updraft properties go the other way.     !
   ! Thus delta-z must reflect this direction. For boundaries, I will use the same value   !
   ! as the valid neighbour, but they should never be used along the routine. The values   !
   ! won't be computed for downdrafts in case this cloud don't have them.                  !
   !---------------------------------------------------------------------------------------!
   !----- Downdrafts ----------------------------------------------------------------------!
   do k=1,mkx-1
      dzd_cld(k) = z_cup(k+1) - z_cup(k)
   end do
   dzd_cld(mkx)  = z_cup(mkx)-z_cup(mkx-1)
   !----- Updrafts ------------------------------------------------------------------------!
   dzu_cld(1) = z_cup(2)-z_cup(1)
   do k=2,mkx
      dzu_cld(k) = z_cup(k)-z_cup(k-1)
   end do
   
   !----- Finding the top height that we allow convection to develop. ---------------------!
   pauseloop: do ktpse=2,mkx-1
      if (z_cup(ktpse) > zmaxtpse) exit pauseloop
   end do pauseloop
   ktpse = ktpse - 1

   return
end subroutine initial_grid_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine simply copies the tendencies to scratch arrays. It is done separated- !
! ly because it is the only place that we must give i and j information.                   !
!------------------------------------------------------------------------------------------!
subroutine initial_tend_grell(m1,tht,tket,rtt,co2t)
   use mem_scratch_grell, only : &
           mkx          & ! intent(in)  - Number of Grell levels.
          ,kgoff        & ! intent(in)  - BRAMS offset related to Grell
          ,lpw          & ! intent(in)  - Lowest thermodynamic point
          ,dthildt      & ! intent(out) - Theta_il tendency
          ,dtkedt       & ! intent(out) - TKE tendency
          ,dqtotdt      & ! intent(out) - Total H2O mixing ratio tendency
          ,dco2dt       ! ! intent(out) - Total CO2 mixing ratio tendency
   use rconstants, only: day_sec
   implicit none
   !------ I/O variables ------------------------------------------------------------------!
   integer, intent(in)                :: m1   ! Grid dimensions
   real   , intent(in), dimension(m1) :: tht  ! Potential temperature tend.
   real   , intent(in), dimension(m1) :: tket ! Turbulent Kinetic Energy tend.
   real   , intent(in), dimension(m1) :: co2t ! Total CO2 mixing ratio tend.
   real   , intent(in), dimension(m1) :: rtt  ! Total H2O mixing ratio tend.
   !------ Local variables ----------------------------------------------------------------!
   integer                            :: k    ! Counter
   integer                            :: kr   ! Counter
   !---------------------------------------------------------------------------------------!

   !----- Here we simply copy the variables, removing the offset --------------------------!
   
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   !write (unit=23,fmt='(a)') '-------------------------------------------------------------'
   !write (unit=23,fmt='(3(a,1x,i5,1x))') 'i=',i,'j=',j,'mkx=',mkx
   !write (unit=23,fmt='(2(1x,a5),3(1x,a13))') adjustr('k'),adjustr('kgoff')                &
   !          ,adjustr('dthildt'),adjustr('dqtotdt'),adjustr('dtkedt'),adjustr('dco2dt')
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

   do k=1,mkx

      kr=k+kgoff
      dthildt(k)      = tht(kr)
      dqtotdt(k)      = rtt(kr)
      dtkedt(k)       = tket(kr)
      dco2dt(k)       = co2t(kr)

      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
      !write (unit=23,fmt='(2(1x,i5),5(1x,es13.6))') k,kgoff                                &
      !      ,day_sec*dthildt(k),day_sec*dqtotdt(k),day_sec*dtkedt(k)                       &
      !      ,day_sec*dthildt_shal(k),day_sec*dqtotdt_shal(k),day_sec*dco2dt_shal(k)
      !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
      !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!

   end do
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
   !write (unit=23,fmt='(a)') '-------------------------------------------------------------'
   !write (unit=23,fmt='(a)') ' '
   !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><!
   !><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!

   return
end subroutine initial_tend_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine initialises Grell's thermodynamic and turbulence fields. If it        !
! is a deep convection call, then it will include the effect of shallow convection, and    !
! make the variables for the case in which no convection happens consistent with this.     !
!------------------------------------------------------------------------------------------!
subroutine initial_thermo_grell(m1,mgmzp,dtime,thp,theta,rtp,co2p,pi0,pp,pc,wp,dn0,tkep    &
                               ,rliq,rice,wstd)

   use mem_scratch_grell, only : &
          dco2dt       & ! intent(in)  - Total CO2 mixing ratio tendency          [  ppm/s]
         ,dthildt      & ! intent(in)  - Temporary theta_il tendency              [    K/s]
         ,dqtotdt      & ! intent(in)  - Total H2O mixing ratio tendency          [kg/kg/s]
         ,dtkedt       & ! intent(in)  - Temporary TKE tendency                   [ J/kg/s]
         ,mkx          & ! intent(in)  - Number of Grell levels.                  [    ---]
         ,kgoff        & ! intent(in)  - BRAMS offset related to Grell            [    ---]
         ,lpw          & ! intent(in)  - Lowest thermodynamic point               [    ---]
         ,z            & ! intent(in)  - Grell's heights                          [      m]
         ,co20         & ! intent(out) - CO2 mixing ratio                         [    ppm]
         ,co2          & ! intent(out) - Forced CO2 mixing ratio                  [    ppm]
         ,co2sur       & ! intent(out) - Sfc. CO2 mixing ratio                    [    ppm]
         ,exner0       & ! intent(out) - Exner function                           [ J/kg/K]
         ,exner        & ! intent(out) - Forced Exner function                    [ J/kg/K]
         ,exnersur     & ! intent(out) - Surface Exner function                   [ J/kg/K]
         ,mconv        & ! intent(out) - Moisture convergence                     [kg/m²/s]
         ,omeg         & ! intent(out) - Lagrangian pressure tendency             [   Pa/s]
         ,p0           & ! intent(out) - Current pressure                         [     Pa]
         ,p            & ! intent(out) - Forced Pressure                          [     Pa]
         ,psur         & ! intent(out) - Surface pressure                         [     Pa]
         ,qtot0        & ! intent(out) - Current Total mixing ratio               [  kg/kg]
         ,qtot         & ! intent(out) - Forced Total mixing ratio                [  kg/kg]
         ,qtotsur      & ! intent(out) - Sfc. Total mixing ratio                  [  kg/kg]
         ,qvap0        & ! intent(out) - Current Vapour mixing ratio              [  kg/kg]
         ,qvap         & ! intent(out) - Water vapour mixing ratio w/ forcing     [  kg/kg]
         ,qvapsur      & ! intent(out) - Sfc. Water vapour mixing ratio           [  kg/kg]
         ,qliq0        & ! intent(out) - Liquid water mixiing ratio               [  kg/kg]
         ,qliq         & ! intent(out) - Liquid water mixing ratio w/ forcing     [  kg/kg]
         ,qliqsur      & ! intent(out) - Sfc. liquid water mixing ratio           [  kg/kg]
         ,qice0        & ! intent(out) - Ice mixing ratio                         [  kg/kg]
         ,qice         & ! intent(out) - Ice mixing ratio with forcing            [  kg/kg]
         ,qicesur      & ! intent(out) - Sfc. ice mixing ratio                    [  kg/kg]
         ,rho0         & ! intent(out) - Density                                  [  kg/m³]
         ,rho          & ! intent(out) - Density with forcing                     [  kg/m³]
         ,rhosur       & ! intent(out) - Sfc. density                             [  kg/m³]
         ,t0           & ! intent(out) - Current Temperature                      [      K]
         ,t            & ! intent(out) - Forced Temperature                       [      K]
         ,tsur         & ! intent(out) - Sfc. Temperature                         [      K]
         ,theiv0       & ! intent(out) - Current ice-vapour equiv. pot. temp.     [      K]
         ,theiv        & ! intent(out) - Forced ice-vapour equiv. pot. temp.      [      K]
         ,theivsur     & ! intent(out) - Sfc. ice-vapour equiv. pot. temp.        [      K]
         ,thil0        & ! intent(out) - Current ice-liquid potential temperature [      K]
         ,thil         & ! intent(out) - Forced ice-liquid potential temperature  [      K]
         ,thilsur      & ! intent(out) - Sfc. ice-liquid potential temperature    [      K]
         ,tke0         & ! intent(out) - Turbulent kinetic energy                 [   J/kg]
         ,tke          & ! intent(out) - Forced Turbulent kinetic energy          [   J/kg]
         ,sigw         & ! intent(out) - Vertical velocity standard deviation     [    m/s]
         ,wwind        ! ! intent(out) - Mean vertical velocity                   [    m/s]
   use rconstants, only : grav    & ! intent(in)
                        , rdry    & ! intent(in)
                        , epi     & ! intent(in)
                        , toodry  & ! intent(in)
                        , sigwmin & ! intent(in)
                        , tkmin   ! ! intent(in)

   !------ External functions -------------------------------------------------------------!
   use therm_lib, only : &
            rslif        & ! Function that finds the saturation mixing ratio
          , thetaeiv     & ! Function that finds Thetae_iv  
          , thil2temp    & ! Function that gives temperature from theta_il, rliq and rice.
          , thil2tqall   & ! Function that finds temperature and condensed phase from thil.
          , idealdens    & ! Function that gives the density for ideal gasses
          , exner2press  & ! Function that converts Exner function into pressure
          , extheta2temp ! ! Function that finds pot. temp. from Exner and temperature
   implicit none
   !------ I/O variables ------------------------------------------------------------------!
   integer, intent(in)                  :: m1     ! Grid dimension               [     ---]
   integer, intent(in)                  :: mgmzp  ! Dim. of the scratch arrays   [     ---]
   real   , intent(in)                  :: dtime  ! Time step (convective scale) [       s]
   real   , intent(in)  , dimension(m1) :: thp    ! Ice-liquid potential temp.   [       K]
   real   , intent(in)  , dimension(m1) :: theta  ! Potential temperature        [       K]
   real   , intent(in)  , dimension(m1) :: rtp    ! Total H2O mixing ratio       [   kg/kg]
   real   , intent(in)  , dimension(m1) :: co2p   ! Total CO2 mixing ratio       [umol/mol]
   real   , intent(in)  , dimension(m1) :: pi0    ! Reference Exner function     [  J/kg/K]
   real   , intent(in)  , dimension(m1) :: pp     ! Current perturbation on pi   [  J/kg/K]
   real   , intent(in)  , dimension(m1) :: pc     ! Future perturbation on pi    [  J/kg/K]
   real   , intent(in)  , dimension(m1) :: dn0    ! Reference density            [   kg/m³]
   real   , intent(in)  , dimension(m1) :: wp     ! Vertical velocity            [     m/s]
   real   , intent(in)  , dimension(m1) :: tkep   ! Turbulent kinetic energy     [    J/kg]
   real   , intent(in)  , dimension(m1) :: rliq   ! Liquid water mixing ratio    [   kg/kg]
   real   , intent(in)  , dimension(m1) :: rice   ! Ice mixing ratio             [   kg/kg]
   real   , intent(in)  , dimension(m1) :: wstd   ! Standard deviation of wp     [     m/s]
   !------ Local variables ----------------------------------------------------------------!
   integer                              :: k      ! Counters                     [     ---]
   integer                              :: kr     ! Counters                     [     ---]
   real                                 :: dq     ! Diff. on vapour mixing ratio [   kg/kg]
   real                                 :: qsat   ! Sat. mixing ratio (scratch)  [   kg/kg]
   !------ Surface variables, copied to a dummy array because of interface check. ---------!
   real                  , dimension(1) :: z1     ! Height                       [       m]
   real                  , dimension(1) :: exner1 ! Exner function               [  J/kg/K]
   real                  , dimension(1) :: p1     ! Pressure                     [      Pa]
   real                  , dimension(1) :: thil1  ! Ice-liquid potential temp.   [       K]
   real                  , dimension(1) :: qtot1  ! Total water mixing ratio     [   kg/kg]
   real                  , dimension(1) :: qliq1  ! Liquid water mixing ratio    [   kg/kg]
   real                  , dimension(1) :: qice1  ! Ice mixing ratio             [   kg/kg]
   real                  , dimension(1) :: qvap1  ! Vapour mixing ratio          [   kg/kg]
   real                  , dimension(1) :: t1     ! Temperature                  [       K]
   real                  , dimension(1) :: theiv1 ! Equiv. ice vapour pot. temp. [       K]
   real                  , dimension(1) :: co21   ! CO2 mixing ratio             [umol/mol] 
   real                  , dimension(1) :: rho1   ! Density                      [   kg/m³]
   !---------------------------------------------------------------------------------------!

   do k=1,mkx
      kr=k+kgoff
      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      !    Find the current state variables, including the effect of shallower cumulus     !
      ! if that is the case.                                                               !
      !------------------------------------------------------------------------------------!
      !------ 1. Exner function -----------------------------------------------------------!
      exner0(k) = pi0(kr)   + pp(kr)
      !------ 2. Pressure. ----------------------------------------------------------------!
      p0(k)     = exner2press(exner0(k))
      !------------------------------------------------------------------------------------!
      ! 3. Temperature and water.                                                          !
      !------------------------------------------------------------------------------------!
      thil0(k)  = thp(kr)
      qtot0(k)  = max(toodry,rtp(kr))
      qliq0(k)  = max(0.,rliq(kr))
      qice0(k)  = max(0.,rice(kr))
      qvap0(k)  = max(toodry,qtot0(k)-qice0(k)-qliq0(k))
      t0(k)     = extheta2temp(exner0(k),theta(kr))
      call grell_sanity_thil2tqall(k,z(k),thil0(k),exner0(k),p0(k),qtot0(k),'thermo_zero')

      !------ 4. Finding the ice-vapour equivalent potential temperature ------------------!
      theiv0(k) = thetaeiv(thil0(k),p0(k),t0(k),qvap0(k),qtot0(k))

      !------ 5. CO2 mixing ratio. --------------------------------------------------------!
      co20(k)   = co2p(kr)
      !------ 6. Turbulent kinetic energy [m²/s²] -----------------------------------------!
      tke0(k)   = tkep(kr)
      !------ 7. Vertical velocity in terms of pressure, or Lagrangian dp/dt [ Pa/s] ------!
      omeg(k)   = -grav*dn0(kr)*.5*( wp(kr)+wp(kr-1) )
      !------ 8. Vertical velocity [m/s], this is staggered, averaging... -----------------!
      wwind(k)  = 0.5 * (wp(kr)+wp(kr-1))
      !------ 9. Standard-deviation of vertical velocity ----------------------------------!
      sigw(k)   = max(wstd(kr),sigwmin)
      !------ 10. Air density -------------------------------------------------------------!
      rho(k)    = idealdens(p0(k),t0(k),qvap0(k),qtot0(k))
      !------------------------------------------------------------------------------------!
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!



      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      !     Find what the state variables will be in the next time, assuming no convec-    !
      ! tion at this point (we will call these forced variables).  Most variables will be  !
      ! updated using the tendency, except for the Exner function and diagnostic vari-     !
      !ables.                                                                              !
      !------------------------------------------------------------------------------------!
      !------ 1. Exner function pc is the future Exner perturbation. ----------------------!
      exner(k) = pi0(kr)   + pc(kr)
      !------ 2. Pressure -----------------------------------------------------------------!
      p(k)     = exner2press(exner(k))
      !------ 3. Ice-liquid potential temperature, with the tendency ----------------------!
      thil(k)  = thp(kr) + dthildt(k)*dtime
      !------ 4. Total mixing ratio, with the tendency ------------------------------------!
      qtot(k)  = max(toodry,rtp(kr)   + dqtotdt(k) * dtime)
      !------ 5. Find the equilibrium state. Temperature 1st guess is simply t0. ----------!
      t(k)     = t0(k)
      call grell_sanity_thil2tqall(k,z(k),thil(k),exner(k),p(k),qtot(k),'thermo_extrap')
      call thil2tqall(thil(k),exner(k),p(k),qtot(k),qliq(k),qice(k),t(k),qvap(k),qsat)
      !------ 6. Finding the ice-vapour equivalent potential temperature ------------------!
      theiv(k) = thetaeiv(thil(k),p(k),t(k),qvap(k),qtot(k))
      !------ 7. CO2 mixing ratio ---------------------------------------------------------!
      co2(k)   = co2p(kr) + dco2dt(k) * dtime
      !------ 8. Turbulent kinetic energy -------------------------------------------------!
      tke(k)   = max(tkmin,tkep(kr) + dtkedt(k) * dtime)
      !------ 9. Air density --------------------------------------------------------------!
      rho(k)   = idealdens(p(k),t(k),qvap(k),qtot(k))
      !------------------------------------------------------------------------------------!
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!

   end do


   ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
   !---------------------------------------------------------------------------------------!
   !     Compute the surface variables. The only one that will be truly different is the   !
   ! Exner function (and consequently, pressure). The other values will be based on the    !
   ! level above. This is going to be just a boundary condition, so they will directly     !
   ! affect the parametrisation.                                                           !
   !---------------------------------------------------------------------------------------!
   !----- 0. Height, always 0. ------------------------------------------------------------!
   z1(1)       = 0.
   !----- 1. Exner function ---------------------------------------------------------------!
   exnersur    = sqrt((pp(lpw-1)+pi0(lpw-1))*(pp(lpw)+pi0(lpw)))
   exner1(1)   = exnersur
   !----- 2. Pressure ---------------------------------------------------------------------!
   psur        = exner2press(exnersur)
   p1(1)       = psur
   !----- 3. Ice liquid potential temperature ---------------------------------------------!
   thilsur     = thp(lpw)
   thil1(1)    = thilsur
   !----- 4. Total mixing ratio -----------------------------------------------------------!
   qtotsur     = max(toodry,rtp(lpw))
   qtot1(1)    = qtotsur
   !----- 5. Liquid water mixing ratio ----------------------------------------------------!
   qliqsur     = max(0.,rliq(lpw))
   qliq1(1)    = qliqsur
   !----- 6. Ice mixing ratio -------------------------------------------------------------!
   qicesur     = max(0.,rice(lpw))
   qice1(1)    = qicesur
   !----- 7. Vapour mixing ratio ----------------------------------------------------------!
   qvapsur     = max(toodry,qtotsur-qliqsur-qicesur)
   qvap1(1)    = qvapsur
   !----- 7. Temperature ------------------------------------------------------------------!
   tsur        = extheta2temp(exnersur,theta(lpw))
   t1(1)       = tsur
   !----- 8. Ice-vapour equivalent potential temperature ----------------------------------!
   theivsur    = thetaeiv(thilsur,psur,tsur,qvapsur,qtotsur)
   theiv1(1)   = theivsur
   !----- 9. CO2 mixing ratio -------------------------------------------------------------!
   co2sur      = co2p(lpw)
   co21(1)     = co2sur
   !----- 10. Air density -----------------------------------------------------------------!
   rhosur      = idealdens(psur,tsur,qvapsur,qtotsur)
   rho1(1)     = rhosur
   call grell_sanity_thil2tqall(1,z1(1),thil1(1),exner1(1),p1(1),qtot1(1),'thermo_surface')
   !---------------------------------------------------------------------------------------!
   !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!


   !---------------------------------------------------------------------------------------!
   !     Check the profiles we will use.                                                   !
   !---------------------------------------------------------------------------------------!
   call grell_sanity_check(mkx,mgmzp,z,p0,exner0,theiv0,thil0,t0,qtot0,qvap0,qliq0         &
                          ,qice0,co20,rho0,'thermo_zero')
   call grell_sanity_check(mkx,mgmzp,z,p,exner,theiv,thil,t,qtot,qvap,qliq,qice,co2        &
                          ,rho,'thermo_extrap')
   call grell_sanity_check(1,1,z1,p1,exner1,theiv1,thil1,t1,qtot1,qvap1,qliq1,qice1,co21   &
                          ,rho1,'thermo_surface')
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find the integrated moisture convergence. This is done outside the loop so we     !
   ! can use vapour mixing ratio at level (k-1) and (k+1).                                 !
   !---------------------------------------------------------------------------------------!
   mconv=0.
   do k=2,mkx-1
      dq          = .5*(qvap(k+1)-qvap(k-1))  ! Delta-q between layers.
      !------ This is moisture convergence, meaning (-1) × moisture divergence. -----------!
      mconv       = mconv+omeg(k)*dq/grav
   end do
   !----- Moisture convergence closure is only interested when convergence is positive. ---!
   mconv     = max(0.,mconv)
   !---------------------------------------------------------------------------------------!

   return
end subroutine initial_thermo_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine intialises the wind information and the previous downdraft mass flux,  !
! which may affect current convection depending on the dynamic control chosen.             !
!------------------------------------------------------------------------------------------!
subroutine initial_winds_grell(prec_cld,m1,m2,m3,i,j,jdim,last_dnmf,ua,va,prev_dnmf)

   use mem_scratch_grell, only: &
            mkx                 & ! intent(in)  - Number of Grell levels.
           ,kgoff               & ! intent(in)  - BRAMS offset related to Grell
           ,lpw                 & ! intent(in)  - Lowest thermodynamic point
           ,p                   & ! intent(in)  - Pressure at Grell's levels
           ,psur                & ! intent(in)  - Pressure at surface
           ,uwind               & ! intent(out) - Zonal wind at thermodynamic point
           ,vwind               & ! intent(out) - Meridional wind at thermodynamic point
           ,z                   ! ! intent(in)  - Height at Grell's levels

   implicit none
   !------ Arguments. ---------------------------------------------------------------------!
   integer                      , intent(in)    :: m1        ! Number of z points
   integer                      , intent(in)    :: m2        ! Number of x points
   integer                      , intent(in)    :: m3        ! Number of y points
   integer                      , intent(in)    :: i,j       ! Current x, y & cld position
   integer                      , intent(in)    :: jdim      ! Dimension in y
   logical                      , intent(in)    :: prec_cld  ! Precipitating cloud (T/F)
   real   , dimension(m2,m3)    , intent(in)    :: last_dnmf ! Last time downdraft
   real   , dimension(m1,m2,m3) , intent(in)    :: ua        ! Zonal wind
   real   , dimension(m1,m2,m3) , intent(in)    :: va        ! Meridional wind
   real   , dimension(1)        , intent(inout) :: prev_dnmf ! Previous downdraft
   !------ Local variables ----------------------------------------------------------------!
   integer            :: k         ! Counter for current Grell level
   integer            :: kr        ! Counter for corresponding BRAMS level
   !---------------------------------------------------------------------------------------!

   
   !------ Initializing scalars -----------------------------------------------------------!
   prev_dnmf(1) = last_dnmf(i,j)
   !---------------------------------------------------------------------------------------!
   !    Transferring the values from BRAMS to Grell's levels, remembering that Grell's     !
   ! grid goes from 1 to m1-1.                                                             !
   !---------------------------------------------------------------------------------------!
   do k = 1, mkx
      !----- kr is the BRAMS level. I included lpw so it works for adaptive coordinate too.!
      kr=k+kgoff

      !------ Wind needs to be interpolated to the thermodynamic point --------------------!
      uwind(k) = .5*( ua(kr,i,j) + ua(kr,i-1,j)    )    ! Zonal wind                [  m/s]
      vwind(k) = .5*( va(kr,i,j) + va(kr,i,j-jdim) )    ! Meridional wind           [  m/s]
   end do
   !---------------------------------------------------------------------------------------!

   return
end subroutine initial_winds_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine computes the relative downdraft and updraft areas. This is used by    !
! some Lagrangian models and it will be also used to feedback the liquid water content.    !
!                                                                                          !
!  The references for this code are:                                                       !
!                                                                                          !
! Fritsch, J.M; Chappell, C. F., 1980: Numerical prediction of convectively driven         !
!      mesoscale pressure systems. Part I: Convective parameterization. J. Atmos. Sci.,    !
!      vol. 37(8), 1722-1733. (fc or FC80).                                                !
!                                                                                          !
! Zhang, D.-L.; Fritsch, J.M., 1986: Numerical simulation of the meso-beta scale structure !
!      and evolution of the 1977 Johnstown flood. Part I: Model description and            !
!      verification. J. Atm. Sci., vol. 43(18). 1913-1943. (zf or ZF86).                   !
!                                                                                          !
!    One important difference from FC80/ZF86 is that here we considered not only the       !
! entrainment, but also the detrainment to estimate the draft velocities, and the entrain- !
! ment of environment kinetic energy (that should have a really small contribution...)     !
!------------------------------------------------------------------------------------------!
subroutine grell_draft_area(m1,mgmzp,kgoff,comp_down,klod,klou,klfc,klnb,ktop,dzu_cld      &
                           ,wwind,tke,sigw,wbuoymin,etad_cld,mentrd_rate,cdd,dbyd,rhod_cld &
                           ,dnmf,etau_cld,mentru_rate,cdu,dbyu,rhou_cld,upmf,areadn,areaup &
                           ,wdndraft,wupdraft)
   use rconstants, only : sigwmin
   implicit none
   
   !----- Input variables, flags ----------------------------------------------------------!
   logical               , intent(in)  :: comp_down   ! Flag for downdraft computation
   !----- Input variables, grid boundaries ------------------------------------------------!
   integer               , intent(in)  :: m1          ! Number of BRAMS levels
   integer               , intent(in)  :: mgmzp       ! Number of scratch levels
   integer               , intent(in)  :: kgoff       ! BRAMS offset
   !----- Input variables, cloud boundaries -----------------------------------------------!
   integer               , intent(in)  :: klod        ! Downdraft origin
   integer               , intent(in)  :: klou        ! Updraft origin
   integer               , intent(in)  :: klfc        ! Level of free convection
   integer               , intent(in)  :: klnb        ! Level of neutral buoyancy
   integer               , intent(in)  :: ktop        ! Level of neutral buoyancy
   !----- Input variables, vertical velocity-related variables ----------------------------!
   real, dimension(mgmzp), intent(in)  :: dzu_cld     ! Bottom-up layer thickness [      m]
   real, dimension(mgmzp), intent(in)  :: wwind       ! Vertical velocity         [    m/s]
   real, dimension(mgmzp), intent(in)  :: tke         ! Turbulent Kinetic Energy  [   J/kg]
   real, dimension(mgmzp), intent(in)  :: sigw        ! Std. deviation of W       [    m/s]
   real                  , intent(in)  :: wbuoymin    ! Minimum buoyant velocity  [    m/s]
   !----- Input variables, downdraft properties -------------------------------------------!
   real, dimension(mgmzp), intent(in)  :: etad_cld    ! Normalized mass flux      [    ---]
   real, dimension(mgmzp), intent(in)  :: mentrd_rate ! Normalized entrainment   .[    1/m]
   real, dimension(mgmzp), intent(in)  :: cdd         ! Normalized detrainment   .[    1/m]
   real, dimension(mgmzp), intent(in)  :: dbyd        ! Buoyancy acceleration    .[   m/s²]
   real, dimension(mgmzp), intent(in)  :: rhod_cld    ! Density                   [  kg/m³]
   real                  , intent(in)  :: dnmf        ! Reference mass flux       [kg/m²/s]
   !----- Input variables, updraf properties ----------------------------------------------!
   real, dimension(mgmzp), intent(in)  :: etau_cld    ! Normalized mass flux      [    ---]
   real, dimension(mgmzp), intent(in)  :: mentru_rate ! Normalized entrainment   .[    1/m]
   real, dimension(mgmzp), intent(in)  :: cdu         ! Normalized detrainment   .[    1/m]
   real, dimension(mgmzp), intent(in)  :: dbyu        ! Buoyancy acceleration    .[   m/s²]
   real, dimension(mgmzp), intent(in)  :: rhou_cld    ! Density                   [  kg/m³]
   real                  , intent(in)  :: upmf        ! Reference mass flux       [kg/m²/s]
   !----- Output variables ----------------------------------------------------------------!
   real                  , intent(out) :: areadn      ! Downdraft relative area   [    ---]
   real                  , intent(out) :: areaup      ! Updraft   relative area   [    ---]
   real                  , intent(out) :: wdndraft    ! Downdraft at LOD          [    m/s]
   real                  , intent(out) :: wupdraft    ! Updraft   at LOU          [    m/s]
   !----- Local variables -----------------------------------------------------------------!
   integer                             :: k         ! Cloud level counter
   integer                             :: kr        ! BRAMS level counter
   real                                :: lske2     ! Scratch large-scale KE * 2  [   J/kg]
   real                                :: wnorm     ! Normalised buoyant velocity [    m/s]
   real                                :: cdfval    ! Scratch with CDF value      [    ---]
   !----- Minimum draft to prevent square root of negative values -------------------------!
   real, parameter                     :: min_downdraft2      = 0.01
   real, parameter                     :: min_wupdraft        = 0.60
   real, parameter                     :: max_wupdraft        = 7.00
   real, parameter                     :: max_areaup          = 0.50
   !----- External functions --------------------------------------------------------------!
   real, external                      :: cdf
   real, external                      :: expected
   !---------------------------------------------------------------------------------------!

   
   !---------------------------------------------------------------------------------------!
   !    Setting the area to zero. This will be replaced by the actual area at the draft    !
   ! layer.                                                                                !
   !---------------------------------------------------------------------------------------!
   areadn   = 0.
   areaup   = 0.
   wdndraft = 0.
   wupdraft = 0.
   
   !---------------------------------------------------------------------------------------!
   !    Downdraft. As in Fritsch and Chappell (1980), I assume downdraft velocity to be    !
   ! zero at the bottom, and integrate backwards to get the profile. Now I am including    !
   ! the effect of entrainment and detrainment, in a simplified way, treating kinetic      !
   ! energy as a thermodynamic variable, which means that it has entrainment and detrain-  !
   ! ment, and also one source (namely the buoyancy), and one sink (friction, defined as   !
   ! in Zhang and Fritsch (1986), as - entrainment*w²).                                    !
   !---------------------------------------------------------------------------------------!
   if (comp_down) then
      wdndraft = 0.
      do k=2,klod
         kr     = k + kgoff
         lske2  = wwind(k-1)*wwind(k-1) + sigw(k-1)*sigw(k-1)
         wdndraft = sqrt(max(min_downdraft2,                                               &
                   (wdndraft*wdndraft*(1.+(1.5*mentrd_rate(k-1)-0.5*cdd(k-1))*dzu_cld(k-1))&
                   -(dbyd(k-1)+dbyd(k))*dzu_cld(k-1)-mentrd_rate(k-1)*lske2)               &
                   / (1.-0.5*(mentrd_rate(k-1)+cdd(k-1))*dzu_cld(k-1))))              
      end do
      areadn = max(1.e-5,dnmf * etad_cld(klod) / (rhod_cld(klod) * wdndraft))
   end if
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !   Updraft. We know that the minimum vertical velocity required at the level of free   !
   ! convection is 0, so we can actually estimate a reference updraft velocity             !
   ! based on the expected value.  Let p(w) be the normal probability density function     !
   ! with average wwind(klou) and standard deviation sigw(klou).  Then the expected value  !
   ! for updrafts wupdraft can be defined in a similar way as the average value, but       !
   ! integrated only between wbuoymin and infinity.  With this we can use the mass flux    !
   ! and the density at the level klou to determine the fraction of area of the cloud.     !
   !---------------------------------------------------------------------------------------!
   wupdraft = min(max_wupdraft,max(expected(wbuoymin,wwind(klou),sigw(klou)),min_wupdraft))
   wnorm    = (wbuoymin - wwind(klou)) / sigw(klou)
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      If there is no convective inhibition, then the FC80 method does not do very well !
   ! because the necessary velocities are so low.  Use use the normalized wind-speed       !
   ! method instead.                                                                       !
   !---------------------------------------------------------------------------------------!
   !areaup = min(1.000,max(1.e-5,upmf / (rhou_cld(klou) * wupdraft)))
   areaup  = min(1.000,max(1.e-5,1. - cdf(wnorm)))
   !---------------------------------------------------------------------------------------!



   return
end subroutine grell_draft_area
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!      This method of calculating the updraft area and updraft velocity assumes vertical   !
! winds in the boundary layer where the updraft starts have a normal distribution defined  !
! by mu_w and sigma_w.  The expectancy is found by integrating  p(w)*w from the wbin to    !
! infinity, divided by the integration of p(w) from the w_bmin to infity.  There are cases !
! where some singularities can occur when w_bin >> mu_w, in which case there is not enough !
! mechanical energy in the boundary layer winds to lift a parcel anyway.                   !
!      Once the expected updraft value is calculated, this follows Fritsch and Chappell    !
! (1980) for the area of the updraft.                                                      !
!                                                                                          !
! J. M. Fritsch and C. F. Chappell., 1980: Numerical prediction of convectively driven     !
!    mesoscale pressure systems. part I: Convective parameterization. J. Atmos. Sci., 37,  !
!    1722--1733. doi:10.1175/1520- 0469(1980)037<1722:NPOCDM>2.0.CO;2.                     !
!------------------------------------------------------------------------------------------!
subroutine grell_mb_updraft(mu_w,sigma_w,w_bmin,mflux_up,rho,ex_w_up,area_up)

   use rconstants,only : pi1        & ! intent(in)
                       , sqrttwopi  & ! intent(in)
                       , sqrthalfpi & ! intent(in)
                       , srtwo      & ! intent(in)
                       , lnexp_min  ! ! intent(in)

   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   real(kind=4), intent(in)  :: mu_w     ! Mean of normal vertical winds          [    m/s]
   real(kind=4), intent(in)  :: sigma_w  ! Std. dev. of normal vertical winds     [    m/s]
   real(kind=4), intent(in)  :: w_bmin   ! Min. vert. vel. to achieve bouyancy    [    m/s]
   real(kind=4), intent(in)  :: mflux_up ! Upward mass flux at draft base         [kg/m2/s]
   real(kind=4), intent(in)  :: rho      ! Density of upward mass flux air        [  kg/m3]
   real(kind=4), intent(out) :: ex_w_up  ! Expected value for the updraft         [    m/s]
   real(kind=4), intent(out) :: area_up  ! Fractional area of this cloud updraft  [      -]
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=4)              :: lnterm   ! Partial solution to ex_w_up
   real(kind=4)              :: ewu_1    ! Partial solution to ex_w_up
   real(kind=4)              :: ewu_2    ! Partial solution to ex_w_up
   real(kind=4)              :: ewu_21   ! Partial solution to ex_w_up
   real(kind=4)              :: ewu_22   ! Partial solution to ex_w_up
   real(kind=4)              :: ewu_3    ! Partial solution to ex_w_up
   real(kind=4)              :: ewu_4    ! Partial solution to ex_w_up
   !----- External functions. -------------------------------------------------------------!
   real(kind=4), external    :: errorfun ! Error function
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Terms to find ex_w_up.                                                            !
   !---------------------------------------------------------------------------------------!
   ewu_1  = 0.5*mu_w
   lnterm = max(- (mu_w-w_bmin) *(mu_w-w_bmin) / ( 2.0 * sigma_w * sigma_w ),lnexp_min )

   ewu_21 = ( sigma_w * sigma_w ) * -exp(lnterm) / ( sqrttwopi * sigma_w)
   ewu_22 = sqrthalfpi * mu_w * sigma_w * errorfun( (mu_w-w_bmin) / ( srtwo * sigma_w ) )  &
          / ( sqrttwopi * sigma_w )
   ewu_2  = ewu_21 - ewu_22

   ewu_3  = 0.5
   ewu_4  = 0.5 * errorfun( (mu_w-w_bmin) / (srtwo*sigma_w) )
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Check whether it can become a singularity.  This only happens when ewu_4 = -0.5,  !
   ! but in this case it would be very unlikely to have convection anyway.                 !
   !---------------------------------------------------------------------------------------!
   if ( ewu_4 < - 0.49999 ) then
      ex_w_up = w_bmin
   else
      ex_w_up = (ewu_1 - ewu_2) / (ewu_3 + ewu_4)
   end if
   !---------------------------------------------------------------------------------------!


   !------ Find the area based on the vertical velocity and mass flux. --------------------!
   area_up = mflux_up / (rho*ex_w_up)
   !---------------------------------------------------------------------------------------!

   return
end subroutine grell_mb_updraft
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This subroutine computes some statistical properties on upward mass flux ensemble.   !
! This is currently used only for CATT.                                                    !
!------------------------------------------------------------------------------------------!
subroutine grell_massflx_stats(m1,icld,itest,dti,maxens_dyn,maxens_lsf,maxens_eff          &
                              ,maxens_cap,inv_ensdim,closure_type,ierr_cap,upmf_ens        &
                              ,sgrell1_3d,sgrell2_3d)

   use rconstants  , only : hr_sec        & ! intent(in)
                          , onethird      ! ! intent(in)
   use mem_ensemble, only : ensemble_vars ! ! type
   use mem_scratch_grell, only: &
           kgoff                & ! intent(in) - BRAMS grid offset
          ,mkx                  ! ! intent(in) - # of cloud grid levels

   implicit none
   integer            , intent(in)    :: m1           ! # of BRAMS levels
   integer            , intent(in)    :: icld         ! Cloud type 
   logical            , intent(in)    :: itest        ! Flag for what to fill in
   real               , intent(in)    :: dti          ! Fraction confrq/frqanl
   integer            , intent(in)    :: maxens_cap   ! # of static controls 
   integer            , intent(in)    :: maxens_dyn   ! # of dynamic controls
   integer            , intent(in)    :: maxens_lsf   ! # of large scale forcings
   integer            , intent(in)    :: maxens_eff   ! # of precip. efficiencies
   real               , intent(in)    :: inv_ensdim   ! 1/(maxens_dyn*maxens_lsf*maxens_eff)
   character(len=2)   , intent(in)    :: closure_type ! Which dynamic control
   integer, dimension(maxens_cap), intent(in) :: ierr_cap
   real, dimension(maxens_dyn,maxens_lsf,maxens_eff,maxens_cap), intent(in) :: upmf_ens
   real, dimension(m1), intent(inout) :: sgrell1_3d ! Lots of things, depends on the index
   real, dimension(m1), intent(inout) :: sgrell2_3d ! Lots of things, depends on the index
   
   real, dimension(maxens_cap) :: upmf_ave_cap  ! Average among all maxens_dyn, all closures
   real, dimension(maxens_cap) :: upmf_ave_cap1 ! Same thing but for Grell only
   real, dimension(maxens_cap) :: upmf_ave_cap2 ! Same thing but for Frank-Cohen (1987)
   real, dimension(maxens_cap) :: upmf_ave_cap3 ! Same thing but for Krishnamurti (1983) 
   real, dimension(maxens_cap) :: upmf_ave_cap4 ! Same thing but for Kain-Fritsch (1990)
   real, dimension(maxens_cap) :: upmf_ave_cap5 ! Same thing but for Arakawa-Shubert (1974)
   real, dimension(maxens_dyn) :: upmf_ave      ! Average for each dynamic test
   real, dimension(maxens_dyn) :: upmf_std      ! Std. Deviation for each dynamic test
   real, dimension(maxens_dyn) :: upmf_ske      ! Skewness for each dynamic test
   real, dimension(maxens_dyn) :: upmf_cur      ! Curtosis for each dynamic test
   real                        :: upmf_tot_ave  ! Overall average
   real                        :: upmf_tot_std  ! Overall standard deviation
   real                        :: upmf_tot_ske  ! Overall skewness
   real                        :: upmf_tot_cur  ! Overall curtosis
   
   real, parameter :: small_number = 1.e-6      ! A small number...
   integer         :: iedt,imbp,idync,icap
   real            :: maxens_effdyn_i, maxens_efflsf_i, maxens_eff_i, maxens_dyn_i
   real            :: maxens_effdynlsf_i,maxens_cap_i
   real            :: thisdiff
   
   !---------------------------------------------------------------------------------------!
   !    Two tests to avoid crashing. BTW, these tests should be done at opspec, but I will !
   ! leave them here for now.                                                              !
   !---------------------------------------------------------------------------------------!
   if (m1 < 30)      call abort_run ('nnzp should be always > 30 for Grell+CATT'           &
                                    ,'grell_massflx_stats','grell_cupar_init.f90'          )
   if (maxens_cap /= 3) call abort_run('maxens_cap should be 3 for Grell + CATT!!!'        &
                                       ,'grell_massflx_stats','grell_cupar_init.f90'       )
   !---------------------------------------------------------------------------------------!


   if (all(ierr_cap(:) /= 0)) then
      select case (icld)
      case (1)
         sgrell1_3d = 0.
      case (2)
         sgrell2_3d = 0.
      end select
      return
   end if

   maxens_effdyn_i    = 1./real(maxens_eff*maxens_dyn)
   maxens_efflsf_i    = 1./real(maxens_eff*maxens_lsf)
   maxens_effdynlsf_i = 1./real(maxens_eff*maxens_lsf*maxens_dyn)
   maxens_eff_i    = 1./real(maxens_eff)
   maxens_cap_i    = 1./real(maxens_cap)

   !---------------------------------------------------------------------------------------!
   !    Flushing all variables to zero in case some closures are not in use.               !
   !---------------------------------------------------------------------------------------!
   upmf_ave_cap (:)= 0.
   upmf_ave_cap1(:)= 0.
   upmf_ave_cap2(:)= 0.
   upmf_ave_cap3(:)= 0.
   upmf_ave_cap4(:)= 0.
   upmf_ave_cap5(:)= 0.
   upmf_ave(:)     = 0.
   upmf_std(:)     = 0.
   upmf_ske(:)     = 0.
   upmf_cur(:)     = 0.
   upmf_tot_ave    = 0. 
   upmf_tot_std    = 0. 
   upmf_tot_ske    = 0. 
   upmf_tot_cur    = 0. 

   !---------------------------------------------------------------------------------------!
   !   The mean for each cap_maxs member                                                   !
   !---------------------------------------------------------------------------------------!
   do imbp=1,maxens_cap
      upmf_ave_cap(icap) = sum(upmf_ens(:,:,:,icap)) * maxens_effdynlsf_i
   end do

   !---------------------------------------------------------------------------------------!
   !   The mean for each closure group. Now this part depends on the closure type.         !
   !---------------------------------------------------------------------------------------!
   select case (trim(closure_type))
   case ('gr')
      do icap=1,maxens_cap
         upmf_ave_cap1(icap)= sum(upmf_ens(1,:,:,icap)) * maxens_efflsf_i
      end do
   case ('lo')
      do icap=1,maxens_cap
         upmf_ave_cap2(icap)= sum(upmf_ens(1,:,:,icap)) * maxens_efflsf_i
      end do
   case ('mc')
      do icap=1,maxens_cap
         upmf_ave_cap3(icap)= sum(upmf_ens(1,:,:,icap)) * maxens_efflsf_i
      end do
   case ('kf')
      do icap=1,maxens_cap
         upmf_ave_cap4(icap)= sum(upmf_ens(1,:,:,icap)) * maxens_efflsf_i
      end do
   case ('as')
      do icap=1,maxens_cap
         upmf_ave_cap5(icap)= sum(upmf_ens(1,:,:,icap)) * maxens_efflsf_i
      end do
   case ('nc','en','qi')
      !------------------------------------------------------------------------------------!
      !  This became quite out of order after I changed the order in which each closure    !
      ! is called. I will keep "disorganized" for back compability.                        !
      ! 1-3: Kain-Fritsch     (upmf_ave_cap4)                                              !
      ! 4-7: Arakawa-Schubert (upmf_ave_cap5)                                              !
      ! 8-10: Grell           (upmf_ave_cap1)                                              !
      ! 11-13: Krishnamurti   (upmf_ave_cap3)                                              !
      ! 14-16: Frank-Cohen    (upmf_ave_cap2)                                              !
      !------------------------------------------------------------------------------------!
      do icap=1,maxens_cap
         upmf_ave_cap4(icap) = sum(upmf_ens(1:3,:,:,icap)) * onethird * maxens_efflsf_i
      end do
      do icap=1,maxens_cap
         upmf_ave_cap5(icap) = sum(upmf_ens(4:7,:,:,icap))*0.250*maxens_efflsf_i
      end do
      if (closure_type == 'nc' .or. closure_type == 'en') then
         do icap=1,maxens_cap
            upmf_ave_cap1(icap) = sum(upmf_ens(8:10,:,:,icap)) * onethird * maxens_efflsf_i
         end do
      end if
      if (closure_type == 'en') then
         do icap=1,maxens_cap
            upmf_ave_cap3(icap) = sum(upmf_ens(11:13,:,:,icap)) * onethird*maxens_efflsf_i
         end do
         do imbp=1,maxens_lsf
            upmf_ave_cap2(icap) = sum(upmf_ens(14:16,:,:,icap)) * onethird*maxens_efflsf_i
         end do
      end if
   end select

   !---------------------------------------------------------------------------------------!
   !   Mass flux average for each closure                                                  !
   !---------------------------------------------------------------------------------------!
   do idync=1,maxens_dyn
      upmf_ave(idync) = sum(upmf_ens(:,:,idync,:)) * maxens_efflsf_i
   end do
   
   upmf_tot_ave = sum(upmf_ens) * inv_ensdim
   
   !---------------------------------------------------------------------------------------!
   !   Computing Standard deviation, skewness, and curtosis                                !
   !---------------------------------------------------------------------------------------!
   do icap=1,maxens_cap
      do iedt=1,maxens_eff
         do imbp=1,maxens_lsf
            do idync=1,maxens_dyn
               thisdiff = max(small_number,upmf_ens(idync,imbp,iedt,icap)-upmf_ave(idync))
               upmf_std(idync) = upmf_std(idync) + thisdiff**2
               upmf_ske(idync) = upmf_ske(idync) + thisdiff**3
               upmf_cur(idync) = upmf_cur(idync) + thisdiff**4
               
               thisdiff = max(small_number,upmf_ens(idync,imbp,iedt,icap)-upmf_tot_ave)
               upmf_tot_std    = upmf_tot_std    + thisdiff**2
               upmf_tot_ske    = upmf_tot_ske    + thisdiff**3
               upmf_tot_cur    = upmf_tot_cur    + thisdiff**4
            end do
         end do
      end do
   end do
   !----- Wrapping up the variables, scaling them back ------------------------------------!
   do idync=1,maxens_dyn
      upmf_std(idync) = sqrt(upmf_std(idync) * maxens_dyn_i)
      if (upmf_std(idync) > small_number) then
         upmf_ske(idync) = upmf_ske(idync) * maxens_dyn_i / (upmf_std(idync)**3)
         upmf_cur(idync) = upmf_cur(idync) * maxens_dyn_i / (upmf_std(idync)**4)
      else
         upmf_ske(idync) = 0.
         upmf_cur(idync) = 0.
      end if
   end do
   if (upmf_tot_std > small_number) then
      upmf_tot_ske = upmf_tot_ske * maxens_dyn_i / (upmf_tot_std**3)
      upmf_tot_cur = upmf_tot_cur * maxens_dyn_i / (upmf_tot_std**4)
   else
      upmf_tot_ske = 0.
      upmf_tot_cur = 0.
   end if
   
   !---------------------------------------------------------------------------------------!
   !   Saving the data into the mixed structure                                            !
   !---------------------------------------------------------------------------------------!

   if (icld == 1 .and. itest) then
      
      sgrell1_3d (2) = upmf_tot_ave
      sgrell1_3d (3) = upmf_tot_std
      sgrell1_3d (4) = upmf_tot_ske
      sgrell1_3d (5) = upmf_tot_cur

      select case (closure_type)
      case ('gr')
         sgrell1_3d(6)  = upmf_ave(1)
         sgrell1_3d(7)  = 0.
         sgrell1_3d(8)  = 0.
         sgrell1_3d(9)  = 0.
         sgrell1_3d(10) = 0.
      case ('lo')
         sgrell1_3d(6)  = 0.
         sgrell1_3d(7)  = upmf_ave(1)
         sgrell1_3d(8)  = 0.
         sgrell1_3d(9)  = 0.
         sgrell1_3d(10) = 0.
      case ('mc')
         sgrell1_3d(6)  = 0.
         sgrell1_3d(7)  = 0.
         sgrell1_3d(8)  = upmf_ave(1)
         sgrell1_3d(9)  = 0.
         sgrell1_3d(10) = 0.
      case ('kf')
         sgrell1_3d(6)  = 0.
         sgrell1_3d(7)  = 0.
         sgrell1_3d(8)  = 0.
         sgrell1_3d(9)  = upmf_ave(1)
         sgrell1_3d(10) = 0.
      case ('as')
         sgrell1_3d(6)  = 0.
         sgrell1_3d(7)  = 0.
         sgrell1_3d(8)  = 0.
         sgrell1_3d(9)  = 0.
         sgrell1_3d(10) = upmf_ave(1)
      case ('nc','en','qi')
         sgrell1_3d(9)  = onethird * sum(upmf_ave(1:3))
         sgrell1_3d(10) = .25      * sum(upmf_ave(4:7))
         if (closure_type == 'nc' .or. closure_type == 'en') then
            sgrell1_3d(6)  = onethird * sum(upmf_ave(8:10))
         end if
         if (closure_type == 'en') then
            sgrell1_3d(7) = onethird * sum(upmf_ave(14:16))
            sgrell1_3d(8) = onethird * sum(upmf_ave(11:13))
         else
            sgrell1_3d(7) = 0.
            sgrell1_3d(8) = 0.
         end if
      end select
      do imbp=1,maxens_lsf
         sgrell1_3d(16+(imbp-1)) = upmf_ave_cap(imbp)
      end do
      
      sgrell2_3d(1)  =  upmf_ave_cap1(1)
      sgrell2_3d(2)  =  upmf_ave_cap1(2)
      sgrell2_3d(3)  =  upmf_ave_cap1(3)
      sgrell2_3d(4)  =  upmf_ave_cap2(1)
      sgrell2_3d(5)  =  upmf_ave_cap2(2)
      sgrell2_3d(6)  =  upmf_ave_cap2(3)
      sgrell2_3d(7)  =  upmf_ave_cap3(1)
      sgrell2_3d(8)  =  upmf_ave_cap3(2)
      sgrell2_3d(9)  =  upmf_ave_cap3(3)
      sgrell2_3d(10) =  upmf_ave_cap4(1)
      sgrell2_3d(11) =  upmf_ave_cap4(2)
      sgrell2_3d(12) =  upmf_ave_cap4(3)
      sgrell2_3d(13) =  upmf_ave_cap5(1)
      sgrell2_3d(14) =  upmf_ave_cap5(2)
      sgrell2_3d(15) =  upmf_ave_cap5(3)

   elseif (icld == 1) then
      select case (closure_type)
      case ('gr')
         sgrell1_3d(11) = upmf_ave(1)*hr_sec
         sgrell1_3d(12) = 0.
         sgrell1_3d(13) = 0.
         sgrell1_3d(14) = 0.
         sgrell1_3d(15) = 0.
      case ('lo')
         sgrell1_3d(11) = 0.
         sgrell1_3d(12) = upmf_ave(1)*hr_sec
         sgrell1_3d(13) = 0.
         sgrell1_3d(14) = 0.
         sgrell1_3d(15) = 0.
      case ('mc')
         sgrell1_3d(11) = 0.
         sgrell1_3d(12) = 0.
         sgrell1_3d(13) = upmf_ave(1)*hr_sec
         sgrell1_3d(14) = 0.
         sgrell1_3d(15) = 0.
      case ('kf')
         sgrell1_3d(11) = 0.
         sgrell1_3d(12) = 0.
         sgrell1_3d(13) = 0.
         sgrell1_3d(14) = upmf_ave(1)*hr_sec
         sgrell1_3d(15) = 0.
      case ('as')
         sgrell1_3d(11) = 0.
         sgrell1_3d(12) = 0.
         sgrell1_3d(13) = 0.
         sgrell1_3d(14) = 0.
         sgrell1_3d(15) = upmf_ave(1)*hr_sec
      case ('nc','en','qi')
         sgrell1_3d(14) = onethird * sum(upmf_ave(1:3)) *hr_sec
         sgrell1_3d(15) = .25      * sum(upmf_ave(4:7)) *hr_sec
         if (closure_type == 'nc' .or. closure_type == 'en') then
            sgrell1_3d(11) = onethird * sum(upmf_ave(8:10))*hr_sec
         end if
         if (closure_type == 'en') then
            sgrell1_3d(12)= onethird * sum(upmf_ave(14:16))*hr_sec
            sgrell1_3d(13)= onethird * sum(upmf_ave(11:13))*hr_sec
         else
            sgrell1_3d(12)= 0.
            sgrell1_3d(13)= 0.
         end if
      end select
      
      !----- Convert from kg/m2/s to kg/m2/hour -------------------------------------------!
      sgrell2_3d(1)   = upmf_ave_cap1(1)*sgrell2_3d(1 )*hr_sec
      sgrell2_3d(2)   = upmf_ave_cap1(2)*sgrell2_3d(2 )*hr_sec
      sgrell2_3d(3)   = upmf_ave_cap1(3)*sgrell2_3d(3 )*hr_sec
      sgrell2_3d(4)   = upmf_ave_cap2(1)*sgrell2_3d(4 )*hr_sec
      sgrell2_3d(5)   = upmf_ave_cap2(2)*sgrell2_3d(5 )*hr_sec
      sgrell2_3d(6)   = upmf_ave_cap2(3)*sgrell2_3d(6 )*hr_sec
      sgrell2_3d(7)   = upmf_ave_cap3(1)*sgrell2_3d(7 )*hr_sec
      sgrell2_3d(8)   = upmf_ave_cap3(2)*sgrell2_3d(8 )*hr_sec
      sgrell2_3d(9)   = upmf_ave_cap3(3)*sgrell2_3d(9 )*hr_sec
      sgrell2_3d(10)  = upmf_ave_cap4(1)*sgrell2_3d(10)*hr_sec
      sgrell2_3d(11)  = upmf_ave_cap4(2)*sgrell2_3d(11)*hr_sec
      sgrell2_3d(12)  = upmf_ave_cap4(3)*sgrell2_3d(12)*hr_sec
      sgrell2_3d(13)  = upmf_ave_cap5(1)*sgrell2_3d(13)*hr_sec
      sgrell2_3d(14)  = upmf_ave_cap5(2)*sgrell2_3d(14)*hr_sec
      sgrell2_3d(15)  = upmf_ave_cap5(3)*sgrell2_3d(15)*hr_sec
      sgrell2_3d(16)  = sgrell2_3d(16) + sgrell2_3d(1 )*dti
      sgrell2_3d(17)  = sgrell2_3d(17) + sgrell2_3d(2 )*dti
      sgrell2_3d(18)  = sgrell2_3d(18) + sgrell2_3d(3 )*dti
      sgrell2_3d(19)  = sgrell2_3d(19) + sgrell2_3d(4 )*dti
      sgrell2_3d(20)  = sgrell2_3d(20) + sgrell2_3d(5 )*dti
      sgrell2_3d(21)  = sgrell2_3d(21) + sgrell2_3d(6 )*dti
      sgrell2_3d(22)  = sgrell2_3d(22) + sgrell2_3d(7 )*dti
      sgrell2_3d(23)  = sgrell2_3d(23) + sgrell2_3d(8 )*dti
      sgrell2_3d(24)  = sgrell2_3d(24) + sgrell2_3d(9 )*dti
      sgrell2_3d(25)  = sgrell2_3d(25) + sgrell2_3d(10)*dti
      sgrell2_3d(26)  = sgrell2_3d(26) + sgrell2_3d(11)*dti
      sgrell2_3d(27)  = sgrell2_3d(27) + sgrell2_3d(12)*dti
      sgrell2_3d(28)  = sgrell2_3d(28) + sgrell2_3d(13)*dti
      sgrell2_3d(29)  = sgrell2_3d(29) + sgrell2_3d(14)*dti
      sgrell2_3d(30)  = sgrell2_3d(30) + sgrell2_3d(15)*dti
     
   else ! (icld == 2)
      select case (closure_type)
      case ('gr')
         sgrell1_3d(21) = upmf_ave(1)
         sgrell1_3d(22) = 0.
         sgrell1_3d(23) = 0.
      case ('kf')
         sgrell1_3d(21) = 0.
         sgrell1_3d(22) = 0.
         sgrell1_3d(23) = upmf_ave(1)
      case ('as')
         sgrell1_3d(21) = 0.
         sgrell1_3d(22) = upmf_ave(1)
         sgrell1_3d(23) = 0.
      case ('nc','qi')
         sgrell1_3d(23) = onethird * sum(upmf_ave(1:3))*hr_sec
         sgrell1_3d(22) = .25      * sum(upmf_ave(4:7))*hr_sec
         if (closure_type == 'nc') then
            sgrell1_3d(21) = onethird * sum(upmf_ave(8:10))*hr_sec
         end if
      end select
   end if
   return
end subroutine grell_massflx_stats
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine checks whether the partial results in the cumulus parametrisation   !
! make sense or not.                                                                       !
!------------------------------------------------------------------------------------------!
subroutine grell_sanity_check(mkx,mgmzp,z,press,exner,theiv,thil,t,qtot,qvap,qliq,qice     &
                             ,co2,rho,which)
   use grell_coms, only : grellmax_zcheck & ! intent(in)
                        , grell_lapse_wet & ! intent(in)
                        , grellmin_t0     & ! intent(in)
                        , grellmax_t0     & ! intent(in)
                        , grellmin_rhv    & ! intent(in)
                        , grellmax_rhv    & ! intent(in)
                        , grellmin_co2    & ! intent(in)
                        , grellmax_co2    ! ! intent(in)
   use therm_lib , only : extemp2theta    & ! intent(in)
                        , eslif           & ! intent(in)
                        , thetaeivs       ! ! intent(in)
   use rconstants, only : ep              & ! intent(in)
                        , t00             & ! intent(in)
                        , gocp            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer                           , intent(in) :: mkx   ! Number of layers
   integer                           , intent(in) :: mgmzp ! Array dimensions
   real            , dimension(mgmzp), intent(in) :: z     ! Height              [       m]
   real            , dimension(mgmzp), intent(in) :: press ! Pressure            [      Pa]
   real            , dimension(mgmzp), intent(in) :: exner ! Exner function      [  J/kg/K]
   real            , dimension(mgmzp), intent(in) :: theiv ! Equiv. pot. temp.   [       K]
   real            , dimension(mgmzp), intent(in) :: thil  ! I.L. Pot. temp.     [       K]
   real            , dimension(mgmzp), intent(in) :: t     ! Temperature         [       K]
   real            , dimension(mgmzp), intent(in) :: qtot  ! Total mixing ratio  [   kg/kg]
   real            , dimension(mgmzp), intent(in) :: qvap  ! Vapour mixing ratio [   kg/kg]
   real            , dimension(mgmzp), intent(in) :: qliq  ! Liquid mixing ratio [   kg/kg]
   real            , dimension(mgmzp), intent(in) :: qice  ! Ice mixing ratio    [   kg/kg]
   real            , dimension(mgmzp), intent(in) :: co2   ! CO2 mixing ratio    [µmol/mol]
   real            , dimension(mgmzp), intent(in) :: rho   ! Air density         [   kg/m3]
   character(len=*)                  , intent(in) :: which ! Which call?
   !----- Local variables. ----------------------------------------------------------------!
   real    :: grellmin_t      ! Minimum temperature                              [       K]
   real    :: grellmax_t      ! Maximum temperature                              [       K]
   real    :: grellmin_thil   ! Minimum ice-liquid potential temperature         [       K]
   real    :: grellmax_thil   ! Maximum ice-liquid potential temperature         [       K]
   real    :: grellmin_theiv  ! Minimum ice-vapour equivalent pot. temperature   [       K]
   real    :: grellmax_theiv  ! Maximum ice-vapour equivalent pot. temperature   [       K]
   real    :: grellmin_pvap   ! Minimum vapour pressure                          [      Pa]
   real    :: grellmax_pvap   ! Maximum vapour pressure                          [      Pa]
   real    :: grellmin_qvap   ! Minimum vapour mixing ratio                      [   kg/kg]
   real    :: grellmax_qvap   ! Maximum vapour mixing ratio                      [   kg/kg]
   real    :: grellmin_qtot   ! Minimum total mixing ratio                       [   kg/kg]
   real    :: grellmax_qtot   ! Maximum total mixing ratio                       [   kg/kg]
   logical :: everything_fine ! This will become false if anything looks wrong   [     T|F]
   integer :: k               ! Counter                                          [      --]
   integer :: l               ! Counter                                          [      --]
   integer :: m               ! Counter                                          [      --]
   !---------------------------------------------------------------------------------------!

   

   !---------------------------------------------------------------------------------------!
   !      Let's be optimistic and assume that everything is fine.                          !
   !---------------------------------------------------------------------------------------!
   everything_fine = .true.
   !---------------------------------------------------------------------------------------!

   !----- Loop over all layers and check for problems. ------------------------------------!
   checkloop: do k = 1,mkx
      if (z(k) > grellmax_zcheck) exit checkloop

      !------------------------------------------------------------------------------------!
      !     Find derived bounds.                                                           !
      !------------------------------------------------------------------------------------!
      !----- Temperature. -----------------------------------------------------------------!
      grellmin_t      = max(t00-120.,grellmin_t0 - gocp * z(k))
      grellmax_t      = grellmax_t0 - grell_lapse_wet * z(k)
      everything_fine = t(k) >= grellmin_t .and. t(k) <= grellmax_t
      !----- Vapour pressure. -------------------------------------------------------------!
      grellmin_pvap   = grellmin_rhv * eslif(grellmin_t)
      grellmax_pvap   = grellmax_rhv * eslif(grellmax_t)
      !----- Vapour mixing ratio. ---------------------------------------------------------!
      grellmin_qvap   = ep * grellmin_pvap / (press(k) - grellmin_pvap)
      grellmax_qvap   = ep * grellmax_pvap / (press(k) - grellmax_pvap)
      everything_fine = qvap(k) >= grellmin_qvap .and. qvap(k) <= grellmax_qvap
      !----- Total mixing ratio.  Minimum is unsaturated, maximum assumed to 2*qvap... ----!
      grellmin_qtot   = grellmin_qvap
      grellmax_qtot   = grellmax_qvap * 2.0
      everything_fine = qtot(k) >= grellmin_qtot .and. qtot(k) <= grellmax_qtot
      !----- Ice-liquid potential temperature. --------------------------------------------!
      grellmin_thil   = extemp2theta(exner(k),grellmin_t)
      grellmax_thil   = extemp2theta(exner(k),grellmax_t)
      everything_fine = thil(k) >= grellmin_thil .and. thil(k) <= grellmax_thil
      !----- Ice-vapour equivalent potential temperature. ---------------------------------!
      grellmin_theiv  = grellmin_thil
      grellmax_theiv  = thetaeivs(grellmax_thil,grellmax_t,grellmax_qvap,grellmax_qvap,0.)
      everything_fine = theiv(k) >= grellmin_theiv .and. theiv(k) <= grellmax_theiv
      !----- CO2 mixing ratio. ------------------------------------------------------------!
      everything_fine = co2(k) >= grellmin_co2 .and. co2(k) <= grellmax_co2
      !------------------------------------------------------------------------------------!

      if (.not. everything_fine) then
         !---------------------------------------------------------------------------------!
         !      This is the problem of being optimistic: more often than not you may be    !
         ! disappointed... But hey, don't hate me messenger, I'm just going to print       !
         ! this so you can work towards a better model.                                    !
         !---------------------------------------------------------------------------------!
         write (unit=*,fmt='(169a)'    ) ('-',m=1,169)
         write (unit=*,fmt='(3(a,1x))' ) ' -> Event: ',trim(which)                         &
                                        ,' has unrealistic thermodynamics...'
         write (unit=*,fmt='(169a)'    ) ('-',m=1,169)
         write (unit=*,fmt='(a)'       ) ' BOUNDS'
         write (unit=*,fmt='(a)'       ) ''
         write (unit=*,fmt='(a,es12.5)') ' Min TEMP  [    degC]: ',grellmin_t     - t00
         write (unit=*,fmt='(a,es12.5)') ' Max TEMP  [    degC]: ',grellmax_t     - t00
         write (unit=*,fmt='(a,es12.5)') ' Min THIL  [       K]: ',grellmin_thil
         write (unit=*,fmt='(a,es12.5)') ' Max THIL  [       K]: ',grellmax_thil
         write (unit=*,fmt='(a,es12.5)') ' Min THEIV [       K]: ',grellmin_theiv
         write (unit=*,fmt='(a,es12.5)') ' Max THEIV [       K]: ',grellmax_theiv
         write (unit=*,fmt='(a,es12.5)') ' Min PVAP  [     hPa]: ',grellmin_pvap  * 0.01
         write (unit=*,fmt='(a,es12.5)') ' Max PVAP  [     hPa]: ',grellmax_pvap  * 0.01
         write (unit=*,fmt='(a,es12.5)') ' Min QVAP  [    g/kg]: ',grellmin_qvap  * 1000.
         write (unit=*,fmt='(a,es12.5)') ' Max QVAP  [    g/kg]: ',grellmax_qvap  * 1000.
         write (unit=*,fmt='(a,es12.5)') ' Min QTOT  [    g/kg]: ',grellmin_qtot  * 1000.
         write (unit=*,fmt='(a,es12.5)') ' Max QTOT  [    g/kg]: ',grellmax_qtot  * 1000.
         write (unit=*,fmt='(a,es12.5)') ' Min CO2   [umol/mol]: ',grellmin_co2
         write (unit=*,fmt='(a,es12.5)') ' Max CO2   [umol/mol]: ',grellmax_co2
         write (unit=*,fmt='(a)'       ) ''
         write (unit=*,fmt='(169a)'    ) ('-',m=1,169)
         write (unit=*,fmt='(13(a,1x))') '       LAYER','      HEIGHT','    PRESSURE'      &
                                        ,'       EXNER','        TEMP','    THETA_IL'      &
                                        ,'  THETAT_EIV','        QVAP','        QLIQ'      &
                                        ,'        QICE','        QTOT','     DENSITY'      &
                                        ,'         CO2'
         write (unit=*,fmt='(13(a,1x))') '         ---','           m','         hPa'      &
                                        ,'      J/kg/K','        degC','           K'      &
                                        ,'           K','        g/kg','        g/kg'      &
                                        ,'        g/kg','        g/kg','       kg/m3'      &
                                        ,'    umol/mol'
         write (unit=*,fmt='(169a)'    ) ('-',m=1,169)
          do l=mkx,1,-1
            write (unit=*,fmt='(i12,1x,12(f12.2))')                                        &
                                         l,z(l),press(l)*0.01,exner(l),t(l)-t00,thil(l)    &
                                        ,theiv(l),qvap(l)*1000.,qliq(l)*1000.              &
                                        ,qice(l)*1000.,qtot(l)*1000.,rho(l),co2(l)
         end do
         write (unit=*,fmt='(169a)'    ) ('-',m=1,169)
         write (unit=*,fmt='(a)'       ) ''
         call brams_fail_whale()
         call abort_run('Unreasonable thermodynamic variables','grell_sanity_check'        &
                       ,'grell_cupar_aux.f90')
         !---------------------------------------------------------------------------------!
      end if
      !------------------------------------------------------------------------------------!
   end do checkloop
   !---------------------------------------------------------------------------------------!

   return
end subroutine grell_sanity_check
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!     This sub-routine checks whether the partial results in the cumulus parametrisation   !
! make sense or not.                                                                       !
!------------------------------------------------------------------------------------------!
subroutine grell_sanity_thil2tqall(k,z,thil,exner,press,qtot,which)
   use grell_coms, only : grellmax_zcheck & ! intent(in)
                        , grell_lapse_wet & ! intent(in)
                        , grellmin_t0     & ! intent(in)
                        , grellmax_t0     & ! intent(in)
                        , grellmin_rhv    & ! intent(in)
                        , grellmax_rhv    & ! intent(in)
                        , grellmin_co2    & ! intent(in)
                        , grellmax_co2    ! ! intent(in)
   use therm_lib , only : extemp2theta    & ! intent(in)
                        , eslif           & ! intent(in)
                        , thetaeivs       ! ! intent(in)
   use rconstants, only : ep              & ! intent(in)
                        , t00             & ! intent(in)
                        , gocp            ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer         , intent(in) :: k     ! Counter                               [      --]
   real            , intent(in) :: z     ! Height                                [       m]
   real            , intent(in) :: thil  ! I.L. Pot. temp.                       [       K]
   real            , intent(in) :: exner ! Exner function                        [  J/kg/K]
   real            , intent(in) :: press ! Pressure                              [      Pa]
   real            , intent(in) :: qtot  ! Total mixing ratio                    [   kg/kg]
   character(len=*), intent(in) :: which ! Which call?
   !----- Local variables. ----------------------------------------------------------------!
   real    :: grellmin_t      ! Minimum temperature                              [       K]
   real    :: grellmax_t      ! Maximum temperature                              [       K]
   real    :: grellmin_thil   ! Minimum ice-liquid potential temperature         [       K]
   real    :: grellmax_thil   ! Maximum ice-liquid potential temperature         [       K]
   real    :: grellmin_pvap   ! Minimum vapour pressure                          [      Pa]
   real    :: grellmax_pvap   ! Maximum vapour pressure                          [      Pa]
   real    :: grellmin_qvap   ! Minimum vapour mixing ratio                      [   kg/kg]
   real    :: grellmax_qvap   ! Maximum vapour mixing ratio                      [   kg/kg]
   real    :: grellmin_qtot   ! Minimum total mixing ratio                       [   kg/kg]
   real    :: grellmax_qtot   ! Maximum total mixing ratio                       [   kg/kg]
   logical :: everything_fine ! This will become false if anything looks wrong   [     T|F]
   integer :: m               ! Counter 
   !---------------------------------------------------------------------------------------!

   !----- Exiting this subroutine because the bounds may be too small... ------------------!
   return
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !      Let's be optimistic and assume that everything is fine.                          !
   !---------------------------------------------------------------------------------------!
   everything_fine = .true.
   !---------------------------------------------------------------------------------------!


   !---------------------------------------------------------------------------------------!
   !     Find derived bounds.                                                              !
   !---------------------------------------------------------------------------------------!
   !----- Temperature. --------------------------------------------------------------------!
   grellmin_t      = max(t00-120.,grellmin_t0 - gocp * z)
   grellmax_t      = grellmax_t0 - grell_lapse_wet * z
   !----- Vapour pressure. ----------------------------------------------------------------!
   grellmin_pvap   = grellmin_rhv * eslif(grellmin_t)
   grellmax_pvap   = grellmax_rhv * eslif(grellmax_t)
   !----- Vapour mixing ratio. ------------------------------------------------------------!
   grellmin_qvap   = ep * grellmin_pvap / (press - grellmin_pvap)
   grellmax_qvap   = ep * grellmax_pvap / (press - grellmax_pvap)
   !----- Total mixing ratio.  Minimum is unsaturated, maximum assumed to 2*qvap... -------!
   grellmin_qtot   = grellmin_qvap
   grellmax_qtot   = grellmax_qvap * 2.0
   everything_fine = qtot >= grellmin_qtot .and. qtot <= grellmax_qtot
   !----- Ice-liquid potential temperature. -----------------------------------------------!
   grellmin_thil   = extemp2theta(exner,grellmin_t)
   grellmax_thil   = extemp2theta(exner,grellmax_t)
   everything_fine = thil >= grellmin_thil .and. thil <= grellmax_thil
   !---------------------------------------------------------------------------------------!

   if (.not. everything_fine) then
      !------------------------------------------------------------------------------------!
      !      This is the problem of being optimistic: more often than not you may be       !
      ! disappointed... But hey, don't hate me messenger, I'm just going to print          !
      ! this so you can work towards a better model.                                       !
      !------------------------------------------------------------------------------------!
      write (unit=*,fmt='(92a)'       ) ('-',m=1,92)
      write (unit=*,fmt='(3(a,1x))'   ) ' -> Event: ',trim(which)                          &
                                       ,' has unrealistic thermodynamics...'
      write (unit=*,fmt='(92a)'       ) ('-',m=1,92)
      write (unit=*,fmt='(a)'         ) ' BOUNDS'
      write (unit=*,fmt='(a)'         ) ''
      write (unit=*,fmt='(a,es12.5)'  ) '  - Min TEMP  [    degC]: ',grellmin_t     - t00
      write (unit=*,fmt='(a,es12.5)'  ) '  - Max TEMP  [    degC]: ',grellmax_t     - t00
      write (unit=*,fmt='(a,es12.5)'  ) '  - Min THIL  [       K]: ',grellmin_thil
      write (unit=*,fmt='(a,es12.5)'  ) '  - Max THIL  [       K]: ',grellmax_thil
      write (unit=*,fmt='(a,es12.5)'  ) '  - Min PVAP  [     hPa]: ',grellmin_pvap  * 0.01
      write (unit=*,fmt='(a,es12.5)'  ) '  - Max PVAP  [     hPa]: ',grellmax_pvap  * 0.01
      write (unit=*,fmt='(a,es12.5)'  ) '  - Min QVAP  [    g/kg]: ',grellmin_qvap  * 1000.
      write (unit=*,fmt='(a,es12.5)'  ) '  - Max QVAP  [    g/kg]: ',grellmax_qvap  * 1000.
      write (unit=*,fmt='(a,es12.5)'  ) '  - Min QTOT  [    g/kg]: ',grellmin_qtot  * 1000.
      write (unit=*,fmt='(a,es12.5)'  ) '  - Max QTOT  [    g/kg]: ',grellmax_qtot  * 1000.
      write (unit=*,fmt='(a)'         ) ''
      write (unit=*,fmt='(a)'         ) ' VALUES'
      write (unit=*,fmt='(a)'         ) ''
      write (unit=*,fmt='(a,1x,i6)'   ) '  - K                   : ',k
      write (unit=*,fmt='(a,1x,f12.2)') '  - THIL      [       K]: ',thil
      write (unit=*,fmt='(a,1x,f12.2)') '  - EXNER     [  J/Kg/K]: ',exner
      write (unit=*,fmt='(a,1x,f12.2)') '  - PRESS     [     hPa]: ',press * 0.01
      write (unit=*,fmt='(a,1x,f12.2)') '  - QTOT      [    g/kg]: ',qtot  * 1000.
      write (unit=*,fmt='(92a)'       ) ('-',m=1,92)
      write (unit=*,fmt='(a)'         ) ''
      call brams_fail_whale()
      call abort_run('Unreasonable thermodynamic variables','grell_sanity_thil2tqall'      &
                    ,'grell_cupar_aux.f90')
      !------------------------------------------------------------------------------------!
   end if
   !---------------------------------------------------------------------------------------!

   return
end subroutine grell_sanity_thil2tqall
!==========================================================================================!
!==========================================================================================!
