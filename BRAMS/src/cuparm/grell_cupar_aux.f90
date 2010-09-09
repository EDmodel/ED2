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
subroutine initial_tend_grell(m1,m2,m3,i,j,tht,tket,rtt,co2t)
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
   integer, intent(in)                      :: m1,m2,m3   ! Grid dimensions
   integer, intent(in)                      :: i,j        ! Current position
   real   , intent(in), dimension(m1,m2,m3) :: tht        ! Potential temperature tend.
   real   , intent(in), dimension(m1,m2,m3) :: tket       ! Turbulent Kinetic Energy tend.
   real   , intent(in), dimension(m1,m2,m3) :: co2t       ! Total CO2 mixing ratio tend.
   real   , intent(in), dimension(m1,m2,m3) :: rtt        ! Total H2O mixing ratio tend.
   !------ Local variables ----------------------------------------------------------------!
   integer                                  :: k,kr       ! Counters
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
      dthildt(k)      = tht(kr,i,j)
      dqtotdt(k)      = rtt(kr,i,j)
      dtkedt(k)       = tket(kr,i,j)
      dco2dt(k)       = co2t(kr,i,j)

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
subroutine initial_thermo_grell(m1,dtime,thp,theta,rtp,co2p,pi0,pp,pc,wp,dn0,tkep,rliq     &
                               ,rice,wstd)

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
   use rconstants, only : cp      & ! intent(in)
                        , cpi     & ! intent(in)
                        , cpor    & ! intent(in)
                        , p00     & ! intent(in)
                        , grav    & ! intent(in)
                        , rdry    & ! intent(in)
                        , epi     & ! intent(in)
                        , toodry  & ! intent(in)
                        , sigwmin & ! intent(in)
                        , tkmin   ! ! intent(in)

   !------ External functions -------------------------------------------------------------!
   use therm_lib, only: &
           rslif        & ! Function that finds the saturation mixing ratio
          ,thetaeiv     & ! Function that finds Thetae_iv  
          ,thil2temp    & ! Function that gives temperature from theta_il, rliq and rice.
          ,idealdens    ! ! Function that gives the density for ideal gasses

   implicit none
   !------ I/O variables ------------------------------------------------------------------!
   integer, intent(in)                  :: m1    ! Grid dimension                 [    ---]
   real   , intent(in)                  :: dtime ! Time step                      [      s]
   real   , intent(in)  , dimension(m1) :: thp   ! Ice-liquid potential temp.     [      K]
   real   , intent(in)  , dimension(m1) :: theta ! Potential temperature          [      K]
   real   , intent(in)  , dimension(m1) :: rtp   ! Total H2O mixing ratio         [  kg/kg]
   real   , intent(in)  , dimension(m1) :: co2p  ! Total CO2 mixing ratio         [    ppm]
   real   , intent(in)  , dimension(m1) :: pi0   ! Reference Exner function       [ J/kg/K]
   real   , intent(in)  , dimension(m1) :: pp    ! Current perturbation on pi     [ J/kg/K]
   real   , intent(in)  , dimension(m1) :: pc    ! Future perturbation on pi      [ J/kg/K]
   real   , intent(in)  , dimension(m1) :: dn0   ! Reference density              [  kg/m³]
   real   , intent(in)  , dimension(m1) :: wp    ! Vertical velocity              [    m/s]
   real   , intent(in)  , dimension(m1) :: tkep  ! Turbulent kinetic energy       [   J/kg]
   real   , intent(in)  , dimension(m1) :: rliq  ! Liquid water mixing ratio      [  kg/kg]
   real   , intent(in)  , dimension(m1) :: rice  ! Ice mixing ratio               [  kg/kg]
   real   , intent(in)  , dimension(m1) :: wstd  ! Standard deviation of wp       [    m/s]
   !------ Local variables ----------------------------------------------------------------!
   integer                              :: k,kr  ! Counters
   real                                 :: dq    ! Diff. on vapour mixing ratio   [  kg/kg]
   real                                 :: qsat  ! Sat. mixing ratio (scratch)    [  kg/kg]
   !---------------------------------------------------------------------------------------!

   do k=1,mkx
      kr=k+kgoff
      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      !    Finding the current state variables, including the effect of shallower cumulus  !
      ! if that is the case.                                                               !
      !------------------------------------------------------------------------------------!
      !------ 1. Exner function -----------------------------------------------------------!
      exner0(k) = pi0(kr)   + pp(kr)
      !------ 2. Pressure. ----------------------------------------------------------------!
      p0(k)     = p00*(cpi*exner0(k))**cpor
      !------------------------------------------------------------------------------------!
      ! 3. Temperature and water.                                                          !
      !------------------------------------------------------------------------------------!
      thil0(k)  = thp(kr)
      qtot0(k)  = max(toodry,rtp(kr))
      qliq0(k)  = max(0.,rliq(kr))
      qice0(k)  = max(0.,rice(kr))
      qvap0(k)  = max(toodry,qtot0(k)-qice0(k)-qliq0(k))
      t0(k)     = cpi * theta(kr) * exner0(k)

      !------ 4. Finding the ice-vapour equivalent potential temperature ------------------!
      theiv0(k) = thetaeiv(thil0(k),p0(k),t0(k),qvap0(k),qtot0(k),5)

      !------ 5. CO2 mixing ratio. --------------------------------------------------------!
      co20(k)   = co2p(kr)
      !------ 6. Turbulent kinetic energy [m²/s²] -----------------------------------------!
      tke0(k)     = tkep(kr)
      !------ 7. Vertical velocity in terms of pressure, or Lagrangian dp/dt [ Pa/s] ------!
      omeg(k)     = -grav*dn0(kr)*.5*( wp(kr)+wp(kr-1) )
      !------ 8. Vertical velocity [m/s], this is staggered, averaging... -----------------!
      wwind(k)    = 0.5 * (wp(kr)+wp(kr-1))
      !------ 9. Standard-deviation of vertical velocity ----------------------------------!
      sigw(k)     = max(wstd(kr),sigwmin)
      !------ 10. Air density -------------------------------------------------------------!
      rho(k)   = idealdens(p0(k),t0(k),qvap0(k),qtot0(k))
      !------------------------------------------------------------------------------------!
      !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!



      ![[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[!
      !------------------------------------------------------------------------------------!
      !     Finding what the state variables will be in the next time, assuming no convec- !
      ! tion at this point (we will call these forced variables).  Most variables will be  !
      ! updated using the tendency, except for the Exner function and diagnostic vari-     !
      !ables.                                                                              !
      !------------------------------------------------------------------------------------!
      !------ 1. Exner function pc is the future Exner perturbation. ----------------------!
      exner(k) = pi0(kr)   + pc(kr)
      !------ 2. Pressure -----------------------------------------------------------------!
      p(k)     = p00*(cpi*exner(k))**cpor
      !------ 3. Ice-liquid potential temperature, with the tendency ----------------------!
      thil(k)  = thp(kr) + dthildt(k)*dtime
      !------ 4. Total mixing ratio, with the tendency ------------------------------------!
      qtot(k)  = max(toodry,rtp(kr)   + dqtotdt(k) * dtime)
      !------ 5. Finding the equilibrium state. Temperature 1st guess is simply t0 --------!
      t(k)     = t0(k)
      call thil2tqall(thil(k),exner(k),p(k),qtot(k),qliq(k),qice(k),t(k),qvap(k),qsat)
      !------ 6. Finding the ice-vapour equivalent potential temperature ------------------!
      theiv(k) = thetaeiv(thil(k),p(k),t(k),qvap(k),qtot(k),6)
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
   !     Computing the surface variables. The only one that will be truly different is the !
   ! Exner function (and consequently, pressure). The other values will be based on the    !
   ! level above. This is going to be just a boundary condition, so they will directly     !
   ! affect the parametrisation.                                                           !
   !---------------------------------------------------------------------------------------!
   !----- 1. Exner function ---------------------------------------------------------------!
   exnersur    = sqrt((pp(lpw-1)+pi0(lpw-1))*(pp(lpw)+pi0(lpw)))
   !----- 2. Pressure ---------------------------------------------------------------------!
   psur        = p00*(cpi*exnersur)**cpor
   !----- 3. Ice liquid potential temperature ---------------------------------------------!
   thilsur     = thp(lpw)
   !----- 4. Total mixing ratio -----------------------------------------------------------!
   qtotsur     = max(toodry,rtp(lpw))
   !----- 5. Liquid water mixing ratio ----------------------------------------------------!
   qliqsur     = max(0.,rliq(lpw))
   !----- 6. Ice mixing ratio -------------------------------------------------------------!
   qicesur     = max(0.,rice(lpw))
   !----- 7. Vapour mixing ratio ----------------------------------------------------------!
   qvapsur     = max(toodry,qtotsur-qliqsur-qicesur)
   !----- 7. Temperature ------------------------------------------------------------------!
   tsur        = cpi*theta(lpw)*exnersur
   !----- 8. Ice-vapour equivalent potential temperature ----------------------------------!
   theivsur    = thetaeiv(thilsur,psur,tsur,qvapsur,qtotsur,7)
   !----- 9. CO2 mixing ratio -------------------------------------------------------------!
   co2sur      = co2p(lpw)
   !---------------------------------------------------------------------------------------!
   !]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]!



   !---------------------------------------------------------------------------------------!
   !     Finding the integrated moisture convergence. This is done outside the loop so we  !
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
   real                         , intent(inout) :: prev_dnmf ! Previous downdraft
   !------ Local variables ----------------------------------------------------------------!
   integer            :: k         ! Counter for current Grell level
   integer            :: kr        ! Counter for corresponding BRAMS level
   !---------------------------------------------------------------------------------------!

   
   !------ Initializing scalars -----------------------------------------------------------!
   prev_dnmf = last_dnmf(i,j)
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
   real, parameter                     :: min_updraft         = 0.60
   real, parameter                     :: min_updraft2        = min_updraft*min_updraft
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
   !   Updraft. We know that the minimum vertical velocity required at the level of free   !
   ! convection is 0, so we can actually estimate a reference updraft velocity             !
   ! based on the expected value.  Let p(w) be the normal probability density function     !
   ! with average wwind(klou) and standard deviation sigw(klou).  Then the expected value  !
   ! for updrafts wupdraft can be defined in a similar way as the average value, but       !
   ! integrated only between wbuoymin and infinity.  With this we can use the mass flux    !
   ! and the density at the level klou to determine the fraction of area of the cloud.     !
   !---------------------------------------------------------------------------------------!
   wupdraft = expected(wbuoymin,wwind(klou),sigw(klou))
   wnorm    = (wbuoymin - wwind(klou)) / sigw(klou)
   areaup   = max(1.e-5,1. - cdf(wnorm))

   return
end subroutine grell_draft_area
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

   use rconstants  , only : hr_sec
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
   case ('nc','en')
      !------------------------------------------------------------------------------------!
      !  This became quite out of order after I changed the order in which each closure    !
      ! is called. I will keep "disorganized" for back compability.                        !
      ! 1-3: Grell            (upmf_ave_cap1)                                              !
      ! 4-7: Arakawa-Schubert (upmf_ave_cap5)                                              !
      ! 8-10: Kain-Fritsch    (upmf_ave_cap4)                                              !
      ! 11-13: Krishnamurti   (upmf_ave_cap3)                                              !
      ! 14-16: Frank-Cohen    (upmf_ave_cap2)                                              !
      !------------------------------------------------------------------------------------!
      do icap=1,maxens_cap
         upmf_ave_cap1(icap) = sum(upmf_ens(1:3,:,:,icap))*0.3333333333*maxens_efflsf_i
      end do
      do icap=1,maxens_cap
         upmf_ave_cap5(icap) = sum(upmf_ens(4:7,:,:,icap))*0.250*maxens_efflsf_i
      end do
      do icap=1,maxens_cap
         upmf_ave_cap4(icap) = sum(upmf_ens(8:10,:,:,icap))                                &
                             * 0.3333333333 * maxens_efflsf_i
      end do
      if (closure_type == 'en') then
         do icap=1,maxens_cap
            upmf_ave_cap3(icap) = sum(upmf_ens(11:13,:,:,icap))                            &
                                * 0.3333333333*maxens_efflsf_i
         end do
         do imbp=1,maxens_lsf
            upmf_ave_cap2(icap) = sum(upmf_ens(14:16,:,:,icap))                            &
                                * 0.3333333333*maxens_efflsf_i
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
      case ('nc','en')
         sgrell1_3d(6)  = .3333333 * sum(upmf_ave(1:3))
         sgrell1_3d(9)  = .3333333 * sum(upmf_ave(8:10))
         sgrell1_3d(10) = .25      * sum(upmf_ave(4:7))
         if (closure_type == 'en') then
            sgrell1_3d(7) = .3333333 * sum(upmf_ave(14:16))
            sgrell1_3d(8) = .3333333 * sum(upmf_ave(11:13))
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
      case ('nc','en')
         sgrell1_3d(11) = .3333333 * sum(upmf_ave(1:3))*hr_sec
         sgrell1_3d(14) = .3333333 * sum(upmf_ave(8:10))*hr_sec
         sgrell1_3d(15) = .25      * sum(upmf_ave(4:7))*hr_sec
         if (closure_type == 'en') then
            sgrell1_3d(12)= .3333333 * sum(upmf_ave(14:16))*hr_sec
            sgrell1_3d(13)= .3333333 * sum(upmf_ave(11:13))*hr_sec
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
      case ('nc')
         sgrell1_3d(21) = .3333333 * sum(upmf_ave(1:3))*hr_sec
         sgrell1_3d(22) = .25      * sum(upmf_ave(4:7))*hr_sec
         sgrell1_3d(23) = .3333333 * sum(upmf_ave(8:10))*hr_sec
      end select
   end if
   return
end subroutine grell_massflx_stats
!==========================================================================================!
!==========================================================================================!
