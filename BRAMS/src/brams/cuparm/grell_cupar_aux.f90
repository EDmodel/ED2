!==========================================================================================!
! grell_cupar_aux.f90                                                                      !
!                                                                                          !
!    This file contains subroutines that will handle the copying from BRAMS to Grell and   !
! vice-versa. This means that all initialization and output routines are here.             !
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This tiny subroutire initializes Grell's grid. This grid goes skips the bottom bound- !
! ary condition and in case of adaptive coordinate, the levels below surface               !
!------------------------------------------------------------------------------------------!
subroutine initialize_grid_grell(m1,deltax,deltay,zt,zm,flpw,rtgt,pblidx)
   use mem_scratch_grell, only: &
           mkx                  & ! intent(out) - Number of Grell levels.
          ,kgoff                & ! intent(out) - BRAMS offset related to Grell
          ,kpbl                 & ! intent(out) - Level of PBL top
          ,lpw                  & ! intent(out) - Lowest thermodynamic point
          ,tscal_kf             & ! intent(out) - Kain-Fritsch(1990) time scale
          ,z                    & ! intent(out) - Height
          ,z_cup                & ! intent(out) - Height at cloud levels
          ,dzu_cld              & ! intent(out) - Delta z for updraft calculations
          ,dzd_cld              ! ! intent(out) - Delta z for downdraft calculations
   implicit none
   !------ I/O variables ------------------------------------------------------------------!
   integer               , intent(in)  :: m1        ! Number of vertical levels
   real                  , intent(in)  :: deltax    ! Zonal resolution
   real                  , intent(in)  :: deltay    ! Meridional resolution
   real   , dimension(m1), intent(in)  :: zt        ! Height at thermodynamic levels
   real   , dimension(m1), intent(in)  :: zm        ! Height at momentum levels
   real                  , intent(in)  :: rtgt      ! Correction for terrain-following
   real                  , intent(in)  :: flpw      ! Level of lowest point above ground
   integer               , intent(in)  :: pblidx    ! Level of PBL top (Nakanishi/Niino)

   !------ Local variables ----------------------------------------------------------------!
   integer                             :: k

  !------ Initializing grid variables -----------------------------------------------------!
   lpw   = nint(flpw) ! This is always 2 for terrain-following, but varies for adaptive.
   kgoff = lpw-1      ! K-offset between BRAMS and Grell's grids.
   mkx   = m1-kgoff   ! Number of Grell levels.
   
   !------ Grid-dependent Kain-Fritsch time scale, based on Betts suggestion --------------!
   tscal_kf = 0.02 * sqrt(deltax*deltay)
   
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

   return
end subroutine initialize_grid_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine simply copies the tendencies to scratch arrays. It is done separated- !
! ly because it is the only place that we must give i and j information.                   !
!------------------------------------------------------------------------------------------!
subroutine initialize_tend_grell(m1,m2,m3,i,j,tht,rtt,ext,thsrc_shal,rtsrc_shal)
   use mem_scratch_grell, only : &
           mkx                   & ! intent(in)  - Number of Grell levels.
          ,kgoff                 & ! intent(in)  - BRAMS offset related to Grell
          ,lpw                   & ! intent(in)  - Lowest thermodynamic point
          ,dthetadt              & ! intent(out) - Temporary temperature tendency
          ,drtdt                 & ! intent(out) - Temporary mixing ratio tendency
          ,dthetadt_shal         & ! intent(out) - Temperature tend. due to shallower clds.
          ,drtdt_shal            & ! intent(out) - Temporary mix. ratio tend. due to shal.
          ,dexnerdt              ! ! Temporary Exner function tendency

   implicit none
   !------ I/O variables ------------------------------------------------------------------!
   integer, intent(in)                      :: m1,m2,m3   ! Grid dimensions
   integer, intent(in)                      :: i,j        ! Grid current position
   real   , intent(in), dimension(m1,m2,m3) :: tht        ! Potential temperature tend.
   real   , intent(in), dimension(m1,m2,m3) :: rtt        ! Total mixing ratio tend.
   real   , intent(in), dimension(m1,m2,m3) :: ext        ! Exner function tend.
   !------ This includes the tendencies due to convection that already happened -----------!
   real   , intent(in), dimension(m1,m2,m3) :: thsrc_shal ! Shallower conv. theta forcing
   real   , intent(in), dimension(m1,m2,m3) :: rtsrc_shal ! Shallower conv. mix.rat. forc.
   !------ Local variables ----------------------------------------------------------------!
   integer                                    :: k,kr     ! Counters

   do k=1,mkx
      kr=k+kgoff
      dthetadt(k)      = tht(kr,i,j)
      drtdt(k)         = rtt(kr,i,j)
      dexnerdt(k)      = ext(kr,i,j)
      dthetadt_shal(k) = thsrc_shal(kr,i,j)
      drtdt_shal(k)    = rtsrc_shal(kr,i,j)
   end do

   return
end subroutine initialize_tend_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine initializes Grell's thermodynamic and turbulence fields. If it        !
! is a deep convection call, then it will include the effect of shallow convection, and    !
! make the variables for the case in which no convection happens consistent with this.     !
!------------------------------------------------------------------------------------------!
subroutine initialize_thermo_grell(m1,dtime,theta,rv,pi0,pp,wa,dn0,tke,rcp)

   use mem_scratch_grell, only : &
          dthetadt               & ! intent(in)  - Temporary temperature tendency
         ,dthetadt_shal          & ! intent(in)  - Temperature tendency due to shallower
         ,drtdt                  & ! intent(in)  - Temporary mixing ratio tendency
         ,drtdt_shal             & ! intent(in)  - Mixing ratio tendency due to shallower
         ,dexnerdt               & ! intent(in)  - Temporary Exner function tendency
         ,mkx                    & ! intent(in)  - Number of Grell levels.
         ,kgoff                  & ! intent(in)  - BRAMS offset related to Grell
         ,lpw                    & ! intent(in)  - Lowest thermodynamic point
         ,z                      & ! intent(in)  - Grell's heights
         ,he0                    & ! intent(out) - Current Moist static energy (MSE)
         ,he                     & ! intent(out) - Forced  MSE
         ,hesur                  & ! intent(out) - Sfc. Moist static energy (MSE)
         ,hes                    & ! intent(out) - Forced Saturation MSE
         ,hes0                   & ! intent(out) - Current Saturation MSE
         ,hessur                 & ! intent(out) - Sfc. Saturation Moist static energy (SMSE)
         ,p0                     & ! intent(out) - Current pressure
         ,p                      & ! intent(out) - Forced Pressure
         ,q0                     & ! intent(out) - Current Specific humidity.
         ,q                      & ! intent(out) - Forced Water vapour mixing ratio
         ,qsur                   & ! intent(out) - Sfc. Water vapour mixing ratio
         ,qes0                   & ! intent(out) - Current Sat. mixing ratio 
         ,qes                    & ! intent(out) - Forced Sat. water vapour mixing ratio
         ,qessur                 & ! intent(out) - Sfc. Saturation water vapour mixing ratio
         ,psur                   & ! intent(out) - Surface pressure
         ,t0                     & ! intent(out) - Current Temperature
         ,t                      & ! intent(out) - Forced Temperature
         ,tsur                   & ! intent(out) - Surface temperature
         ,omeg                   & ! intent(out) - Lagrangian pressure tendency 
         ,mconv                  & ! intent(out) - Moisture convergence
         ,tkeg                   & ! intent(out) - Turbulent kinetic energy
         ,rcpg                   & ! intent(out) - Cloud droplet mixing ratio
         ,rho                    & ! intent(out) - Density
         ,wwind                  ! ! intent(out) - Vertical velocity

   use rconstants, only : cp,cpi,cpor,p00,g,alvl,rgas

   implicit none
   !------ I/O variables ------------------------------------------------------------------!
   integer, intent(in)                  :: m1
   real   , intent(in)                  :: dtime
   real   , intent(in)  , dimension(m1) :: theta,rv,pi0,pp,dn0,wa,tke,rcp

   !------ Local variables ----------------------------------------------------------------!
   integer                              :: k,kr
   real                                 :: exner,theta_shal,theta_forced
   real                                 :: dq,next_pres,next_exner
   real, external                       :: rslf,virtt

   do k=1,mkx
      kr=k+kgoff
      
      !------ Finding current and future Exner functions ----------------------------------!
      exner      = pi0(kr)+pp(kr)
      next_exner = exner + dexnerdt(k)*dtime

      
      !------ Finding current and future pressure. Both in Pa. Now Grell uses Pa, not hPa -!
      p0(k)   = p00*(cpi*exner)**cpor
      p(k)    = p00*(cpi*next_exner)**cpor


      !------ Temperature with shallow cumulus effect (if deep convection call) -----------!
      theta_shal   = theta(kr)+dtime*+dthetadt_shal(k)
      t0(k)        = cpi*theta_shal*exner
      
      
      !------ Temperature with the tendency term ------------------------------------------!
      ! (If shallow convection happened, dthetadt already contains this forcing)           !
      theta_forced = theta(kr)+dthetadt(k)*dtime
      t(k)         = cpi*theta_forced*next_exner
      if (t0(k) < 200.) t(k) = t0(k)
      
      
      !------ Mixing ratios with shallow cumulus effect (if deep convection call) ---------!
      qes0(k)      = max(1.e-8,rslf(p0(k),t0(k)))  ! This is the saturation mixing ratio
      q0(k)        = min(max(1.e-8,rv(kr)+drtdt_shal(k)*dtime),qes0(k))


      !------ Expected mixing ratios at the next time step --------------------------------!
      qes(k)     = max(1.e-8,rslf(p(k),t(k)))
      q(k)       = min(max(1.e-8, rv(kr) +  drtdt(k) * dtime),qes(k))


      !------ Finding moist static energies -----------------------------------------------!
      hes0(k)      = g * z(k) + cp * t0(k) + alvl * qes0(k)
      he0(k)       = g * z(k) + cp * t0(k) + alvl * q0(k)


      !------ Finding moist static energies with the tendency term added ------------------!
      hes(k)      = g * z(k) + cp * t(k) + alvl * qes(k)
      he(k)       = g * z(k) + cp * t(k) + alvl * q(k)


      !------ Vertical velocity in terms of pressure, or Lagrangian dp/dt [ Pa/s] ---------!
      omeg(k)     = -g*dn0(kr)*.5*( wa(kr)+wa(k) )


      !------ Turbulent kinetic energy [m²/s²] --------------------------------------------!
      tkeg(k)     = tke(kr)
      
      
      !------ Cloud droplet mixing ratio [kg/kg] ------------------------------------------!
      rcpg(k)     = rcp(kr)
      
      !------ Vertical velocity [m/s], this is staggered, averaging... --------------------!
      wwind(k)    = 0.5 * (wa(kr)+wa(kr-1))

      !------ Air density [kg/m³] ---------------------------------------------------------!
      rho(k) = p(k) / (rgas * virtt(t(k),q(k)))

   end do
   
   !----- Computing the surface variables -------------------------------------------------!
   psur        = .5*p00*( (cpi * (pp(lpw-1)+pi0(lpw-1)))**cpor +  &
                          (cpi * (pp(lpw)  +pi0(lpw)  ))**cpor )
   tsur        = cpi*theta(lpw)*(pi0(lpw)+pp(lpw)) 
   qessur      = max(1.e-8,rslf(psur,tsur))
   qsur        = rv(lpw)
   hessur      = cp * tsur + alvl * qessur
   hesur       = cp * tsur + alvl * qsur
   
   !------ Finding the integrated moisture convergence ------------------------------------!
   mconv=0.
   do k=2,mkx-1
      dq          = .5*(q(k+1)-q(k-1))  ! Delta-q between layers.
      mconv       = mconv+omeg(k)*dq/g  ! Moisture convergence, or (-1)*moisture divergence.
   end do
   mconv     = max(0.,mconv)            ! Throwing away negative convergence (aka divergence).

   return
end subroutine initialize_thermo_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!   This subroutine intializes the upstream wind information, that may affect convection   !
! at the current grid cell. This is done only if the user asked for it.                    !
!------------------------------------------------------------------------------------------!
subroutine initialize_upstream_grell(comp_down,iupstrm,m1,m2,m3,i,j,jdim,last_dnmf         &
                                    ,last_ierr,ua,va)

   use mem_scratch_grell, only: &
            mkx                 & ! intent(in)  - Number of Grell levels.
           ,kgoff               & ! intent(in)  - BRAMS offset related to Grell
           ,lpw                 & ! intent(in)  - Lowest thermodynamic point
           ,p                   & ! intent(in)  - Pressure at Grell's levels
           ,psur                & ! intent(in)  - Pressure at surface
           ,dens_curr           & ! intent(out) - Downdraft mass flux as density current
           ,prev_dnmf           & ! intent(out) - Downdraft mass flux at the previous call
           ,uwind               & ! intent(out) - Zonal wind at thermodynamic point
           ,vwind               & ! intent(out) - Meridional wind at thermodynamic point
           ,upconv              & ! intent(out) - Flag for upstream convection.
           ,z                   ! ! intent(in)  - Height at Grell's levels
   
   use rconstants       , only: cpi,cpor,p00,g,onerad

   use grell_coms       , only: pbotmean, ptopmean,vspeedmin
   
   implicit none
   !------ I/O variables ------------------------------------------------------------------!
   integer                      , intent(in) :: iupstrm   ! Consider upstream downdrafts
   integer                      , intent(in) :: m1,m2,m3  ! Number of z,x, and y points
   integer                      , intent(in) :: i,j       ! Current x and y position
   integer                      , intent(in) :: jdim      ! Dimension in y
   logical                      , intent(in) :: comp_down ! Computing downdrafts (T/F)
   real   , dimension   (m2,m3) , intent(in) :: last_dnmf ! Last time downdraft
   real   , dimension   (m2,m3) , intent(in) :: last_ierr ! Last time error flag
   real   , dimension(m1,m2,m3) , intent(in) :: ua        ! Zonal wind
   real   , dimension(m1,m2,m3) , intent(in) :: va        ! Meridional wind

   !------ Local variables ----------------------------------------------------------------!
   integer             :: k         ! Counter for current Grell level
   integer             :: kr        ! Counter for corresponding BRAMS level
   integer             :: iwest     ! Neighbour point to the west
   integer             :: ieast     ! Neighbour point to the east
   integer             :: jsouth    ! Neighbour point to the south
   integer             :: jnorth    ! Neighbour point to the north
   integer             :: myoctant  ! Octant of upstream flow 
   real                :: dp        ! Layer thickness in terms of pressure          [   Pa]
   real                :: dz        ! Layer thickness in terms of height            [    m]
   real                :: pthick    ! Total thickness in terms of pressure          [   Pa]
   real                :: zthick    ! Total thickness in terms of height            [    m]
   real                :: vspeed    ! Average wind speed                            [  m/s]
   real                :: direction ! Average wind direction                        [  deg]
   real                :: umean     ! Average zonal wind                            [  m/s]
   real                :: vmean     ! Average meridional wind                       [  m/s]
   
   !------ Some local parameters. ---------------------------------------------------------!
   integer, parameter :: nocta    =  8  ! Number of octants
   real,    parameter :: sizeocta = 45. ! Size of each octant (360/8=45.)
   !---------------------------------------------------------------------------------------!

   
   !------ Initializing scalars -----------------------------------------------------------!
   prev_dnmf = last_dnmf(i,j)
   dens_curr = 0.
   upconv    = .false.
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

   if ((.not. comp_down) .or. iupstrm == 0) return

   !------ Initializing neighbour indices -------------------------------------------------!
   iwest  = i-1
   ieast  = i+1
   jsouth = j-1
   jnorth = j+1

   !---------------------------------------------------------------------------------------!
   !   Computing some column means. Here I switched the pressure average by height average,!
   ! it is more intuitive in my opinion. I left the tools to switch back to pressure       !
   ! average, though.                                                                      !
   !---------------------------------------------------------------------------------------!
   !------ Initializing column means and scales -------------------------------------------!
   umean       = 0.
   vmean       = 0.
   pthick      = 0.
   zthick      = 0.
   
   !------ Looping through layers, and integrating column. --------------------------------!
   do k=2,mkx-1

      !------ Wind integral done just at the mid atmosphere -------------------------------!
      if((psur-p(k)) > pbotmean .and. p(k) < ptopmean) then
         dp       = -.5*(p(k+1)-p(k-1))             ! Old method
         dz       =  .5*(z(k+1)-z(k-1))
         umean    = umean+uwind(k)*dz
         vmean    = vmean+vwind(k)*dz
         pthick   = pthick+dp
         zthick   = zthick+dz
      end if
   end do
   
   !------ Normalizing variables ----------------------------------------------------------!
   umean     = umean/zthick  ! Scaling the integral so we have mean zonal wind.
   vmean     = vmean/zthick  ! Scaling the integral so we have mean meridional wind.
   vspeed    = sqrt(umean*umean+vmean*vmean)
   !------ Compute wind direction only if wind is enough strong ---------------------------!
   if(vspeed < vspeedmin) then
      ! Wind is too weak to define upstream, so the "upstream"  is here... ----------------!
      myoctant=-999
      upconv = last_ierr(i,j) == 0.
      if (iupstrm > 0) then
         dens_curr = prev_dnmf
      end if
      return
   else
      !----- Defining direction of the mean wind, and the corresponding quadrant ----------!
      direction = modulo(270.-atan2(vmean,umean)*onerad,360.)
      myoctant  = mod(nint(direction/sizeocta),nocta)
   end if
   !---------------------------------------------------------------------------------------!


   !----- Now check for upstream convection depending on the octant -----------------------!
   select case (myoctant)
   case (0) ! Northerly
     upconv = last_ierr(i,jnorth) == 0.
     if (iupstrm == 2) prev_dnmf=max(last_dnmf(i,j),last_dnmf(i,jnorth))

   case (1) ! Northeasterly flow
      upconv = last_ierr(ieast,jnorth) == 0.
      if (iupstrm == 2) prev_dnmf=max(last_dnmf(i,j),last_dnmf(ieast,jnorth))
      
   case (2) ! Easterly flow
      upconv = last_ierr(ieast,j) == 0.
      if (iupstrm == 2) prev_dnmf=max(last_dnmf(i,j),last_dnmf(ieast,j))
      
   case (3) ! Southeasterly flow
      upconv = last_ierr(ieast,jsouth) == 0.
      if (iupstrm == 2) prev_dnmf=max(last_dnmf(i,j),last_dnmf(ieast,jsouth))
      
   case (4) ! Southerly flow
      upconv = last_ierr(i,jsouth) == 0.
      if (iupstrm == 2) prev_dnmf=max(last_dnmf(i,j),last_dnmf(i,jsouth))
      
   case (5) ! Southwesterly flow
      upconv = last_ierr(iwest,jsouth) == 0.
      if (iupstrm == 2) prev_dnmf=max(last_dnmf(i,j),last_dnmf(iwest,jsouth))
      
   case (6) ! Westerly flow
      upconv = last_ierr(iwest,j) == 0.
      if (iupstrm == 2) prev_dnmf=max(last_dnmf(i,j),last_dnmf(iwest,j))
      
   case (7) ! Northwesterly flow
      upconv = last_ierr(iwest,jnorth) == 0.
      if (iupstrm == 2) prev_dnmf=max(last_dnmf(i,j),last_dnmf(iwest,jnorth))
      
   end select
   
   !---------------------------------------------------------------------------------------!
   !    Dens_curr is simulating the impact of a previous downdraft lifting the air         !
   ! downstream like a density current. It will only do it if iupstrm is > 0               !
   !---------------------------------------------------------------------------------------!
   if (iupstrm > 0) then
      dens_curr = prev_dnmf
   end if

   return
end subroutine initialize_upstream_grell
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
!    This subroutine organizes the variables that go to the output, namely the heating and !
! moistening rates due to convection.                                                      !
!------------------------------------------------------------------------------------------!
subroutine grell_output(comp_down,m1,mgmzp,rtgt,zt,zm,pi0,pp,dnmf,upmf,xierr,zjmin,zk22    &
                       ,zkbcon,zkdt,zktop,conprr,thsrc,rtsrc,areadn,areaup,cupcond)
   use rconstants, only: cp
   use mem_scratch_grell, only: &
           cdd                  & ! intent(in) - Normalized downdraft detr. rate  [    1/m]
          ,cdu                  & ! intent(in) - Normalized updraft detr. rate    [    1/m]
          ,dzd_cld              & ! intent(in) - Top-down layer thickness         [      m]
          ,dzu_cld              & ! intent(in) - Bottom-up layer thickness        [      m]
          ,etad_cld             & ! intent(in) - Normalized downdraft mass flux   [    ---]
          ,etau_cld             & ! intent(in) - Normalized updraft mass flux     [    ---]
          ,ierr                 & ! intent(in) - Error flag
          ,jmin                 & ! intent(in) - Downdraft originating level
          ,k22                  & ! intent(in) - Updraft originating level
          ,kbcon                & ! intent(in) - Cloud base
          ,kdet                 & ! intent(in) - Top of downdraft detrainment layer
          ,kgoff                & ! intent(in) - BRAMS grid offset
          ,ktop                 & ! intent(in) - Cloud top
          ,mentrd_rate          & ! intent(in) - Normalized downdraft entr. rate  [    1/m]
          ,mentru_rate          & ! intent(in) - Normalized updraft entr. rate    [    1/m]
          ,mkx                  & ! intent(in) - # of cloud grid levels
          ,outq                 & ! intent(in) - Vapour mixing ratio forcing      [kg/kg/s]
          ,outqc                & ! intent(in) - Liq. water mixing ratio forcing  [kg/kg/s]
          ,outt                 & ! intent(in) - Temperature forcing              [    K/s]
          ,p_cup                & ! intent(in) - Pressure @ Grell's levels        [     Pa]
          ,q_cup                & ! intent(in) - Mixing ratio @ Grell's levels    [  kg/kg]
          ,qd_cld               & ! intent(in) - Downdraft vapour mixing ratio    [  kg/kg]
          ,qu_cld               & ! intent(in) - Updraft vapour mixing ratio      [  kg/kg]
          ,qld_cld              & ! intent(in) - Downdraft liquid mixing ratio    [  kg/kg]
          ,qlu_cld              & ! intent(in) - Updraft liquid mixing ratio      [  kg/kg]
          ,precip               & ! intent(in) - Precipitation rate               [kg/m²/s]
          ,rho                  & ! intent(in) - Air density                      [  kg/m³]
          ,tkeg                 & ! intent(in) - TKE at the model levels          [   J/kg]
          ,t_cup                & ! intent(in) - Temperature @ Grell's levels     [      K]
          ,td_cld               & ! intent(in) - Downdraft temperature            [      K]
          ,tu_cld               & ! intent(in) - Updraft temperature              [      K]
          ,wwind                ! ! intent(in) - Vertical velocity                [    m/s]

   implicit none
   
   logical            , intent(in)  :: comp_down ! I will compute downdraft stuff
   integer            , intent(in)  :: m1        ! Number of levels
   integer            , intent(in)  :: mgmzp     ! Number of Grell's levels
   real               , intent(in)  :: rtgt      ! Correction to get the heights  [   ----]
   real, dimension(m1), intent(in)  :: zt        ! Height at thermodynamic levels [      m]
   real, dimension(m1), intent(in)  :: zm        ! Height at momentum levels      [      m]
   real, dimension(m1), intent(in)  :: pi0       ! Exner function reference       [ J/kg/K]
   real, dimension(m1), intent(in)  :: pp        ! Exner function perturbation    [ J/kg/K]
   real               , intent(in)  :: dnmf      ! Reference downdraft mass flux  [kg/m²/s]
   real               , intent(in)  :: upmf      ! Reference updraft mass flux    [kg/m²/s]
   
   real, dimension(m1), intent(out) :: thsrc     ! Potential temperature feedback [    K/s]
   real, dimension(m1), intent(out) :: rtsrc     ! Total mixing ratio feedback    [kg/kg/s]
   real, dimension(m1), intent(out) :: cupcond   ! Updraft condensed water        [  kg/kg]
   real, dimension(m1), intent(out) :: areadn    ! Fractional downdraft area      [    ---]
   real, dimension(m1), intent(out) :: areaup    ! Fractional updraft area        [    ---]
   real               , intent(out) :: conprr    ! Rate of convective precip.     [kg/m²/s]
   real               , intent(out) :: xierr     ! Error flag
   real               , intent(out) :: zjmin     ! Downdraft originating level    [      m]
   real               , intent(out) :: zk22      ! Cloud base                     [      m]
   real               , intent(out) :: zkbcon    ! Level of free c                [      m]
   real               , intent(out) :: zkdt      ! Top of the dndraft detr. layer [      m]
   real               , intent(out) :: zktop     ! Cloud top                      [      m]
   
   integer                          :: k         ! Cloud level counter
   integer                          :: kr        ! BRAMS level counter
   real                             :: exner     ! Exner fctn. for tend. convers. [ J/kg/K]

   !---------------------------------------------------------------------------------------!
   !    Flushing all variables to zero in case convection didn't happen.                   !
   !---------------------------------------------------------------------------------------!
   thsrc   = 0.
   rtsrc   = 0.
   areadn  = 0.
   areaup  = 0.
   conprr  = 0.
   cupcond = 0.
   zkdt    = 0.
   zk22    = 0.
   zkbcon  = 0.
   zjmin   = 0.
   zktop   = 0.
   
   !---------------------------------------------------------------------------------------!
   !   Copying the error flag. This should not be zero in case convection failed.          !
   !---------------------------------------------------------------------------------------!
   xierr = real(ierr)

   !---------------------------------------------------------------------------------------!
   !    Fixing the levels, here I will add back the offset so the output will be consist-  !
   ! ent. I will return these variables even when no cloud developed for debugging         !
   ! purposes. When the code is running fine, then I should return them only when          !
   ! convection happens.                                                                   !
   !---------------------------------------------------------------------------------------!
   zkdt   = (zt(kdet  + kgoff)-zm(kgoff))*rtgt
   zk22   = (zt(k22   + kgoff)-zm(kgoff))*rtgt
   zkbcon = (zt(kbcon + kgoff)-zm(kgoff))*rtgt
   zjmin  = (zt(jmin  + kgoff)-zm(kgoff))*rtgt
   zktop  = (zt(ktop  + kgoff)-zm(kgoff))*rtgt


   if (ierr /= 0) return !----- No cloud, no need to return anything ----------------------!
   !---------------------------------------------------------------------------------------!
   !    Precipitation is simply copied, it could even be output directly from the main     !
   ! subroutine, brought here just to be together with the other source terms.             !
   !---------------------------------------------------------------------------------------!
   conprr = precip
   
   do k=1,ktop
      kr    = k + kgoff
      exner = pi0(kr)+pp(kr)
      !------------------------------------------------------------------------------------!
      !    I think this is okay as long as pressure is not directly changed by the cumulus !
      ! parameterization.                                                                  !
      !------------------------------------------------------------------------------------!
      thsrc(kr)   = cp * outt(k) / exner
      rtsrc(kr)   = outq(k) + outqc(k)
   end do


   !---------------------------------------------------------------------------------------!
   !   Computing the relative area covered by downdrafts and updrafts.                     !
   !---------------------------------------------------------------------------------------!
   call grell_draft_area(comp_down,m1,mgmzp,kgoff,jmin,k22,kbcon,ktop,rho,wwind,tkeg,p_cup &
                        ,q_cup,t_cup,dzd_cld,etad_cld,mentrd_rate,cdd,qd_cld,td_cld        &
                        ,dzu_cld,etau_cld,mentru_rate,cdu,qu_cld,tu_cld,dnmf,upmf,areadn   &
                        ,areaup)


   !---------------------------------------------------------------------------------------!
   !   Finally I compute the cloud condensed mixing ratio. This is used by Harrington when !
   ! cumulus feedback is requested, so I rescale the liquid water at the downdrafts and    !
   ! updrafts by their area (roughly stretching the cloud to the entire grid cell and find-!
   ! ing the "equivalent stratus".                                                         !
   !---------------------------------------------------------------------------------------!
   do k=1,ktop
      kr=k+kgoff
      cupcond(kr) = max(0.,qld_cld(k) * areadn(kr) + qlu_cld(k) * areaup(kr))
   end do

   return
end subroutine grell_output
!==========================================================================================!
!==========================================================================================!






!==========================================================================================!
!==========================================================================================!
subroutine grell_draft_area(comp_down,m1,mgmzp,kgoff,jmin,k22,kbcon,ktop,rho,wwind,tkeg    &
                           ,p_cup,q_cup,t_cup,dzd_cld,etad_cld,mentrd_rate,cdd,qd_cld      &
                           ,td_cld,dzu_cld,etau_cld,mentru_rate,cdu,qu_cld,tu_cld,dnmf     &
                           ,upmf,areadn,areaup)
!------------------------------------------------------------------------------------------!
!    This subroutine computes the relative downdraft and updraft areas. This is used by    !
! some Lagrangian models and it will be also used to feedback the liquid water content.    !
!                                                                                          !
!  The references for this code are:                                                       !
!                                                                                          !
! Anthes, R.A.; A cumulus parameterization scheme utilizing a one-dimensional cloud model. !
!      Mon. Wea. Rev., vol. 105, 270-286 (an or AN77).                                     !
!                                                                                          !
! Fritsch, J.M; Chappell, C. F., 1980: Numerical prediction of convectively driven         !
!      mesoscale pressure systems. Part I: Convective parameterization. J. Atmos. Sci.,    !
!      vol. 37(8), 1722-1733. (fc or FC80).                                                !
!                                                                                          !
! Zhang, D.-L.; Fritsch, J.M., 1986: Numerical simulation of the meso-beta scale structure !
!      and evolution of the 1977 Johnstown flood. Part I: Model description and            !
!      verification. J. Atm. Sci., vol. 43(18). 1913-1943. (zf or ZF86).                   !
!                                                                                          !
!    One important difference from AN77/ZF86 is that here we considered not only the       !
! entrainment, but also the detrainment to estimate the draft velocities.                  !
!------------------------------------------------------------------------------------------!

   use rconstants, only : cp, g, ep, rgas

   implicit none
   
   !----- Input variables, flags ----------------------------------------------------------!
   logical               , intent(in)  :: comp_down   ! Flag for downdraft computation
   !----- Input variables, grid boundaries ------------------------------------------------!
   integer               , intent(in)  :: m1          ! Number of BRAMS levels
   integer               , intent(in)  :: mgmzp       ! Number of scratch levels
   integer               , intent(in)  :: kgoff       ! BRAMS offset
   !----- Input variables, cloud boundaries -----------------------------------------------!
   integer               , intent(in)  :: jmin        ! Downdraft origin
   integer               , intent(in)  :: k22         ! Updraft origin
   integer               , intent(in)  :: kbcon       ! Level of free convection
   integer               , intent(in)  :: ktop        ! Cloud top
   !----- Input variables, thermodynamic and dynamic properties ---------------------------!
   real, dimension(mgmzp), intent(in)  :: rho         ! Air density               [  kg/m³]
   real, dimension(mgmzp), intent(in)  :: wwind       ! Vertical velocity         [    m/s]
   real, dimension(mgmzp), intent(in)  :: tkeg        ! Turbulent Kinetic Energy  [   J/kg]
   real, dimension(mgmzp), intent(in)  :: p_cup       ! Environment pressure      [     Pa]
   real, dimension(mgmzp), intent(in)  :: q_cup       ! Environment mixing ratio  [  kg/kg]
   real, dimension(mgmzp), intent(in)  :: t_cup       ! Environment temperature   [      K]
   real, dimension(mgmzp), intent(in)  :: dzd_cld     ! Top-down layer thickness  [      m]
   real, dimension(mgmzp), intent(in)  :: etad_cld    ! Normalized downdraft flux [    ---]
   real, dimension(mgmzp), intent(in)  :: mentrd_rate ! Normalized downdraft entr.[    ---]
   real, dimension(mgmzp), intent(in)  :: cdd         ! Normalized downdraft detr.[    ---]
   real, dimension(mgmzp), intent(in)  :: qd_cld      ! Downdraft mixing ratio    [  kg/kg]
   real, dimension(mgmzp), intent(in)  :: td_cld      ! Downdraft temperature     [      K]
   real, dimension(mgmzp), intent(in)  :: dzu_cld     ! Bottom-up layer thickness [      m]
   real, dimension(mgmzp), intent(in)  :: etau_cld    ! Normalized updraft flux   [    ---]
   real, dimension(mgmzp), intent(in)  :: mentru_rate ! Normalized updraft entr.  [    ---]
   real, dimension(mgmzp), intent(in)  :: cdu         ! Normalized updraft detr.  [    ---]
   real, dimension(mgmzp), intent(in)  :: qu_cld      ! Updraft mixing ratio      [  kg/kg]
   real, dimension(mgmzp), intent(in)  :: tu_cld      ! Updraft temperature       [      K]
   real                  , intent(in)  :: dnmf        ! Downdraft ref. mass flux  [kg/m²/s]
   real                  , intent(in)  :: upmf        ! Updraft ref. mass flux    [kg/m²/s]
   !----- Output variables ----------------------------------------------------------------!
   real, dimension(m1)   , intent(out) :: areadn      ! Downdraft relative area   [    ---]
   real, dimension(m1)   , intent(out) :: areaup      ! Updraft   relative area   [    ---]
   !----- Local variables -----------------------------------------------------------------!
   integer                             :: k         ! Cloud level counter
   integer                             :: kr        ! BRAMS level counter
   real                                :: buoyancy  ! Buoyancy acceleration       [   m/s²]
   real                                :: rho_cld   ! Scratch draft density       [  kg/m³]
   real                                :: tv        ! Scratch virtual temperature [      K]
   real                                :: tv_cld    ! Scratch draft virtual temp. [      K]
   real                                :: w_cld     ! Scratch draft vert. veloc.  [    m/s]
   real                                :: mudz_cld  ! Scratch draft µe-µd         [    ---]
   !----- Constant from AN77, to mitigate the neglection of non-hydrostatic effect. ------!
   real, parameter                     :: virt_mass_an77 = 0.5    ! "Virtual mass"
   !----- Minimum draft to prevent zero downdraft -----------------------------------------!
   real, parameter                     :: min_downdraft      = 0.10
   real, parameter                     :: min_updraft        = 0.60
   !----- Function to compute virtual temperature -----------------------------------------!
   real, external                      :: virtt
   !---------------------------------------------------------------------------------------!

   
   !---------------------------------------------------------------------------------------!
   !    Setting the area to zero. This will be replaced by the actual area at the draft    !
   ! layer.                                                                                !
   !---------------------------------------------------------------------------------------!
   areadn=0.
   areaup=0.
   
   !---------------------------------------------------------------------------------------!
   !    Downdraft. As in Fritsch and Chappell (1980), I assume downdraft velocity to be    !
   ! zero at the bottom, and integrate backwards to get the profile.                       !
   !---------------------------------------------------------------------------------------!
   if (comp_down) then
      w_cld = 0.
      do k=2,jmin
         kr=k+kgoff
         tv_cld   = 0.5* (virtt(td_cld(k),qd_cld(k)) + virtt(td_cld(k-1),qd_cld(k-1)))
         tv       = 0.5* (virtt(t_cup(k),q_cup(k))   + virtt(t_cup(k-1),q_cup(k-1)))
         buoyancy = g  * (tv_cld-tv) / ((1 + virt_mass_an77)*tv)
         mudz_cld = 0.5* dzu_cld(k)* (sqrt(cdd(k-1)**2+mentrd_rate(k-1)**2)                &
                                     +sqrt(cdd(k)**2+mentrd_rate(k)**2))
         w_cld    = sqrt(max((w_cld*rho(k)/rho(k-1))**2                                    &
                           ,(w_cld*w_cld*(1-mudz_cld)-2.*buoyancy*dzu_cld(k))/(1+mudz_cld)))
         if (w_cld == 0) w_cld = min_downdraft
         rho_cld    = sqrt(p_cup(k)*p_cup(k-1))  / (rgas * tv_cld)
         areadn(kr) = dnmf * etad_cld(k) / (rho_cld * w_cld)
      end do
   end if
   
   !---------------------------------------------------------------------------------------!
   !   Updraft. Here I am assuming that the velocity at the updraft origin depends on the  !
   ! TKE. Many parameterizations define this velocity as w0+[2.*TKE(k22)^½], where w0 is   !
   ! often assumed 1. I will try to use the large scale velocity as my w0 and see if the   !
   ! area does not get too large. If it gets, then I will increase w0.                     !
   ! With this information, I go all the way to the top and go back to the updraft origin. !
   !---------------------------------------------------------------------------------------!
   !----- Updraft at the updraft origin ---------------------------------------------------!
   w_cld   = max(min_updraft,0.5*(wwind(k22)-wwind(k22-1)) + sqrt(tkeg(k22)+tkeg(k22-1)))
   tv_cld  = 0.5*(virtt(tu_cld(kbcon-1),qu_cld(kbcon-1))+virtt(tu_cld(kbcon),qu_cld(kbcon)))
   rho_cld = sqrt(p_cup(kbcon)*p_cup(kbcon-1)) / (rgas*tv_cld)
   areaup(k22+kgoff) = upmf * etau_cld(kbcon) / (rho_cld * w_cld)
   !----- Updraft for other levels --------------------------------------------------------!
   do k=k22+1,ktop
      kr=k+kgoff
      tv_cld   = 0.5* (virtt(tu_cld(k-1),qu_cld(k-1)) + virtt(tu_cld(k),qu_cld(k)))
      tv       = 0.5* (virtt(t_cup(k-1),q_cup(k-1))   + virtt(t_cup(k),q_cup(k)))
      buoyancy = g  * (tv_cld-tv) / ((1 + virt_mass_an77)*tv)
      mudz_cld = 0.5* dzu_cld(k)* (sqrt(cdu(k-1)**2+mentru_rate(k-1)**2)                   &
                                  +sqrt(cdu(k)**2+mentru_rate(k)**2))
      w_cld     = sqrt(max((w_cld *rho(k)/rho(k-1))**2                                     &
                         ,(w_cld * w_cld*(1-mudz_cld)+ 2*buoyancy*dzu_cld(k))/(1+mudz_cld)))
      rho_cld   = sqrt(p_cup(k)*p_cup(k-1)) / (rgas * tv_cld)
      areaup(kr)= upmf * etau_cld(k) / (rho_cld * w_cld)
   end do


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
                              ,maxens_cap,inv_ensdim,closure_type,upmf_ens,sgrell1_3d      &
                              ,sgrell2_3d)

   use rconstants, only : hr_sec
   use mem_scratch_grell, only: &
           ierr                 & ! intent(in) - Error flag
          ,jmin                 & ! intent(in) - Downdraft originating level
          ,k22                  & ! intent(in) - Updraft originating level
          ,kbcon                & ! intent(in) - Cloud base
          ,kdet                 & ! intent(in) - Top of downdraft detrainment layer
          ,kgoff                & ! intent(in) - BRAMS grid offset
          ,ktop                 & ! intent(in) - Cloud top
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
   
   real, dimension(maxens_dyn,maxens_lsf,maxens_eff,maxens_cap), intent(in) &
                                      :: upmf_ens ! Ensemble

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


   if (ierr /= 0) then
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
         upmf_ave_cap4(icap) = sum(upmf_ens(8:10,:,:,icap))*0.3333333333*maxens_efflsf_i
      end do
      if (closure_type == 'en') then
         do icap=1,maxens_cap
            upmf_ave_cap3(icap) = sum(upmf_ens(11:13,:,:,icap))*0.3333333333*maxens_efflsf_i
         end do
         do imbp=1,maxens_lsf
            upmf_ave_cap2(icap) = sum(upmf_ens(14:16,:,:,icap))*0.3333333333*maxens_efflsf_i
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
