!
!   ##########################################################################
subroutine URBAN(PTS_TOWN, PEMIS_TOWN, PALB_TOWN,             &
     PT_CANYON, PQ_CANYON, PU_CANYON,                         &
     PTS_ROOF,PTS_ROAD,PTS_WALL,PTI_ROAD,PTI_BLD,             &
     PT_ROOF, PT_ROAD, PT_WALL, PWS_ROOF,PWS_ROAD,            &
     PPS, PPA, PEXNS, PEXNA, PTA, PQA, PRHOA,                 &
     PLW_RAD, PDIR_SW_RAD, PSCA_SW_RAD, PTANZEN,              &
     PRR,                                                     &
     PZREF, PDIRCOSZW, PUSLOPE, PVSLOPE, PVMOD,               &
     PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
     PTSTEP,                                                  &
     PZ0_TOWN,                                                &
     PBLD,PBLD_HEIGHT,PWALL_O_HOR,PCAN_HW_RATIO,              &
     PALB_ROOF, PEMIS_ROOF,                                   &
     PHC_ROOF,PTC_ROOF,PD_ROOF,                               &
     PALB_ROAD, PEMIS_ROAD, PSVF_ROAD,                        &
     PHC_ROAD,PTC_ROAD,PD_ROAD,                               &
     PALB_WALL, PEMIS_WALL, PSVF_WALL,                        &
     PHC_WALL,PTC_WALL,PD_WALL,                               &
     PRN_ROOF, PH_ROOF, PLE_ROOF, PGFLUX_ROOF,                &
     PRUNOFF_ROOF,               		              &
     PRN_ROAD, PH_ROAD, PLE_ROAD, PGFLUX_ROAD,                &
     PRUNOFF_ROAD,                                            &
     PRN_WALL, PH_WALL, PLE_WALL, PGFLUX_WALL,                &
     PRN_TOWN, PH_TOWN, PLE_TOWN, PGFLUX_TOWN,                &
     PRUNOFF_TOWN, PSFU_TOWN, PSFV_TOWN, PCH_TOWN             )
  !   ##########################################################################
  !
  !!****  *URBAN*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the evolution of prognostic variables and the fluxes
  !     over artificial surfaces as towns, taking into account the canyon like
  !     geometry of urbanized areas.
  !         
  !     
  !!**  METHOD
  !     ------
  !
  !     The prognostic variables are:
  !       - the surface temperature for roofs, roads, and walls
  !       - the water reservoir, whose maximum value is 10mm
  !
  !
  !    1 : computation of input solar radiation on each surface
  !        ****************************************************
  !
  !
  !
  !    2 : drag coefficient for momentum 
  !        *****************************
  !
  !
  !    3 : aerodynamical resistance for heat transfers
  !        *******************************************
  !
  !
  !    4 : equation for evolution of Ts_roof
  !        *********************************
  !
  !
  !       Rn = (dir_Rg + sca_Rg) (1-a) + emis * ( Rat - sigma Ts**4 (t+dt) )
  !
  !       H  = rho Cp CH V ( Ts (t+dt) - Tas )
  !
  !       LE = rho Lv CH V ( qs (t+dt) - qas )
  !
  !      where the as subscript denotes atmospheric values at ground level
  !      (and not at first half level)
  !
  !
  !    5 : equations for evolution of Ts_wall and Ts_road simultaneously
  !        *************************************************************
  !
  !  Equation 12 from masson (2000)
  !
  !   Rn_w = abs_Rg_w 
  !  + emis_w                 *      SVF_w                * Rat            (1)

  !  - sigma * emis_w                                     * Ts_w**4 (t+dt) (2)

  !  + sigma * emis_w * emis_r *     SVF_w                * Ts_r**4 (t+dt) (3)

  !  + sigma * emis_w * emis_w * (1-2*SVF_w)              * Ts_w**4 (t+dt) (4)

  !  + sigma * emis_w * emis_w * (1-emis_w) * (1-2*SVF_w) * (1-2*SVF_w) * Ts_w**4 (t+dt) (7)

  !  + sigma * emis_w * emis_w * (1-emis_r) *      SVF_w  * (1-  SVF_r) * Ts_w**4 (t+dt) (8)

  !  + sigma * emis_w * emis_r * (1-emis_w) *      SVF_w  * (1-2*SVF_w) * Ts_r**4 (t+dt) (9)


  !  Equation 11 from masson (2000)

  !   Rn_r = abs_Rg_r
  !  - sigma * emis_r                                    * Ts_r**4 (t+dt) (2)
  !  +         emis_r                       *    SVF_r   * Rat
  !  + sigma * emis_r * emis_w              * (1-SVF_r)  * Ts_w**4 (t+dt) (3)
  !  + sigma * emis_r * emis_w * (1-emis_w) * (1-SVF_r)  * (1-2*SVF_w) * Ts_w**4 (t+dt) (5)

  !  + sigma * emis_r * emis_r * (1-emis_w) * (1-SVF_r)  *      SVF_w  * Ts_r**4 (t+dt) (6)
  !
  !  H_w  = rho Cp CH V ( Ts_w (t+dt) - Ta_canyon )
  !
  !  LE_w = rho Lv CH V ( qs_w (t+dt) - qa_canyon )
  !
  !  H_r  = rho Cp CH V ( Ts_r (t+dt) - Ta_canyon )
  !
  !  LE_r = rho Lv CH V ( qs_r (t+dt) - qa_canyon )
  !
  ! with again
  !                AC_can * Swall/Sroad * Twall + AC_can * Troad + AC_top * Ta + H_traffic/Cp/rho/Sroad
  !   Ta_canyon = --------------------------------------   (Eq 26, Masson, 2000)
  !                AC_can * Swall/Sroad         + AC_can         + AC_top
  !
  !
  !                 AC_can * delt_road * Hu_road * qsat(Troad) + AC_top * qa + LE_traffic/Lv/rho/Sroad
  !   qa_canyon = ------------------------------------    (Eq 27, Masson, 2000)
  !                 AC_can * delt_road                        + AC_top
  !
  !
  !
  !
  !    6 : computation of fluxes for each surface type
  !        *******************************************
  !
  !
  !    7 : averaging of the fluxes
  !        ***********************
  !
  !   This is done on the total exchange surface (roof + wall + road),
  !  which is bigger than the horizontal surface (roof+road), leading
  !  to bigger fluxes.
  !
  !   The fluxes due to industrial activity are directly added into the 
  !  atmosphere
  !
  !
  !    8 : road reservoir evolution
  !        ************************
  !
  !   The roof reservoir runoff goes directly into the road reservoir.
  !
  !   Runoff occurs for road reservoir (too much water), as well as drainage
  !   (evacuation system, typical time scale: 1 day)
  !
  !
  !
  !
  !!      
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original   23/01/98 
  !!                 13/05/2002 Edmilson Freitas: changing from modules to 
  !!                            normal subroutines. Inclusion of the tebconst.h
  !!                            file, which has the constants used by TEB.!
  !!                            Elimination of declarations that are not needed.
  !!                 09/01/2004 Edmilson Freitas: Elimination of all snow 
  !!                            treatment and inclusion of the scheme in BRAMS.
  !-----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  !
  use teb_vars_const
  use therm_lib, only: rslif        !RAMS function to calculate the saturation vapor pressure
  implicit none

  !
  !*      0.1    declarations of arguments
  !
  !INPUT/OUTPUT VARIABLES
  !
  real :: PTS_TOWN      ! town surface temperature
  real :: PEMIS_TOWN    ! town equivalent emissivity
  real :: PALB_TOWN     ! town equivalent albedo
  real :: PT_CANYON     ! canyon air temperature
  real :: PQ_CANYON     ! canyon air specific humidity
  real :: PU_CANYON     ! canyon hor. wind
  real :: PTS_ROOF      ! roof surface temperature (= to PT_ROOF(1))
  real :: PTS_ROAD      ! road surface temperature (= to PT_ROAD(1))
  real :: PTS_WALL      ! wall surface temperature (= to PT_WALL(1))
  real :: PTI_ROAD      ! road deep temperature
  real :: PTI_BLD       ! inside building temperature
  real, dimension(3) :: PT_ROOF     ! roof layers temperatures
  real, dimension(3) :: PT_ROAD     ! road layers temperatures
  real, dimension(3) :: PT_WALL     ! wall layers temperatures
  real :: PWS_ROOF      ! roof water reservoir
  real :: PWS_ROAD      ! road water reservoir
  real :: PPS           ! pressure at the surface (roof level)
  real :: PPA           ! pressure at the first atmospheric level
  real :: PEXNS         ! surface exner function (roof level)
  real :: PEXNA         ! exner function at the first atmospheric level
  ! at the lowest level
  real :: PTA           ! temperature at the lowest level
  real :: PQA           ! specific humidity
  ! at the lowest level
  real :: PVMOD         ! module of the horizontal wind
  real :: PRHOA         ! air density
  ! at the lowest level
  real :: PLW_RAD       ! atmospheric infrared radiation
  real :: PDIR_SW_RAD   ! incoming direct solar radiation
  ! on an horizontal surface
  real :: PSCA_SW_RAD   ! scattered incoming solar rad.
  real :: PTANZEN       !  of zenithal angle
  real :: PRR           ! rain rate
  real :: PH_TRAFFIC    ! anthropogenic sensible
  !                     ! heat fluxes due to traffic
  real :: PLE_TRAFFIC   ! anthropogenic latent
  !                     ! heat fluxes due to traffic
  real :: PH_INDUSTRY   ! anthropogenic sensible
  !                     ! heat fluxes due to factories
  real :: PLE_INDUSTRY  ! anthropogenic latent
  !                     ! heat fluxes due to factories
  real :: PZREF         ! reference height of the first
  ! atmospheric level (temperature)
  real :: PDIRCOSZW     ! Director Cosinus along z
  !                     ! directions at surface w-point
  real :: PUSLOPE       ! tangential wind component
  ! at first atmospheric level
  real :: PVSLOPE       ! tangential wind component
  ! at first atmospheric level
  real :: PTSTEP        ! time step
  real :: PZ0_TOWN      ! town roughness length
  ! for momentum
  real :: PBLD          ! fraction of buildings
  real :: PBLD_HEIGHT   ! buildings h
  real :: PWALL_O_HOR   ! wall surf. / hor. surf.
  real :: PCAN_HW_RATIO ! canyon    h/W
  real :: PALB_ROOF     ! roof albedo
  real :: PEMIS_ROOF    ! roof emissivity
  real, dimension(3)  :: PHC_ROOF      ! heat capacity for roof layers
  real, dimension(3)  :: PTC_ROOF      ! thermal conductivity for roof layers
  real, dimension(3)  :: PD_ROOF       ! depth of roof layers
  real :: PALB_ROAD      ! road albedo
  real :: PEMIS_ROAD     ! road emissivity
  real, dimension(3)  :: PHC_ROAD      ! heat capacity for road layers
  real, dimension(3)  :: PTC_ROAD      ! thermal conductivity for road layers
  real, dimension(3)  :: PD_ROAD       ! depth of road layers
  real :: PSVF_ROAD     ! road sky view factor
  real :: PALB_WALL     ! wall albedo
  real :: PEMIS_WALL    ! wall emissivity
  real, dimension(3)  :: PHC_WALL      ! heat capacity for wall layers
  real, dimension(3)  :: PTC_WALL      ! thermal conductivity for wall layers
  real, dimension(3)  :: PD_WALL       ! depth of wall layers
  real :: PSVF_WALL     ! wall sky view factor
  !
  !OUTPUT VARIABLES
  !
  real :: PRN_ROOF     ! net radiation over roof
  real :: PH_ROOF      ! sensible heat flux over roof
  real :: PLE_ROOF     ! latent heat flux over roof
  real :: PGFLUX_ROOF  ! flux through the roof
  real :: PRUNOFF_ROOF ! runoff over the ground
  real :: PRN_ROAD     ! net radiation over road
  real :: PH_ROAD      ! sensible heat flux over road
  real :: PLE_ROAD     ! latent heat flux over road
  real :: PGFLUX_ROAD  ! flux through the road
  real :: PRUNOFF_ROAD ! runoff over the ground
  real :: PRN_WALL     ! net radiation over wall
  real :: PH_WALL      ! sensible heat flux over wall
  real :: PLE_WALL     ! latent heat flux over wall
  real :: PGFLUX_WALL  ! flux through the wall
  !
  !
  real :: PRN_TOWN     ! net radiation over town
  real :: PH_TOWN      ! sensible heat flux over town
  real :: PLE_TOWN     ! latent heat flux over town
  real :: PGFLUX_TOWN  ! flux through the ground
  real :: PRUNOFF_TOWN ! runoff over the ground
  real :: PSFU_TOWN    ! U momentum flux over town
  real :: PSFV_TOWN    ! V momentum flux over town
  real :: PCH_TOWN     ! town averaged heat transfer
  !                    ! coefficient
  !
  !*      0.2    declarations of local variables
  !
  logical   :: GTI_EVOL       ! true --> internal temperature
  !                           !          of buildings evolves
  !                           ! false--> it is fixed
  !
  !
  ! Local variables for the calculation of QSATP and QSATZ
  real :: QSATP     ! qsat(Ta) previous
  real :: QSATZ     ! qsat(Ta)
  !
  real :: ZWS_ROOF_MAX   ! maximum deepness of roof
  real :: ZWS_ROAD_MAX   ! and road water reservoirs
  !
  !
  real :: ZDELT_ROOF     ! water fraction
  real :: ZDELT_ROAD     ! on the surface
  !
  real :: ZQSAT_ROAD     ! q_sat(Ts)
  real :: ZQSAT_ROOF     ! q_sat(Ts)
  real :: ZAC_ROOF       ! surface conductance
  !                      ! for heat transfers
  !                      ! above roofs
  real :: ZAC_ROOF_WAT   ! surface conductance
  !                      ! for heat transfers
  !                      ! above roofs (for water)
  real :: ZAC_WALL       ! surface conductance
  !                      ! for heat transfers
  !                      ! between wall and canyon air
  real :: ZAC_ROAD       ! surface conductance
  !                      ! for heat transfers
  !                      ! between road and canyon air
  real :: ZAC_ROAD_WAT   ! surface conductance
  !                      ! for heat transfers
  !                      ! inside canyon (for water)
  real :: ZAC_TOP        ! surface conductance
  !                      ! for heat transfers
  !                      ! canyon top and atm.
  real :: ZAC_BLD        ! surface conductance inside
  ! the building itself
  real :: ZCD            ! drag coefficient
  real :: ZTA            ! air temperature extrapolated at roof level
  real :: ZQA            ! air humidity extrapolated at roof level
  !
  real :: ZROOF          ! roof, wall and
  real :: ZWALL          ! road fractions
  real :: ZROAD          ! of exchange surf.

  !EDF**************************************************************************
  real :: ZTOTS_O_HORS   ! total canyon+roof surface   
  !verificar no teste se essa variavel eh necessaria
  !aparentemente, nao.
  !EDF**************************************************************************
  !
  !                      ! over horizontal surface
  real :: ZWALL_O_ROAD   ! wall surface over road surface
  !     
  ! absorbed solar and infra-red radiation by road, wall and roof
  !                                                      
  real :: ZABS_SW_ROAD
  real :: ZABS_SW_WALL
  real :: ZABS_SW_ROOF
  real :: ZABS_LW_ROAD
  real :: ZABS_LW_WALL
  real :: ZABS_LW_ROOF
  !
  !
  !-----------------------------------------------------------------------------
  !
  !*      1.     initializations
  !              ---------------
  !
  !*      1.1    water reservoirs
  !              ----------------
  !
  ZWS_ROOF_MAX =  1. ! (1mm) maximum deepness of roof water reservoir
  ZWS_ROAD_MAX =  1. ! (1mm) maximum deepness of road water reservoir
  !
  !*      1.2    surfaces relative fractions
  !              ---------------------------
  !
  ZTOTS_O_HORS= 1. + PWALL_O_HOR
  !
  ZROOF=PBLD        /ZTOTS_O_HORS 
  ZWALL=PWALL_O_HOR /ZTOTS_O_HORS
  ZROAD=(1.-PBLD)   /ZTOTS_O_HORS
  !
  ZWALL_O_ROAD = ZWALL / ZROAD
  !
  !*       1.3    option for internal temperature of buildings
  !               --------------------------------------------
  !
  GTI_EVOL = .true.
  !
  !-----------------------------------------------------------------------------
  !
  !*      2.     snow-covered surfaces relative effects
  !              --------------------------------------
  !               Not included in BRAMS version
  !          
  !
  !*      3.     SOLAR RADIATIONS
  !              ----------------
  !
  !* for sake of simplicity, the radiative properties of the vegetation (if any)
  !  situated in the canyon (and then influencing h/w) is replaced by the road
  !  properties. This influences only the radiative part due to reflections,
  !  which are less important when h/w is low (i.e. large streets with garden),
  !  than when h/w is large (from 0.5 and above approximately). This case
  !  (h/w large) occurs in city downtowns, and then does not occur in presence
  !  of significant vegetation area.
  !
  call URBAN_SOLAR_ABS(PDIR_SW_RAD, PSCA_SW_RAD, PTANZEN, &
       PBLD, PWALL_O_HOR, PCAN_HW_RATIO,                  &
       PALB_ROOF,                                         &
       PALB_ROAD, PSVF_ROAD, PALB_WALL, PSVF_WALL,        &
       ZABS_SW_ROOF, ZABS_SW_ROAD, ZABS_SW_WALL,          &
       PALB_TOWN                                          )
  !
  !-----------------------------------------------------------------------------
  !
  !*      4.     drag accounting for the effets of town temperature and humidity
  !              ---------------------------------------------------------------
  !
  call URBAN_DRAG(PTSTEP, PT_CANYON, PQ_CANYON,             &
       PTS_ROOF, PTS_ROAD, PTS_WALL,                        &
       PEXNS, PEXNA, PTA, PQA, PPS, PRHOA,                  &
       PZREF,  PDIRCOSZW, PVMOD,                            &
       PZ0_TOWN, PBLD, PBLD_HEIGHT, PCAN_HW_RATIO,          &
       ZWALL_O_ROAD,                                        &
       PWS_ROOF, PWS_ROAD,                                  &
       ZWS_ROOF_MAX, ZWS_ROAD_MAX,                          &
       ZQSAT_ROOF, ZQSAT_ROAD, ZDELT_ROOF, ZDELT_ROAD,      &
       ZCD, ZAC_ROOF, ZAC_ROOF_WAT,                         &
       ZAC_WALL, ZAC_ROAD, ZAC_ROAD_WAT, ZAC_TOP,           &
       PU_CANYON                                            )

  !
  !* area-averaged heat transfer coefficient
  !
  PCH_TOWN =        (       PBLD  * ZAC_ROOF      &
       + (1.-PBLD) * ZAC_TOP     ) &
       / PVMOD
  !
  !-----------------------------------------------------------------------------
  !
  !*      5.     Extrapolation of atmospheric T and q at roof level 
  !              (for fluxes computation)
  !              --------------------------------------------------
  !  Equations without number from Masson (2000). pg 373 item 2.9.4
  ZTA = PTA * PEXNS / PEXNA

  !
  !*****************************************************************************
  !ZQA = PQA * QSAT(ZTA,PPS) / QSAT(PTA,PPA)
  QSATP=rslif(PPS,ZTA)

  QSATZ=rslif(PPA,PTA)

  !
  !modified
  ZQA = PQA * QSATP / QSATZ
  !
  !-----------------------------------------------------------------------------
  !
  !*      6.     Evolution of interior building air temperature
  !              ----------------------------------------------
  !
  !
  call BLD_E_BUDGET(GTI_EVOL, PTSTEP, PBLD, PWALL_O_HOR,             &
       PRHOA, PT_ROOF, PT_WALL, PTI_BLD, ZAC_BLD        )
  !
  !-----------------------------------------------------------------------------
  !
  !*      7.     Snow mantel model
  !              -----------------
  !
  !* not used in BRAMS
  !
  !

  !-----------------------------------------------------------------------------
  !
  !*      8.     Roof Ts computation
  !              -------------------
  !
  !* ts_roof and qsat_roof are updated
  !
  call ROOF_LAYER_E_BUDGET(PTS_ROOF, PT_ROOF, ZQSAT_ROOF,&
       ZTA, ZQA, PPS,                                    &
       PLW_RAD, PTSTEP,                                  &
       PEMIS_ROOF, PHC_ROOF, PTC_ROOF, PD_ROOF,          &
       PTI_BLD, ZAC_BLD, ZDELT_ROOF,                     &
       PRHOA, ZAC_ROOF, ZAC_ROOF_WAT,                    &
       ZABS_SW_ROOF, ZABS_LW_ROOF                        )
  !
  !-----------------------------------------------------------------------------
  !
  !
  !*      9.    Road and wall Ts computations
  !             -----------------------------
  !
  !* ts_road, ts_wall, qsat_road, t_canyon and q_canyon are updated
  !
  !
  call ROAD_WALL_LAYER_E_BUDGET(PTS_ROAD, PT_ROAD, &
       PTS_WALL, PT_WALL, ZQSAT_ROAD,              &
       PT_CANYON, PQ_CANYON,                       &
       ZTA, ZQA, PPS,                              &
       PLW_RAD,  PTSTEP,                           &
       PH_TRAFFIC, PLE_TRAFFIC,                    &
       PBLD, ZWALL_O_ROAD,                         &
       PEMIS_ROAD, PSVF_ROAD,                      &
       PHC_ROAD, PTC_ROAD, PD_ROAD,                &
       PEMIS_WALL, PSVF_WALL,                      &
       PHC_WALL, PTC_WALL, PD_WALL,                &
       PTI_BLD, ZAC_BLD, PTI_ROAD,                 &
       ZDELT_ROAD,                                 &
       PRHOA, ZAC_WALL,                            &
       ZAC_ROAD, ZAC_ROAD_WAT, ZAC_TOP,            &
       ZABS_SW_ROAD, ZABS_SW_WALL,                 &
       ZABS_LW_ROAD, ZABS_LW_WALL                  )
  !
  !-----------------------------------------------------------------------------
  !
  !*     10.     Fluxes
  !              ------
  !
  !*      10.8.     Fluxes at snow-free roofs
  !              -------------------------
  !
  !                                            net radiation
  !
  PRN_ROOF = ZABS_SW_ROOF + ZABS_LW_ROOF
  !
  !                                            sensible heat flux
  !
  PH_ROOF =(PTS_ROOF - PTA) &
       * ZAC_ROOF * PRHOA * XCPD
  !
  !                                            latent heat of evaporation from
  !                                            the ground
  !
  PLE_ROOF =(ZQSAT_ROOF - PQA)                     &
       * ZAC_ROOF_WAT * ZDELT_ROOF * PRHOA * XLVTT
  !
  !-----------------------------------------------------------------------------
  !
  !*      10.9.     Fluxes at snow-free roads
  !              -------------------------
  !
  !                                            net radiation
  !
  PRN_ROAD = ZABS_SW_ROAD + ZABS_LW_ROAD
  !
  !                                            sensible heat flux
  !
  PH_ROAD =  (PTS_ROAD - PT_CANYON) * ZAC_ROAD * PRHOA * XCPD 
  !
  !                                            latent heat of evaporation from
  !                                            the ground
  !
  PLE_ROAD = ( ZQSAT_ROAD - PQ_CANYON )                   &
       *   ZAC_ROAD_WAT * ZDELT_ROAD * PRHOA * XLVTT
  !
  !
  !-----------------------------------------------------------------------------
  !
  !*     10.10.     Fluxes at walls
  !              ---------------
  !
  !                                            net radiation
  !
  PRN_WALL = ZABS_SW_WALL + ZABS_LW_WALL
  !
  !                                            sensible heat flux
  !
  PH_WALL =  (PTS_WALL - PT_CANYON) * ZAC_ROAD * PRHOA * XCPD
  !
  !                                            latent heat of evaporation from
  !                                            the ground
  !
  PLE_WALL = 0.
  !
  !                                            heat flux into the ground
  !
  PGFLUX_WALL = PRN_WALL - PH_WALL - PLE_WALL
  !
  !-----------------------------------------------------------------------------
  !
  !*     10.12.     Snow-free and snow-covered surfaces averaging
  !              ---------------------------------------------
  !
  !*     10.12.1    roads
  !              -----
  !
  !                                            heat flux into the ground
  !
  !
  PGFLUX_ROAD =(PRN_ROAD - PH_ROAD - PLE_ROAD )
  !
  !
  !
  !*     12.2    roofs
  !              -----
  !
  !                                            heat flux into the ground
  !
  !
  PGFLUX_ROOF  =   (PRN_ROOF - PH_ROOF - PLE_ROOF )
  !
  !
  !
  !-----------------------------------------------------------------------------
  !
  !*     10.13.     Momentum fluxes
  !              ---------------
  !
  !
  PSFU_TOWN = - ZCD * PVMOD*PUSLOPE
  !
  PSFV_TOWN = - ZCD * PVMOD*PVSLOPE
  !
  !-----------------------------------------------------------------------------
  !
  !*     10.14.     Averaged fluxes
  !              ---------------
  !
  PRN_TOWN= ZROOF * PRN_ROOF  * ZTOTS_O_HORS  &
       + ZROAD * PRN_ROAD  * ZTOTS_O_HORS     &
       + ZWALL * PRN_WALL  * ZTOTS_O_HORS
  !
  PH_TOWN = ZROOF * PH_ROOF   * ZTOTS_O_HORS  &
       + ZROAD * PH_ROAD   * ZTOTS_O_HORS     &
       + ZWALL * PH_WALL   * ZTOTS_O_HORS     &
       + PH_TRAFFIC                           &
       + PH_INDUSTRY 
  !
  PLE_TOWN= ZROOF * (PLE_ROOF) * ZTOTS_O_HORS &
       + ZROAD * (PLE_ROAD) * ZTOTS_O_HORS    &
       + ZWALL * (PLE_WALL) * ZTOTS_O_HORS    &
       + PLE_TRAFFIC                          &
       + PLE_INDUSTRY
  !
  !
  PGFLUX_TOWN=    ZROOF * PGFLUX_ROOF * ZTOTS_O_HORS &
       + ZROAD * PGFLUX_ROAD * ZTOTS_O_HORS          &
       + ZWALL * PGFLUX_WALL * ZTOTS_O_HORS
  !
  !-----------------------------------------------------------------------------
  !
  !*     10.15.     Averaged temperature (for radiation)
  !              ------------------------------------
  !
  PTS_TOWN = sqrt(sqrt( (                                           &
       PBLD                  * PEMIS_ROOF  * PTS_ROOF    **4        &
       + (1.-PBLD) *    PSVF_ROAD  * PEMIS_ROAD  * PTS_ROAD    **4  &
       + (1.-PBLD) *(1.-PSVF_ROAD) * PEMIS_WALL  * PTS_WALL    **4  &
       ) / PEMIS_TOWN                                               &
       )    )
  !


  !-----------------------------------------------------------------------------


  !*     11.     Roof and road reservoirs evolution
  !              ----------------------------------
  !
  call URBAN_HYDRO(ZWS_ROOF_MAX,ZWS_ROAD_MAX, PWS_ROOF, PWS_ROAD, &
       PRR, PTSTEP, PBLD, PLE_ROOF, PLE_ROAD,                     &
       PRUNOFF_ROOF,                                              &
       PRUNOFF_ROAD,                                              &
       PRUNOFF_TOWN                                               )
  !
  !-----------------------------------------------------------------------------
  !
  return
end subroutine urban
!
!   ##########################################################################
subroutine URBAN_SOLAR_ABS(PDIR_SW, PSCA_SW,  TZEN,&
     PBLD, PWALL_O_HOR, PCAN_HW_RATIO,             &
     PALB_ROOF,                                    &
     PALB_ROAD, PSVF_ROAD, PALB_WALL, PSVF_WALL,   &
     PABS_SW_ROOF, PABS_SW_ROAD, PABS_SW_WALL,     &
     PALB_TOWN                                     )
  !   ##########################################################################
  !
  !!****  *URBAN_SOLAR_ABS*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the solar radiation flux absorbed by roofs, roads and walls.
  !     The absorption by roofs is trivial.
  !         
  !     
  !!**  METHOD
  !     ------
  !
  !
  !        computation of input solar radiation on each surface
  !        ****************************************************
  !
  !    direct fluxes:
  !    -------------
  !    Equation 13 from masson(2000)
  !
  !    dir_Rg_road (Wm-2) =   S [ 2*theta0/pi
  !                              - 2/pi* h/W *tan(zen) * (1-cos(theta0))]
  !
  !    Equation 14 from masson(2000)
  !
  !    dir_Rg_wall (Wm-2) =   S [ W/h * (1/2 -theta0/pi) 
  !                              + tan(zen)/pi * (1-cos(theta0))]
  !
  !   where zen      is the zenithal angle, from zenit
  !         h/W      is the aspect ratio of the canyon
  !         S        is the direct solar radiation flux on a horizontal surface
  !
  !         theta0 = arcsin(min(W/h * 1./tan(zen),1))
  !
  !   The surfaces will keep (1-a) times these fluxes, and reflect the
  !   remaining
  !
  !    scattered fluxes:
  !    ----------------
  !
  !   sca_Rg_road = sca_Rg * SVF_road
  !
  !   sca_Rg_wall = sca_Rg * SVF_wall
  !
  !
  !    solar flux and isotropic reflections :
  !    ------------------------------------
  !
  !  after 0 reflection, the absorbed part of the flux is:
  !
  !      ARg_r(0) = (1-a_r) (sca_Rg_road + dir_Rg_road)
  !
  !      ARg_w(0) = (1-a_w) (sca_Rg_wall + dir_Rg_wall)
  !  
  !    and the reflected parts are
  !
  !      RRg_r(0) = a_r (sca_Rg_road + dir_Rg_road)
  !
  !      RRg_w(0) = a_w (sca_Rg_wall + dir_Rg_wall)
  !
  !  after n reflection:
  !
  !      ARg_r(n) = ARg_r(n-1) + RRg_w(n-1) * (1-  SVF_r)(1-a_r)
  !
  !      ARg_w(n) = ARg_w(n-1) + RRg_r(n-1) *      SVF_w (1-a_w)
  !                            + RRg_w(n-1) * (1-2*SVF_w)(1-a_w)
  !
  !      RRg_r(n) = (1- SVF_r) a_r RRg_w(n-1)
  !
  !      RRg_w(n) =     SVF_w  a_w RRg_r(n-1)
  !                +(1-2SVF_w) a_w RRg_w(n-1)
  !
  !
  !   i.e.
  !                                               n-1
  !      ARg_r(n) = ARg_r(0) + (1-  SVF_r)(1-a_r) SUM RRg_w(k)
  !                                               k=0
  !
  !                                               n-1
  !      ARg_w(n) = ARg_w(0) +      SVF_w (1-a_w) SUM RRg_r(k)
  !                                               k=0
  !                                               n-1
  !                          + (1-2*SVF_w)(1-a_w) SUM RRg_w(k)
  !                                               k=0
  !
  ! with
  !
  !     n                             n-1
  !    SUM RRg_r(k) = (1-  SVF_r) a_r SUM RRg_w(k)      +  RRg_r(0)
  !    k=0                            k=0
  !
  !     n                             n-1
  !    SUM RRg_w(k) =      SVF_w  a_w SUM RRg_r(k) 
  !    k=0                            k=0
  !                                   n-1
  !                  +(1-2*SVF_w) a_w SUM RRg_w(k)      +  RRg_w(0)
  !                                   k=0
  !
  !
  !   Then
  !
  !     n                                        n-1
  !    SUM RRg_w(k) =  (1-2*SVF_w)       a_w     SUM RRg_w(k)
  !    k=0                                       k=0
  !                                              n-2
  !                  + (1-  SVF_r) SVF_w a_w a_r SUM RRg_w(k) 
  !                                              k=0
  !
  !                  + RRg_w(0) + SVF_w a_w RRg_r(0)
  !
  !
  !
  !
  !  solving this system, lead after an infinity of reflections/absorptions:
  !
  !
  !    inf            RRg_r(0) + (1-  SVF_r) a_r ( RRg_w(0)  + a_w SVF_w RRg_r(0) )
  !    SUM RRg_r(k) = ----------------------------------------------------------
  !    k=0             1 - (1-2*SVF_w) a_w - (1-  SVF_r) SVF_w a_w a_r
  !
  !    inf                      RRg_w(0) + SVF_w a_w RRg_r(0)
  !    SUM RRg_w(k) = ----------------------------------------------------
  !    k=0             1 - (1-2*SVF_w) a_w - (1-  SVF_r) SVF_w a_w a_r
  !
  !
  ! ARg_r(n) and ARg_w(n) follow
  !
  !   ARg_r(n) = 
  !
  !
  !!    EXTERNAL
  !!    --------
  !!
  !!
  !!
  !!      
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original    23/01/98 
  !!                  21/11/00 (V. Masson) bug in reflections for roads
  !-----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  use teb_vars_const
  implicit none
  !
  !*      0.1    declarations of arguments
  !
  ! INPUT VARIABLES
  !
  real,intent(IN)  :: PDIR_SW           ! incoming direct solar radiation
  ! on an horizontal surface
  real,intent(IN)  :: TZEN              ! tangente of zenithal angle
  real,intent(IN)  :: PSCA_SW           ! scattered incoming solar rad.
  real,intent(IN)  :: PBLD              ! buildings fraction
  real,intent(IN)  :: PCAN_HW_RATIO     ! canyon    h/W
  real,intent(IN)  :: PWALL_O_HOR       ! wall surf. / hor. surf
  real,intent(IN)  :: PALB_ROOF         ! roof albedo
  real,intent(IN)  :: PALB_ROAD         ! road albedo
  real,intent(IN)  :: PSVF_ROAD         ! road sky view factor
  real,intent(IN)  :: PALB_WALL         ! wall albedo
  real,intent(IN)  :: PSVF_WALL         ! wall sky view factor
  !
  !OUTPUT VARIABLES
  !
  real,intent(OUT)  :: PABS_SW_ROOF       ! solar radiation absorbed
  !                                       ! by roofs
  real,intent(OUT)  :: PABS_SW_ROAD       ! solar radiation absorbed
  !                                       ! by roads
  real,intent(OUT)  :: PABS_SW_WALL       ! solar radiation absorbed
  !                                       ! by walls
  real,intent(OUT)  :: PALB_TOWN        ! town equivalent albedo
  !
  !*      0.2    declarations of local variables
  !
  !
  real  :: THETA0         ! canyon angle for
  !                       ! which solar
  !                       ! radiation
  !                       ! reaches the road
  !
  real  :: ZDIR_SW_ROAD
  real  :: ZDIR_SW_WALL
  real  :: ZREF0_SW_ROAD
  real  :: ZREF0_SW_WALL
  real  :: ZSREF_SW_ROAD
  real  :: ZSREF_SW_WALL
  real  :: ZROOF_SW   ! roof, wall and
  real  :: ZWALL_SW   ! road fractions of SW
  real  :: ZROAD_SW   ! interacting surf.
  !
  !-----------------------------------------------------------------------------
  !

  PABS_SW_ROOF = 0.
  PABS_SW_ROAD = 0.
  PABS_SW_WALL = 0.
  !
  !
  !*      1.     SOLAR RADIATIONS FOR ROOFS
  !              --------------------------
  !
  if (PDIR_SW+PSCA_SW>0.)then
     PABS_SW_ROOF      = (PDIR_SW+PSCA_SW) * (1. - PALB_ROOF  )
  else
     PABS_SW_ROOF      = 0.
  endif
  !
  !-----------------------------------------------------------------------------
  !
  !*      2.     SOLAR RADIATIONS FOR ROADS AND WALLS
  !              ------------------------------------
  !
  if (PCAN_HW_RATIO>0. .and. PDIR_SW+PSCA_SW>0.)then
     !
     !*      2.1    radiation coefficients
     !              ----------------------
     !
     !   Equation w/n just before eq 13 from Masson(2000)
     !
     THETA0 = asin(min(abs(1./TZEN)/PCAN_HW_RATIO,1.))
     !
     !
     !*      2.2    direct solar radiation received by roads
     !               ---------------------------------------
     !
     !   Equation 13 from Masson (2000)
     !
     ZDIR_SW_ROAD =  PDIR_SW * 2. * THETA0 / XPI                  &
          - PDIR_SW * 2. * abs(TZEN) /XPI              &
          * PCAN_HW_RATIO * (1.-cos(THETA0))
     !
     !*      2.3    direct solar radiation received by walls
     !              ----------------------------------------
     !
     !   Equation 14 from Masson (2000)
     !
     ! ZDIR_SW_WALL =  PDIR_SW * ABS(TZEN) /XPI                &
     !                                   * (1.-COS(ZTHETA0))                &
     !                    + PDIR_SW / PCAN_HW_RATIO *(0.5-ZTHETA0/XPI)
     !
     !   Equation 14 modified by using eq 13 from Masson (2000)
     !
     ZDIR_SW_WALL = ( PDIR_SW - ZDIR_SW_ROAD )  &
          * 0.5 / PCAN_HW_RATIO
     !
     !*      2.4    averaged albedos when snow is present (not the case)
     !              -------------------------------------
     !              Not needed as snow is never present
     !
     !*      2.5    first solar radiation reflection
     !              --------------------------------
     !
     !   Equation w/n just after eq 17 from Masson (2000). The scattered solar radiation received 
     !   by the surfaces (PSCA_SW) is directly deduced from the sky-view factors (PSVF_ROAD/PSVF_WALL).
     !
     ZREF0_SW_ROAD = PALB_ROAD * (PSCA_SW * PSVF_ROAD + ZDIR_SW_ROAD)
     !
     ZREF0_SW_WALL = PALB_WALL  * (PSCA_SW * PSVF_WALL + ZDIR_SW_WALL)
     !
     !*      2.6    sum of solar radiation reflected
     !              --------------------------------
     ! Equation 17 from Masson (2000). NOTE: there was a difference between the
     ! equation described in the paper and the one below. Following the paper,
     ! it was modified here. The original equation was keept to be checked later
     !
     !
     ! ZSREF_SW_WALL = (ZREF0_SW_WALL + PSVF_WALL * PALB_WALL * ZREF0_SW_ROAD) &
     !     / (1. - (1.-2.*PSVF_WALL) * PALB_WALL - (1.- PSVF_ROAD) * PSVF_WALL &
     !     * PALB_WALL * PALB_ROAD )
     !
     ZSREF_SW_WALL = (ZREF0_SW_WALL + PSVF_WALL * PALB_WALL * ZREF0_SW_ROAD) &
          / (1. - (1.-2.*PSVF_WALL) * PALB_WALL + (1.- PSVF_ROAD) * PSVF_WALL &
          * PALB_WALL * PALB_ROAD )
     !
     ! Equation 16 from Masson (2000). NOTE: The equation below was also different
     ! from the one in the paper. The original way in the programm was keept
     ! to be checked later
     !
     ! ZSREF_SW_ROAD = ( (1.- PSVF_ROAD) * PALB_ROAD * ZREF0_SW_WALL      &
     !                +(1.- PSVF_ROAD) * PALB_ROAD * PSVF_WALL          &
     !                * PALB_WALL * ZREF0_SW_ROAD          )             &
     ! / (1. - (1.-2.*PSVF_WALL) * PALB_WALL - (1.-   PSVF_ROAD) * PSVF_WALL &
     ! * PALB_WALL  * PALB_ROAD) + ZREF0_SW_ROAD

     ZSREF_SW_ROAD = ( (1.- PSVF_ROAD) * PALB_ROAD * ZREF0_SW_WALL      &
          +(1.- PSVF_ROAD) * PALB_ROAD * PSVF_WALL          &
          * PALB_WALL * ZREF0_SW_ROAD   + ZREF0_SW_ROAD)     &
          / (1. - (1.-2.*PSVF_WALL) * PALB_WALL +            &
          (1.-   PSVF_ROAD) * PSVF_WALL * PALB_WALL * PALB_ROAD)
     !
     !*      2.6    total solar radiation received by roads
     !              ---------------------------------------
     !
     !     Equation 18 from Masson(2000)
     !
     PABS_SW_ROAD = (1.-PALB_ROAD)  *     &         ! absorvity
          (ZDIR_SW_ROAD         &         ! direct solar radiation
          + PSCA_SW * PSVF_ROAD    &         ! scaterred solar radiation
          + ZSREF_SW_WALL * (1.- PSVF_ROAD)) ! reflected solar radiation
     !
     !
     !*      2.7    total solar radiation received by walls
     !              ---------------------------------------
     !
     !     Equation 19 from Masson (2000)
     !
     PABS_SW_WALL = (1.-PALB_WALL)                              &
          * (   ZDIR_SW_WALL                       &
          + PSCA_SW * PSVF_WALL                &
          + ZSREF_SW_WALL * (1.-2.*PSVF_WALL)  &
          + ZSREF_SW_ROAD * PSVF_WALL         ) 
     !
  endif
  !
  !-----------------------------------------------------------------------------
  !
  !*      3.     TRIVIAL CASES
  !              -------------
  !
  if(PCAN_HW_RATIO==0. .and. PDIR_SW+PSCA_SW>0.)then
     PABS_SW_ROAD      = (1.-PALB_ROAD  ) * ( PSCA_SW + PDIR_SW )
     PABS_SW_WALL      = 0.
  endif
  !
  !-----------------------------------------------------------------------------
  !
  !*      4.     ALBEDO
  !              ------
  !
  if (PDIR_SW + PSCA_SW >0.)then
     !
     ZROOF_SW = PBLD
     ZWALL_SW = PWALL_O_HOR
     ZROAD_SW = 1.- PBLD
     !
     PALB_TOWN = 1.                                    &
          - (  ZROOF_SW* (  PABS_SW_ROOF )      &
          + ZROAD_SW* (  PABS_SW_ROAD )      &
          + ZWALL_SW * PABS_SW_WALL    )     &
          !             /( ZROOF_SW + ZWALL_SW + ZROAD_SW    ) VERIFICAR 09/01/04 &
     /( PDIR_SW + PSCA_SW )
  else
     PALB_TOWN =            PBLD  * PALB_ROOF    &
          + (1. - PBLD) * PALB_ROAD    
  endif
  !
  !-----------------------------------------------------------------------------
  !
  return
end subroutine URBAN_SOLAR_ABS
!
!   ##########################################################################
subroutine URBAN_DRAG(PTSTEP, PT_CANYON, PQ_CANYON,    &
     PTS_ROOF, PTS_ROAD, PTS_WALL,                     &
     PEXNS, PEXNA, PTA, PQA, PPS, PRHOA,               &
     PZREF,  PDIRCOSZW, PVMOD,                         &
     PZ0_TOWN, PBLD, PBLD_HEIGHT, PCAN_HW_RATIO,       &
     PWALL_O_ROAD,                                     &
     PWS_ROOF, PWS_ROAD,                               &
     PWS_ROOF_MAX, PWS_ROAD_MAX,                       &
     PQSAT_ROOF, PQSAT_ROAD, PDELT_ROOF, PDELT_ROAD,   &
     PCD, PAC_ROOF, PAC_ROOF_WAT,                      &
     PAC_WALL, PAC_ROAD, PAC_ROAD_WAT, PAC_TOP,        &
     PU_CAN                                            )
  !   ##########################################################################
  !
  !!****  *URBAN_DRAG*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the surface drag over artificial surfaces as towns, 
  !     taking into account the canyon like geometry of urbanized areas.
  !         
  !     
  !!**  METHOD
  !!    ------
  !
  !!    REFERENCE
  !!    ---------
  !!
  !!      
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original    20/01/98 
  !       Edmilson Freitas (IAG-USP) 06/05/2001
  !-----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  !
  use teb_vars_const
  use therm_lib, only: rslif
  implicit none
  !*      0.1    declarations of arguments
  !
  !INPUT VARIABLES
  !
  real :: PTSTEP         ! time-step
  real :: PT_CANYON      ! canyon air temperature
  real :: PQ_CANYON      ! canyon air specific humidity.
  real :: PTS_ROOF       ! surface temperature
  real :: PTS_ROAD       ! surface temperature
  real :: PTS_WALL       ! surface temperature
  real :: PEXNS          ! surface exner function
  real :: PTA            ! temperature at the lowest level
  real :: PQA            ! specific humidity
  ! at the lowest level
  real :: PVMOD          ! module of the horizontal wind
  real :: PPS            ! pressure at the surface
  real :: PEXNA          ! exner function
  ! at the lowest level
  real :: PRHOA          ! air density
  real :: PZREF          ! reference height of the first
  ! atmospheric level (temperature)
  real :: PDIRCOSZW      ! Director Cosinus along z
  !                      ! directions at surface w-point
  real :: PZ0_TOWN       ! roughness length for momentum
  real :: PBLD           ! fraction of buildings
  real :: PBLD_HEIGHT    ! h
  real :: PCAN_HW_RATIO  ! h/W
  real :: PWALL_O_ROAD   ! wall surf. / road surf.
  !
  real :: PWS_ROOF       ! roof water content (kg/m2)
  real :: PWS_ROAD       ! road water content (kg/m2)
  real :: PWS_ROOF_MAX   ! maximum deepness of roof
  real :: PWS_ROAD_MAX   ! and water reservoirs (kg/m2)
  !
  !OUTPUT VARIABLES
  !
  real :: PQSAT_ROOF     ! qsat(Ts)
  real :: PQSAT_ROAD     ! qsat(Ts)
  real :: PDELT_ROOF     ! water fraction on
  real :: PDELT_ROAD     ! snow-free surfaces
  real :: PCD            ! drag coefficient
  real :: PAC_ROOF       ! aerodynamical conductance
  real :: PAC_ROOF_WAT   ! aerodynamical conductance (for water)
  real :: PAC_WALL       ! aerodynamical conductance
  !                      ! between canyon air and
  !                      ! walls 
  real :: PAC_ROAD       ! aerodynamical conductance
  !                      ! between canyon air and
  !                      ! roads
  real :: PAC_ROAD_WAT   ! aerodynamical conductance
  !                      ! between canyon air and
  !                      ! road (for water)
  real :: PAC_TOP        ! aerodynamical conductance
  !                      ! between canyon top and atm.
  real :: PU_CAN         ! hor. wind in canyon
  !
  !*      0.2    declarations of local variables
  !
  real :: ZZ0_O_Z0H      ! urban value for ratio between
  !                      ! z0 and z0h (Voogt & Grimmond 2000) =200.

  !***************************************************************************!  
  !
  real :: ZTS_TOWN       ! town averaged temp.
  real :: ZQ_TOWN        ! town averaged hum.
  real :: ZAVDELT_ROOF   ! averaged water frac.
  real :: ZQ_ROOF        ! roof spec. hum.
  real :: ZZ0_ROOF       ! roof roughness length
  real :: ZZ0_ROAD       ! road roughness length
  real :: ZW_CAN         ! ver. wind in canyon
  real :: ZRI            ! Richardson number
  real :: ZCDN           ! neutral drag coefficient
  real :: ZLE_MAX        ! maximum latent heat flux available
  real :: ZLE            ! actual latent heat flux
  real :: zwake
  !
  !-----------------------------------------------------------------------------
  !
  !
  ZZ0_O_Z0H = 200.
  !
  ZZ0_ROOF    = 0.15     ! z0 for roofs is equal to 01cm
  ZZ0_ROAD    = 0.05     ! z0 for roads
  !
  !-----------------------------------------------------------------------------
  !
  !*      1.     roof and road saturation specific humidity
  !              ------------------------------------------
  !*       1.    COMPUTE SATURATION VAPOR PRESSURE
  !              ---------------------------------
  !
  !*       2.    COMPUTE SATURATION HUMIDITY
  !              ---------------------------
  PQSAT_ROOF =rslif(PPS,PTS_ROOF)
  !
  !
  PQSAT_ROAD =rslif(PPS,PTS_ROAD)
  !
  !-----------------------------------------------------------------------------
  !
  !*      2.     fraction of water on roofs
  !              --------------------------
  !
  !
  !*      2.1    general case
  !              ------------
  !
  ! The expression below is the snow free fraction of water occupied by 
  ! liquid water given by Noilham and Planton (1989) where PWS_ROOF_MAX is
  ! the maximum water amount on the roof (see Masson (2000), pg 366 2nd paragraph.
  !
  PWS_ROOF_MAX=1.
  if (PQSAT_ROOF >= PQA )then
     PDELT_ROOF = (PWS_ROOF/PWS_ROOF_MAX)**(2./3.)
     !
  else
     !*      2.2    dew deposition on roofs (PDELT_ROOF=1)
     !              -----------------------
     PDELT_ROOF=1.
     !
  endif
  !-----------------------------------------------------------------------------
  !
  !*      3.     fraction of water on roads
  !              --------------------------
  !
  !
  !*      3.1    general case
  !              ------------
  ! The expression below is the snow free fraction of water occupied by 
  ! liquid water given by Noilham and Planton (1989) where PWS_ROAD_MAX is
  ! the maximum water amount on the road (see Masson (2000), pg 366 2nd paragraph.
  !
  PWS_ROAD_MAX=1.
  if (PQSAT_ROAD >= PQ_CANYON )then
     PDELT_ROAD = (PWS_ROAD/PWS_ROAD_MAX)**(2./3.)


  else        !It means that the surface is wet (see Masson(2000) pg. 366

     PDELT_ROAD=1.

  endif
  !
  !
  !-----------------------------------------------------------------------------
  !
  !*      4.     Drag coefficient for momentum between roof level and atmosphere
  !              ---------------------------------------------------------------
  !
  !
  !*      4.1    Averaged temperature at roof level
  !              ----------------------------------
  !
  ZTS_TOWN = PBLD * PTS_ROOF + (1.-PBLD) * PT_CANYON
  !
  !*      4.2    Averaged water fraction on roofs
  !              -------------------------------
  !
  ZAVDELT_ROOF = PDELT_ROOF   !Eliminar variavel mais tarde
  !
  !*      4.3    Roof specific humidity
  !              ----------------------
  !***************************************************************************
  ! Edmilson 08/05/2002
  ! the equation below was modified since PQSAT_ROOF was already calculated up
  !
  !ZQ_ROOF = QSAT(PTS_ROOF,PPS) * ZAVDELT_ROOF
  !
  ZQ_ROOF = PQSAT_ROOF * ZAVDELT_ROOF
  !***************************************************************************!
  !*      4.4    Averaged Saturation specific humidity
  !              -------------------------------------
  !
  ZQ_TOWN =       PBLD  * ZQ_ROOF &
       + (1.-PBLD) * PQ_CANYON 
  !
  !
  !-----------------------------------------------------------------------------
  !
  !*      5.     Momentum drag coefficient
  !              -------------------------
  !
  ! Mofification in calling the subroutine below. This subroutine is now located
  ! at the end of the subroutine urban_drag (this one)
  ! The subroutines are:SURFACE_RI,SURFACE_CD and SURFACE_AERO_COND
  !
  call SURFACE_RI(ZTS_TOWN, ZQ_TOWN, PEXNS, PEXNA, PTA, PQA,        &
       PZREF + PBLD_HEIGHT/3.,    &
       PDIRCOSZW, PVMOD, ZRI                             )
  !
  call SURFACE_CD(ZRI, PZREF + PBLD_HEIGHT/3., &
       PZ0_TOWN, PZ0_TOWN/ZZ0_O_Z0H, PCD, ZCDN             )
  !
  !-----------------------------------------------------------------------------
  !
  !*      6.     Drag coefficient for heat fluxes between roofs and atmosphere
  !              -------------------------------------------------------------
  !
  call SURFACE_RI(PTS_ROOF, ZQ_ROOF, PEXNS, PEXNA, PTA, PQA,        &
       PZREF, PDIRCOSZW, PVMOD, ZRI               )
  !
  call SURFACE_AERO_COND(ZRI, PZREF, PVMOD, ZZ0_ROOF,        &
       ZZ0_ROOF/ZZ0_O_Z0H, PAC_ROOF               )
  !
  !
  ZLE_MAX     = PWS_ROOF / PTSTEP * XLVTT
  ZLE         =(PQSAT_ROOF - PQA)                     &
       * PAC_ROOF * PDELT_ROOF * XLVTT * PRHOA
  !
  PAC_ROOF_WAT = PAC_ROOF
  !
  if (PDELT_ROOF==0.)PAC_ROOF_WAT=0.
  !
  if (ZLE>0.)PAC_ROOF_WAT = PAC_ROOF * min ( 1. , ZLE_MAX/ZLE )
  !-----------------------------------------------------------------------------
  !
  !*      7.     Drag coefficient for heat fluxes between canyon and atmosphere
  !              --------------------------------------------------------------
  !

  call SURFACE_RI(PT_CANYON, PQ_CANYON, PEXNS, PEXNA, PTA, PQA,        &
       PZREF + PBLD_HEIGHT/3.,                              &
       PDIRCOSZW, PVMOD, ZRI                                )
  !
  !Chamada abaixo estava diferente da versao atualizada
  call SURFACE_AERO_COND(ZRI, PZREF, PVMOD, PZ0_TOWN,                  &
       PZ0_TOWN/ZZ0_O_Z0H, PAC_TOP                   )
  !
  !CALL SURFACE_AERO_COND(ZRI, PZREF, PVMOD, PZ0_TOWN,                  &
  !                       PZ0_TOWN, PAC_TOP                             )
  !
  !
  !-----------------------------------------------------------------------------
  !
  !*      8.     Drag coefficient for heat fluxes between walls, road and canyon
  !              ---------------------------------------------------------------
  !
  !
  !*      8.1    Horizontal wind speed in canyon
  !              -------------------------------
  !
  !* skimming flow for h/w>1 (maximum effect of direction on wind in the canyon);
  !* isolated flow for h/w<0.5 (wind is the same in large streets for all dir.)
  !* wake flow between.
  !
  ZWAKE= 1. + (2./XPI-1.) * 2. * (PCAN_HW_RATIO-0.5)
  ZWAKE= max(min(ZWAKE,1.),2./XPI)
  !
  !* Estimation of canyon wind speed from wind just above roof level
  !  (at 1.33h). Wind at 1.33h is estimated using the log law.

  PU_CAN =  ZWAKE*exp(-PCAN_HW_RATIO/4.) * PVMOD                 &
       * log( (             2* PBLD_HEIGHT/3.) / PZ0_TOWN)     &
       / log( (PZREF      + 2* PBLD_HEIGHT/3.) / PZ0_TOWN)
  !
  !*      8.2    Vertical wind speed in canyon = ustar
  !              -----------------------------
  !
  ZW_CAN = sqrt (PCD) * PVMOD
  !
  !*      8.3    aerodynamical conductance for roads or walls
  !              --------------------------------------------
  !
  PAC_WALL = ( 11.8 + 4.2 * sqrt(PU_CAN**2 + ZW_CAN**2) ) &
       / XCPD / PRHOA
  !
  call SURFACE_RI(PTS_ROAD , PQ_CANYON, PEXNS, PEXNS, &
       PT_CANYON, PQ_CANYON,               &
       PBLD_HEIGHT/2.,      &
       PDIRCOSZW, PU_CAN, ZRI              )
  !
  call SURFACE_AERO_COND(ZRI,  PBLD_HEIGHT/2.,  &
       PU_CAN, ZZ0_ROAD,                     &
       ZZ0_ROAD/ZZ0_O_Z0H, PAC_ROAD          )
  !
  !PAC_CAN =  ( PWALL_O_ROAD * ZAC_WALL + ZAC_ROAD ) &
  !            / ( PWALL_O_ROAD               + 1.          )
  !
  !
  !*      8.4    aerodynamical conductance for water limited by available water
  !              --------------------------------------------------------------
  !
  ZLE_MAX     = PWS_ROAD / PTSTEP * XLVTT
  ZLE         = ( PQSAT_ROAD - PQ_CANYON )                   &
       *   PAC_ROAD * PDELT_ROAD * XLVTT * PRHOA
  !
  PAC_ROAD_WAT = PAC_ROAD
  !
  if (PDELT_ROAD==0.)PAC_ROAD_WAT= 0.
  !
  if (ZLE>0.) PAC_ROAD_WAT = PAC_ROAD * min ( 1. , ZLE_MAX/ZLE )
  !
  !-----------------------------------------------------------------------------
  !
  return
end subroutine URBAN_DRAG

!   ######################################################################
subroutine SURFACE_RI(PTS, PQS, PEXNS, PEXNA, PTA, PQA,   &
     PZREF, PDIRCOSZW, PVMOD, PRI )
  !   ######################################################################
  !
  !!****  *SURFACE_RI*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the richardson number near the ground
  !!      
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original    22/09/98 
  !!   Inclusion in the subroutine urban_drag by Edmilson Freitas
  !!   08/05/2002
  !-----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  !
  use teb_vars_const
  implicit none
  !
  !*      0.1    declarations of arguments
  !
  !
  !INPUT VARIABLES
  !
  real    :: PTS      ! surface temperature
  real    :: PQS      ! surface specific humidity
  real    :: PEXNS    ! surface exner function
  real    :: PTA      ! temperature at the lowest level
  real    :: PQA      ! specific humidity
  ! at the lowest level
  real    :: PEXNA    ! exner function
  ! at the lowest level
  real    :: PVMOD    ! module of the horizontal wind
  !
  real    :: PZREF    ! reference height of the first
  ! atmospheric level
  real    :: PDIRCOSZW! Cosine of the angle between
  !                   ! the normal to the surface and
  !                   ! the vertical
  !
  ! OUTPUT VARIABLES
  real    :: PRI      ! Richardson number
  !
  !*      0.2    declarations of local variables
  !
  !
  real    :: ZTHVA, ZTHVS
  !-----------------------------------------------------------------------------
  !
  !       1.     Richardson number
  !              -----------------
  !                
  !                                                 virtual potential        
  !                                                 temperature at the 
  !                                                 first atmospheric level and
  !                                                 at the surface
  !
  ZTHVA=PTA/PEXNA*( 1.+(XRV/XRD-1.)*PQA )   
  ZTHVS=PTS/PEXNS*( 1.+(XRV/XRD-1.)*PQS)
  !                                                 
  !
  ! Richardson's number
  PRI = XG * PDIRCOSZW * PZREF 		  &
       * (ZTHVA-ZTHVS) / (0.5 * (ZTHVA+ZTHVS) )  &
       / (PVMOD*PVMOD) 


  !PRI = XG * PDIRCOSZW * PZREF * PZREF              &
  !        * (ZTHVA-ZTHVS) / (0.5 * (ZTHVA+ZTHVS) )  &
  !        / (PVMOD*PVMOD) /PZREF
  !
  !-----------------------------------------------------------------------------
  !
  !
  return
end subroutine SURFACE_RI
!
!   #################################################################
subroutine SURFACE_CD(PRI, PZREF,  PZ0EFF, PZ0H,   &
     PCD, PCDN)
  !   #################################################################
  !
  !!****  *SURFACE_CD*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the drag coefficients for momentum near the ground
  !         
  !     
  !!**  METHOD
  !!    ------
  !
  !
  !
  !    1 and 2 : computation of relative humidity near the ground
  !
  !    3 : richardson number
  !
  !    4 : the aerodynamical resistance for heat transfers is deduced
  !
  !    5 : the drag coefficient for momentum ZCD is computed
  !
  !
  !!      
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original    20/01/98 
  !!                  02/04/01 (P Jabouille) limitation of Z0 with 0.5 PZREF
  !!                  08/05/02  Edmilson Freitas: Inclusion in the subroutine
  !!                            urban_drag.f90. Elimination of the declarations
  !!                            DIMENSION, INTENT. Inclusion of the tebconst.h.
  !!                            Elimination of the functions CMSTAR and PM.
  !-----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  use teb_vars_const
  implicit none
  !
  !*      0.1    declarations of arguments
  !
  !
  !INPUT VARIABLES
  real   :: PRI      ! Richardson number
  real   :: PZREF    ! reference height of the first
  ! atmospheric level
  real   :: PZ0EFF   ! roughness length for momentum
  ! with subgrid-scale orography
  real   :: PZ0H     ! roughness length for heat
  !

  !OUTPUT VARIABLES
  real   :: PCD      ! drag coefficient for momentum
  real   :: PCDN     ! neutral drag coefficient for momentum
  !
  !*      0.2    declarations of local variables
  !
  !
  real :: ZZ0EFF, ZZ0H, ZMU,ZCMSTAR, ZPM, ZCM, ZFM
  !-----------------------------------------------------------------------------
  !
  !*       1.     Drag coefficient for momentum transfers
  !               ---------------------------------------
  !
  ZZ0EFF = min(PZ0EFF,PZREF*0.5)
  ZZ0H   = min(ZZ0EFF,PZ0H)
  !
  ZMU = log( min(ZZ0EFF/ZZ0H,200.) )
  !
  PCDN = (XKARMAN/log(PZREF/ZZ0EFF))**2   !(a^2) Eq. 6 from Marcart et al 1995
  !it is the drag coefficient for neutral con
  !conditions
  !
  !
  ZCMSTAR=6.8741 + 2.6933*ZMU - 0.3601*ZMU*ZMU + 0.0154*ZMU*ZMU*ZMU
  !
  ZPM = 0.5233 - 0.0815*ZMU + 0.0135*ZMU*ZMU - 0.0010*ZMU*ZMU*ZMU
  !
  ZCM = 10.*ZCMSTAR*PCDN*( PZREF/ZZ0EFF )**ZPM    !10 is probably the parameter b
  !in equation 12 from 
  !Mascart et al 1995
  !PCDN is the parameter a^2
  !
  if( PRI > 0.0 )then
     ZFM = 1. + 10.*PRI / sqrt( 1.+5.*PRI )
     ZFM = 1. / ZFM
  else
     ZFM = 1. - 10.*PRI / ( 1.+ZCM*sqrt(-PRI) )
  endif
  !
  PCD = PCDN*ZFM
  !
  return
end subroutine SURFACE_CD
!
!RETURN
!END
!   ######################################################################
subroutine SURFACE_AERO_COND(PRI, PZREF, PVMOD, PZ0,&
     PZ0H, PAC                     )
  !   ######################################################################
  !
  !!****  *SURFACE_AERO_COND*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the drag coefficients for heat and momentum near the ground
  !         
  !     
  !!**  METHOD
  !!    ------
  !
  !
  !
  !    1 and 2 : computation of relative humidity near the ground
  !
  !    3 : richardson number
  !
  !    4 : the aerodynamical resistance for heat transfers is deduced
  !
  !    5 : the drag coefficient for momentum ZCD is computed
  !
  !
  !!      
  !!    REFERENCE
  !!    ---------
  !!
  !!      
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original    20/01/98 
  !!                  02/04/01 (P Jabouille) limitation of Z0 with 0.5 PZREF
  !!                  08/05/02  Edmilson Freitas: Inclusion in the subroutine
  !!                            urban_drag.f90. Elimination of the declarations
  !!                            DIMENSION, INTENT. Inclusion of the tebconst.h.
  !!                            Elimination of the functions CHSTAR and PH.
  !-----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  !
  use teb_vars_const
  implicit none
  !
  !*      0.1    declarations of arguments
  !
  !INPUT VARIABLES
  real    :: PRI      ! Richardson number
  real    :: PVMOD    ! module of the horizontal wind
  real    :: PZREF    ! reference height of the first
  ! atmospheric level
  real    :: PZ0      ! roughness length for momentum
  real    :: PZ0H     ! roughness length for heat
  !
  !OUTPUT VARIABLES
  real    :: PAC      ! aerodynamical conductance
  !
  !*      0.2    declarations of local variables
  !
  !
  real    :: ZZ0, ZZ0H,ZMU,ZFH,ZCHSTAR, ZPH,ZCDN,ZSTA,ZDI
  !-----------------------------------------------------------------------------
  !
  !*       4.     Surface aerodynamic resistance for heat transfers
  !               -------------------------------------------------
  !
  ZZ0 = min(PZ0,PZREF*0.5)
  ZZ0H = min(ZZ0,PZ0H)
  !
  !Equation 13 from Marcart. et al. 1995

  ZMU = max( log( ZZ0/ZZ0H ), 0.0 )

  !Equation without number from Mascart. et al. 1995 (pg 336)

  ZFH = log( PZREF/ZZ0 )/ log( PZREF/ ZZ0H )
  !
  !Equation 16 from Marcart. et al. 1995

  ZCHSTAR = 3.2165 + 4.3431*ZMU + 0.5360*ZMU*ZMU - 0.0781*ZMU*ZMU*ZMU
  !
  !Equation 17 from Marcart. et al. 1995

  ZPH = 0.5802 - 0.1571*ZMU + 0.0327*ZMU*ZMU - 0.0026*ZMU*ZMU*ZMU
  ! 

  ZCDN = (XKARMAN/log(PZREF/ZZ0))**2.

  ZSTA = PRI*PVMOD*PVMOD
  !
  !
  if ( PRI < 0.0 )then
     ZDI= 1. / ( PVMOD +ZCHSTAR*ZCDN*15.*(PZREF/ZZ0H)**ZPH*ZFH * sqrt(-ZSTA) )
     PAC = ZCDN*(PVMOD-15.*ZSTA*ZDI)*ZFH
  else
     ZDI = sqrt(PVMOD * PVMOD + 5. * ZSTA )
     PAC = ZCDN*PVMOD/(1.+15.*ZSTA*ZDI / PVMOD /PVMOD /PVMOD )*ZFH
  endif

  !
  !
  return
  !   ######################################################################
end subroutine SURFACE_AERO_COND
!   ######################################################################

!   ##########################################################################
subroutine BLD_E_BUDGET(OTI_EVOL, PTSTEP, PBLD, PWALL_O_HOR,           &
     PRHOA, PT_ROOF, PT_WALL, PTI_BLD, PAC_BLD      )
  !   ##########################################################################
  !
  !!****  *BLD_E_BUDGET*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the evoultion of the temperature of inside building air
  !         
  !     
  !!**  METHOD
  !     ------
  !
  !     The resistance term between the surfaces and the room is given
  !     by a standard value, which mimics both the convection
  !     and the radiative interactions in the room.
  !     This explains the very low resistance. It is used to compute
  !     the evolution of the surfaces only.
  !     This resistance value is 0.123 Km/W  (typical for inside surfaces).
  !     (ENVIRONMENTAL SCIENCE IN BUILDING, 3rd Edition, Randall McMullan,
  !      THE MACMILLAN PRESS Limited).
  !
  !
  !
  !     On the contrary, the evolution of the air temperature is mainly
  !     governed by the convection (considering the low radiative absorption
  !     of the air itself).
  !     In order to have a simple formulation, a diurnal cycle is assumed,
  !     with a force restore formulation.
  !
  !     The floor temperature is fixed
  !
  !!      
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original    24/08/00 
  !-----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  !
  use teb_vars_const
  implicit none
  !
  !*      0.1    declarations of arguments
  !
  ! INPUT VARIABLES
  logical :: OTI_EVOL      ! true --> internal temp. of
  !                        !      of buildings evolves
  !                        ! false--> it is fixed
  real  :: PTSTEP        ! time step
  real  :: PBLD          ! building fraction
  real  :: PWALL_O_HOR   ! wall surf. / hor. surf.
  real  :: PRHOA         ! air density at the lowest level
  real, dimension(3)  :: PT_ROOF       ! roof layers temperatures
  real, dimension(3)  :: PT_WALL       ! wall layers temperatures
  !
  !INPUT AND OUTPUT VARIABLE
  real  :: PTI_BLD       ! building air temperature

  !OUTPUT VARIABLE
  real  :: PAC_BLD       ! aerodynamical conductance
  ! inside building itself
  !
  !*      0.2    declarations of local variables
  !
  real  :: ZT_FLOOR     ! floor temperature
  !
  real  :: ZTAU         ! temporal filter period
  !
  integer :: IROOF        ! number of roof layers
  integer :: IWALL        ! number of wall layers
  !-----------------------------------------------------------------------------
  !
  !*      1.   initializations
  !            ---------------
  !
  IROOF = size(PT_ROOF)
  IWALL = size(PT_WALL)
  !
  ZT_FLOOR  = TMINBLD + XTT
  !
  !
  !*      2.   inside conductance FOR SURFACES
  !            ------------------
  !
  !* (normalized by rho Cp for convenience)
  !
  PAC_BLD = 1. / 0.123 / (XCPD * PRHOA)
  !
  !*      3.   no evolution of interior temperature if OTI_EVOL=.FALSE.
  !            --------------------------------------------------------
  !
  if (.not. OTI_EVOL) return
  !
  !*      4.   evolution of the internal temperature
  !            -------------------------------------
  !
  ZTAU=XDAY
  !
  PTI_BLD = PTI_BLD * (ZTAU-PTSTEP)/ZTAU                    &
       + (   PT_ROOF(IROOF) * PBLD                   &
       + PT_WALL(IWALL) * PWALL_O_HOR            &
       + ZT_FLOOR       * PBLD                )  &
       /(  2. * PBLD  +  PWALL_O_HOR ) * PTSTEP / ZTAU
  !
  !
  !*      5.   internal temperature set to a minimum value (heating)
  !            -----------------------------------------------------
  !
  PTI_BLD = max( PTI_BLD , TMINBLD + XTT )
  !
  !-----------------------------------------------------------------------------
  !
  return
end subroutine BLD_E_BUDGET

!   ##########################################################################
subroutine  ROOF_LAYER_E_BUDGET(PTS_ROOF, PT_ROOF, ZQSAT_ROOF,            &
     ZTA, ZQA, PPS,                            &
     PLW_RAD, PTSTEP,                          &
     PEMIS_ROOF, PHC_ROOF, PTC_ROOF, PD_ROOF,  &
     PTI_BLD, ZAC_BLD, ZDELT_ROOF,             &
     PRHOA, ZAC_ROOF, ZAC_ROOF_WAT,            &
     ZABS_SW_ROOF, ZABS_LW_ROOF                )
  !
  !!****  *ROOF_LAYER_E_BUDGET*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the evolution of surface temperature of roofs
  !         
  !     
  !!**  METHOD
  !     ------
  !
  !
  !
  !
  !    5 : equation for evolution of Ts_roof
  !        *********************************
  !
  !     dTt_1(t) / dt = 1/(dt_1*Ct_1) * (  Rn - H - LE
  !                                      - 2*Kt_1*(Tt_1-Tt_2)/(dt_1 +dt_2)       )
  !
  !     dTt_k(t) / dt = 1/(dt_k*Ct_k) * (- 2*Kt_k-1*(Tt_k-Tt_k-1)/(dt_k-1 +dt_k) 
  !                                      - 2*Kt_k  *(Tt_k-Tt_k+1)/(dt_k+1 +dt_k) )
  !
  !       with
  !
  !       K*_k  = (d*_k+ d*_k+1)/(d*_k/k*_k+ d*_k+1/k*_k+1)
  !
  !       Rn = (dir_Rg + sca_Rg) (1-a) + emis * ( Rat - sigma Ts**4 (t+dt) )
  !
  !       H  = rho Cp CH V ( Ts (t+dt) - Tas )
  !
  !       LE = rho Lv CH V ( qs (t+dt) - qas )
  !
  !      where the as subscript denotes atmospheric values at ground level
  !      (and not at first half level)
  !
  !      The tridiagonal system is linearized with
  !
  !       using      Ts**4 (t+dt) = Ts**4 (t) + 4*Ts**3 (t) * ( Ts(t+dt) - Ts(t) )
  !
  !       and  qs (t+dt) = Hu(t) * qsat(t) + Hu(t) dqsat/dT * ( Ts(t+dt) - Ts(t) )
  !
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original    23/01/98 
  !!                  13/05/2002  Edmilson Freitas (IAG-USP): Change of the module
  !!                              to a normal subroutine. Inclusion the tebconst.h
  !!                              file, which has the constants needed by TEB.
  !!                              Elimination of the declarations that are not
  !!                              needed.
  !-----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  !
  use teb_vars_const
  use therm_lib, only : rslif,rslifp
  implicit none
  !
  !*      0.1    declarations of arguments
  !
  !INPUT / OUTPUT VARIABLES

  real :: PTS_ROOF       ! roof surface temperature
  real, dimension(3) :: PT_ROOF      ! roof layers temperatures
  real :: ZQSAT_ROOF     ! q_sat(Ts)

  ! Only INPUT VARIABLES
  real  :: ZTA            ! air temperature at roof level
  real  :: ZQA            ! air specific humidity
  ! at roof level
  real  :: PPS            ! pressure at the surface
  real  :: PLW_RAD        ! atmospheric infrared radiation
  real  :: PTSTEP         ! time step
  real  :: PEMIS_ROOF     ! roof emissivity
  real, dimension(3) :: PHC_ROOF       ! heat capacity for roof layers
  real, dimension(3) :: PTC_ROOF       ! thermal conductivity for roof layers
  real, dimension(3) :: PD_ROOF        ! depth of roof layers
  real  :: PTI_BLD        ! inside building temp.
  real  :: ZAC_BLD        ! aerodynamical resistance
  ! inside building itself
  real  :: ZDELT_ROOF     ! fraction of water
  real  :: PRHOA          ! air density
  real  :: ZAC_ROOF       ! aerodynamical conductance
  real  :: ZAC_ROOF_WAT   ! aerodynamical conductance (for water)
  real  :: ZABS_SW_ROOF   ! absorbed solar radiation

  ! Only OUTPUT VARIABLES
  real  :: ZABS_LW_ROOF   ! absorbed infra-red rad.
  !
  !*      0.2    declarations of local variables
  !
  real :: ZIMPL = 0.5        ! implicit coefficient
  real :: ZEXPL = 0.5        ! explicit coefficient
  !
  real, dimension(3) :: ZA ! lower diag.
  real, dimension(3) :: ZB ! main  diag.
  real, dimension(3) :: ZC ! upper diag.
  real, dimension(3) :: ZY ! r.h.s.
  real, dimension(3) :: ZX ! solution
  !
  real :: ZRHO_AC_ROOF ! conductance * rho
  real :: ZRHO_AC_ROOF_WAT ! ! conductance * rho (for water)
  real :: ZDQSAT_ROOF  ! dq_sat/dTs
  real, dimension(size(PT_ROOF)) :: ZMTC_O_D
  ! mean thermal conductivity over distance between 2 layers
  real, dimension(size(PT_ROOF)) :: ZHC_D_ROOF
  ! thermal capacity times layer depth
  !
  integer :: IROOF_LAYER           ! number of roof layers
  integer :: JLAYER                ! loop counter
  !-----------------------------------------------------------------------------
  !
  ZABS_LW_ROOF = 0.
  !
  !*      1.     Layer thermal properties
  !              ------------------------
  !
  IROOF_LAYER = size(PT_ROOF)
  ZMTC_O_D(:) = 0.
  ZHC_D_ROOF(:) = 0.
  !
  do JLAYER=1,IROOF_LAYER-1
     ZMTC_O_D(JLAYER) = 2./(   PD_ROOF(JLAYER  )/PTC_ROOF(JLAYER  ) &
          + PD_ROOF(JLAYER+1)/PTC_ROOF(JLAYER+1) )
     ZHC_D_ROOF(JLAYER) = PHC_ROOF(JLAYER) * PD_ROOF(JLAYER)
  end do
  !
  ZMTC_O_D(IROOF_LAYER) = 2. * PTC_ROOF(IROOF_LAYER) &
       / PD_ROOF (IROOF_LAYER)
  ZMTC_O_D(IROOF_LAYER) = 1./(  1./ZMTC_O_D(IROOF_LAYER)      &
       + 1./(XCPD*PRHOA*ZAC_BLD)   )
  !
  ZHC_D_ROOF(IROOF_LAYER) = PHC_ROOF(IROOF_LAYER) &
       * PD_ROOF (IROOF_LAYER)
  !
  !-----------------------------------------------------------------------------
  !
  !*      2.     Roof Ts coefficients
  !              --------------------
  !
  !
  ZRHO_AC_ROOF     = PRHOA * ZAC_ROOF    
  ZRHO_AC_ROOF_WAT = PRHOA * ZAC_ROOF_WAT
  !
  !
  !*      2.1    dqsat/dTs, and humidity for roofs
  !              ---------------------------------
  ! Modification here. ZDQSAT_ROOF will be calculated locally
  !ZDQSAT_ROOF = DQSAT(PTS_ROOF,PPS,ZQSAT_ROOF)
  !
  !*       2.1.1.    COMPUTE SATURATION VAPOR PRESSURE
  !              ---------------------------------
  !
  !
  !*       2.1.2.    DERIVATION ACCORDING TO TEMPERATURE
  !              -----------------------------------
  !
  ZDQSAT_ROOF=(rslifp(PPS,PTS_ROOF))
  !
  !
  !
  !*      2.2    coefficients
  !              ------------
  !
  ZA(1) =   0.

  ZB(1) =   ZHC_D_ROOF(1)/PTSTEP                                      &
       + ZIMPL * ( (4. * PEMIS_ROOF * XSTEFAN * PTS_ROOF**3              &
       + ZRHO_AC_ROOF * XCPD                                  &
       + ZRHO_AC_ROOF_WAT * XLVTT * ZDELT_ROOF * ZDQSAT_ROOF  &
       )                                 &
       + ZMTC_O_D(1) )

  ZC(1) = ZIMPL * ( - ZMTC_O_D(1) )
  !
  ZY(1) =   ZHC_D_ROOF(1)/PTSTEP * PT_ROOF(1)                         &
       + (ZABS_SW_ROOF + PEMIS_ROOF*PLW_RAD                      &
       + ZRHO_AC_ROOF * XCPD * ZTA                               &
       - ZRHO_AC_ROOF_WAT * XLVTT * ZDELT_ROOF * ZQSAT_ROOF      &
       + ZRHO_AC_ROOF_WAT * XLVTT * ZDELT_ROOF * ZQA)            &
       + ZIMPL * ( 3. * PEMIS_ROOF * XSTEFAN * PTS_ROOF** 4              &
       + ZRHO_AC_ROOF_WAT * XLVTT * ZDELT_ROOF             &
       * ZDQSAT_ROOF * PTS_ROOF)         & 
       + ZEXPL * ((- PEMIS_ROOF * XSTEFAN * PTS_ROOF** 4                 &
       - ZRHO_AC_ROOF * XCPD * PTS_ROOF)                     &
       - ZMTC_O_D(1) * PT_ROOF(1)                            &
       + ZMTC_O_D(1) * PT_ROOF(2)  )
  !
  !-----------------------------------------------------------------------------
  !
  !*      3.     Other layers coefficients
  !              -------------------------
  !
  do JLAYER=2,IROOF_LAYER-1
     ZA(JLAYER) = ZIMPL * ( - ZMTC_O_D(JLAYER-1))

     ZB(JLAYER) =   ZHC_D_ROOF(JLAYER)/PTSTEP                         &
          + ZIMPL * (   ZMTC_O_D(JLAYER-1)                            &
          + ZMTC_O_D(JLAYER  ))

     ZC(JLAYER) = ZIMPL * ( - ZMTC_O_D(JLAYER  ))
     !
     ZY(JLAYER) =   ZHC_D_ROOF(JLAYER)/PTSTEP * PT_ROOF(JLAYER)       &
          + ZEXPL * (   ZMTC_O_D(JLAYER-1) * PT_ROOF(JLAYER-1)        &
          - ZMTC_O_D(JLAYER-1) * PT_ROOF(JLAYER  )        &
          - ZMTC_O_D(JLAYER  ) * PT_ROOF(JLAYER  )        &
          + ZMTC_O_D(JLAYER  ) * PT_ROOF(JLAYER+1) )
  end do
  !
  !-----------------------------------------------------------------------------
  !
  !*      4.     Inside layer coefficients
  !              -------------------------
  !
  ZA(IROOF_LAYER) =                                                      &
       ZIMPL * ( - ZMTC_O_D(IROOF_LAYER-1))

  ZB(IROOF_LAYER) =   ZHC_D_ROOF(IROOF_LAYER) / PTSTEP                   &
       + ZIMPL * (   ZMTC_O_D(IROOF_LAYER-1)                        &
       + ZMTC_O_D(IROOF_LAYER  ))

  ZC(IROOF_LAYER) =   0.
  !
  ZY(IROOF_LAYER) =     ZHC_D_ROOF(IROOF_LAYER)/PTSTEP                   &
       * PT_ROOF(IROOF_LAYER)   &
       + ZMTC_O_D(IROOF_LAYER) * PTI_BLD                &
       + ZEXPL * (   ZMTC_O_D(IROOF_LAYER-1)                        &
       * PT_ROOF(IROOF_LAYER-1)                        &
       - ZMTC_O_D(IROOF_LAYER-1)                        &
       * PT_ROOF(IROOF_LAYER  )                        &
       - ZMTC_O_D(IROOF_LAYER  )                        &
       * PT_ROOF(IROOF_LAYER  ) )
  !
  !
  !-----------------------------------------------------------------------------
  !
  !*      5.     Tri-diagonal system resolution
  !              ------------------------------
  !ok
  call TRID(ZX,ZA,ZB,ZC,ZY,3)
  !ZA, lower diag.
  !ZB, main  diag.
  !ZC, upper diag.
  !ZY, r.h.s.
  !ZX  solution
  PT_ROOF(:) = ZX(:)
  !
  !-----------------------------------------------------------------------------
  !
  !*      6.     Roof surface temperature
  !              ------------------------
  !
  PTS_ROOF = PT_ROOF(1)
  !
  !-----------------------------------------------------------------------------
  !
  !*      7.     Infra-red radiation absorbed by roofs
  !              -------------------------------------
  !
  ZABS_LW_ROOF = PEMIS_ROOF * (PLW_RAD - XSTEFAN * PTS_ROOF** 4)
  !
  !-----------------------------------------------------------------------------
  !
  !*      8.     New saturated specified humidity near the roof surface
  !              ------------------------------------------------------
  !
  ZQSAT_ROOF=rslif(PPS,PTS_ROOF)

  !
  !-----------------------------------------------------------------------------
  !
  return
end subroutine ROOF_LAYER_E_BUDGET
!   ##########################################################################
!   ##########################################################################
subroutine ROAD_WALL_LAYER_E_BUDGET(PTS_ROAD, PT_ROAD,                     &
     PTS_WALL, PT_WALL, PQSAT_ROAD,           &
     PT_CANYON, PQ_CANYON,                    &
     PTA, PQA, PPS,                           &
     PLW_RAD,  PTSTEP,                        &
     PH_TRAFFIC, PLE_TRAFFIC,                 &
     PBLD, PWALL_O_ROAD,                      &
     PEMIS_ROAD, PSVF_ROAD,                   &
     PHC_ROAD,PTC_ROAD,PD_ROAD,               &
     PEMIS_WALL, PSVF_WALL,                   &
     PHC_WALL,PTC_WALL,PD_WALL,               &
     PTI_BLD, PAC_BLD, PTI_ROAD,              &
     PDELT_ROAD,                              &
     PRHOA, PAC_WALL,                         &
     PAC_ROAD, PAC_ROAD_WAT, PAC_TOP,         &
     PABS_SW_ROAD, PABS_SW_WALL,              &
     PABS_LW_ROAD, PABS_LW_WALL               )
  !   ##########################################################################
  !
  !!****  *ROAD_WALL_LAYER_E_BUDGET*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the evoultion of roads and walls surface temperatures
  !         
  !     
  !!**  METHOD
  !     ------
  !
  !    6 : equations for evolution of Ts_road and Ts_wall simultaneously
  !        *************************************************************
  !
  !     dTw_k(t) / dt = 1/(dw_k*Cw_k) * (- 2*Kw_k-1*(Tw_k-Tw_k-1)/(dw_k-1 +dw_k) 
  !                                      - 2*Kw_k  *(Tw_k-Tw_k+1)/(dw_k+1 +dw_k) )
  !
  !     dTw_1(t) / dt = 1/(dw_1*Cw_1) * (  Rn_w - H_w - LE_w 
  !                                      - 2*Kw_1*(Tw_1-Tw_2)/(dw_1 +dw_2)       )
  !
  !     dTr_1(t) / dt = 1/(dr_1*Cr_1) * (  Rn_r - H_r - LE_r 
  !                                      - 2*Kr_1*(Tr_1-Tr_2)/(dr_1 +dr_2)       )
  !
  !     dTr_k(t) / dt = 1/(dr_k*Cr_k) * (- 2*Kr_k-1*(Tr_k-Tr_k-1)/(dr_k-1 +dr_k) 
  !                                      - 2*Kr_k  *(Tr_k-Tr_k+1)/(dr_k+1 +dr_k) )
  !
  !       with
  !
  !   K*_k  = (d*_k+ d*_k+1)/(d*_k/k*_k+ d*_k+1/k*_k+1)
  !
  !   Rn_w = abs_Rg_w 
  !  - sigma * emis_w                                                   * Ts_w**4 (t+dt)
  !  +         emis_w                       *      SVF_w                * LWR
  !  + sigma * emis_w * emis_r              *      SVF_w                * Ts_r**4 (t+dt)
  !  + sigma * emis_w * emis_w              * (1-2*SVF_w)               * Ts_w**4 (t+dt)
  !  +         emis_w            (1-emis_r) *      SVF_r  *      SVF_w  * LWR
  !  +         emis_w            (1-emis_w) *      SVF_w  * (1-2*SVF_w) * LWR
  !  + sigma * emis_w * emis_w * (1-emis_r) *      SVF_w  * (1-  SVF_r) * Ts_w**4 (t+dt)
  !  + sigma * emis_w * emis_w * (1-emis_w) * (1-2*SVF_w) * (1-2*SVF_w) * Ts_w**4 (t+dt)
  !  + sigma * emis_w * emis_r * (1-emis_w) *      SVF_w  * (1-2*SVF_w) * Ts_r**4 (t+dt)
  !
  !   Rn_r = abs_Rg_r
  !  - sigma * emis_r                                                   * Ts_r**4 (t+dt)
  !  +         emis_r                       *    SVF_r                  * LWR
  !  + sigma * emis_r * emis_w              * (1-SVF_r)                 * Ts_w**4 (t+dt)
  !  +         emis_r            (1-emis_w) * (1-SVF_r)   *      SVF_w  * LWR
  !  + sigma * emis_r * emis_w * (1-emis_w) * (1-SVF_r)   * (1-2*SVF_w) * Ts_w**4 (t+dt)
  !  + sigma * emis_r * emis_r * (1-emis_w) * (1-SVF_r)   *      SVF_w  * Ts_r**4 (t+dt)
  !
  !  H_w  = rho Cp CH V ( Ts_w (t+dt) - Ta_canyon )
  !
  !  LE_w = rho Lv CH V ( qs_w (t+dt) - qa_canyon )
  !
  !  H_r  = rho Cp CH V ( Ts_r (t+dt) - Ta_canyon )
  !
  !  LE_r = rho Lv CH V ( qs_r (t+dt) - qa_canyon )
  !
  ! with again
  !                AC_can * Swall/Sroad * Twall + AC_can * Troad + AC_top * Ta + H_traffic/rho/Cp/Sroad
  !   Ta_canyon = -------------------------------------------------------------------------------------
  !                AC_can * Swall/Sroad         + AC_can         + AC_top
  !
  !
  !                 AC_can * delt_road * qsat(Troad) + AC_top * qa + LE_traffic/rho/Lv/Sroad
  !   qa_canyon = --------------------------------------------------------------------------
  !                 AC_can * delt_road               + AC_top
  !
  !
  ! where H_traffic and LE_traffic are scaled to road area.
  !
  !
  ! The system is implicited (or semi-implicited).
  !
  ! ZIMPL=1    ---> implicit system
  ! ZIMPL=0.5  ---> semi-implicit system
  ! ZIMPL=0    ---> explicit system
  !
  !!      
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original    23/01/98 
  !-----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  !
  use teb_vars_const
  use therm_lib, only: rslif,rslifp
  implicit none
  !
  !*      0.1    declarations of arguments
  !
  !INPUT/OUTPUT VARIABLES
  !
  real               :: PTS_ROAD     ! road surface temperature
  real, dimension(3) :: PT_ROAD      ! road layers temperatures
  real               :: PTS_WALL     ! wall surface temperature
  real, dimension(3) :: PT_WALL      ! wall layers temperatures
  real               :: PQSAT_ROAD   ! q_sat(Ts)
  real               :: PTI_ROAD     ! road deep temperature
  !
  !ONLY INPUT VARIABLES
  !
  real               :: PTA          ! atmospheric air temperature
  real               :: PQA          ! and specific humidity at roof level
  real               :: PPS          ! pressure at the surface
  real               :: PLW_RAD      ! atmospheric infrared radiation
  real               :: PTSTEP       ! time step
  real               :: PH_TRAFFIC   ! anthropogenic sensible
  !                                  ! heat fluxes due to traffic
  real               :: PLE_TRAFFIC  ! anthropogenic latent
  !                                  ! heat fluxes due to traffic
  real               :: PBLD         ! fraction of buildings
  real               :: PWALL_O_ROAD ! wall Surf. / road Surf.
  real               :: PEMIS_ROAD   ! road emissivity
  real, dimension(3) :: PHC_ROAD     ! heat capacity for road layers
  real, dimension(3) :: PTC_ROAD     ! thermal conductivity for road layers
  real, dimension(3) :: PD_ROAD      ! depth of road layers
  real               :: PSVF_ROAD    ! road sky view factor
  real               :: PEMIS_WALL   ! road emissivity
  real, dimension(3) :: PHC_WALL     ! heat capacity for wall layers
  real, dimension(3) :: PTC_WALL     ! thermal conductivity for wall layers
  real, dimension(3) :: PD_WALL      ! depth of wall layers
  real               :: PSVF_WALL    ! road sky view factor
  real               :: PTI_BLD      ! inside building temperatur
  real               :: PAC_BLD      ! aerodynamical conductance
  ! inside the building itself
  real               :: PDELT_ROAD   ! fraction of water
  real               :: PRHOA        ! rho
  real               :: PAC_WALL     ! aerodynamical conductance
  !                                  ! between wall and canyon
  real               :: PAC_ROAD     ! aerodynamical conductance
  !                                  ! between road and canyon
  real               :: PAC_ROAD_WAT ! aerodynamical conductance
  !                                  ! between road and canyon
  !                                  ! (for water)
  real               :: PAC_TOP      ! aerodynamical conductance
  !                                  ! between atmosphere and
  !                                  ! canyon top
  real               :: PABS_SW_ROAD ! absorbed solar radiation
  real               :: PABS_SW_WALL ! absorbed solar radiation
  !
  !Only OUTPUT VARIABLES
  !
  real               :: PABS_LW_ROAD ! absorbed infrared rad.
  real               :: PABS_LW_WALL ! absorbed infrared rad.
  real               :: PT_CANYON    ! air canyon temperature
  real               :: PQ_CANYON    ! and specific humidity
  real               :: PQ_PCANYON    ! previousspecific humidity 
  !
  !*      0.2    declarations of local variables
  !
  !
  real :: ZIMPL=0.5      ! implicit coefficient
  real :: ZEXPL=0.5      ! explicit coefficient
  !
  real, dimension(size(PT_ROAD)+size(PT_WALL)) ::  ZA,& ! lower diag.
       ZB,& ! main  diag.
       ZC,& ! upper diag.
       ZY,& ! r.h.s.
       ZX   ! solution

  !
  real :: ZLW_W_TO_W  ! L.W. interactions
  real :: ZLW_R_TO_W  ! from first Temp.
  real :: ZLW_W_TO_R  ! on second Temp.
  real :: ZLW_R_TO_R  !
  real :: ZLW_S_TO_W  ! idem. but from
  real :: ZLW_S_TO_R  ! sky rad.


  real :: ZDQSAT_ROAD ! dq_sat/dTs
  real :: ZRHO_AC_W   ! rho * conductance (for walls)
  real :: ZRHO_AC_R   ! rho * conductance (for roads)
  real :: ZRHO_ACF_R  ! rho * conductance
  real :: ZRHO_ACF_R_WAT  ! rho * conductance for water
  real :: ZSAC_T      ! weighted sum 
  !                   ! of conductances
  !                   ! for PT_CANYON
  real :: ZSAC_Q      ! weighted sum 
  !                   ! of conductances
  !                   ! for PQ_CANYON
  real, dimension(size(PT_ROAD)) :: ZMTC_O_D_ROAD ! mean thermal 
  !                                               ! conductivity over distance 
  !                                               ! between 2 layers
  real, dimension(size(PT_WALL)) :: ZMTC_O_D_WALL ! mean thermal 
  !                                               ! conductivity over distance
  !                                               ! between 2 layers
  real, dimension(size(PT_ROAD)) :: ZHC_D_ROAD    ! thermal capacity
  !                                               ! times layer depth
  real, dimension(size(PT_WALL)) :: ZHC_D_WALL    ! thermal capacity 
  !                                               ! times layer depth
  !
  integer :: IROAD_LAYER           ! number of road layers
  integer :: IWALL_LAYER           ! number of wall layers
  integer :: ILAYER                ! current layer
  integer :: JLAYER                ! loop counter
  !-----------------------------------------------------------------------------
  !
  !
  !
  PABS_LW_ROAD = 0.
  PABS_LW_WALL = 0.
  !
  !*      1.     Layer thermal properties
  !              ------------------------
  !
  !*      1.1    Roads
  !              -----
  !
  IROAD_LAYER = size(PT_ROAD)
  ZMTC_O_D_ROAD(:) = 0.
  !
  do JLAYER=1,IROAD_LAYER-1
     ZMTC_O_D_ROAD(JLAYER) = 2./(  PD_ROAD(JLAYER  )/PTC_ROAD(JLAYER  ) &
          + PD_ROAD(JLAYER+1)/PTC_ROAD(JLAYER+1) )
     ZHC_D_ROAD   (JLAYER) = PHC_ROAD(JLAYER) * PD_ROAD (JLAYER)
  end do
  !
  ZMTC_O_D_ROAD(IROAD_LAYER) = 2. * PTC_ROAD(IROAD_LAYER) &
       / PD_ROAD (IROAD_LAYER)
  ZHC_D_ROAD   (IROAD_LAYER) = PHC_ROAD(IROAD_LAYER) &
       * PD_ROAD (IROAD_LAYER)
  !
  !*      1.2    Walls
  !              -----
  !
  IWALL_LAYER = size(PT_WALL)
  ZMTC_O_D_WALL(:) = 0.
  !
  do JLAYER=1,IWALL_LAYER-1
     ZMTC_O_D_WALL(JLAYER) = 2./(  PD_WALL(JLAYER  )/PTC_WALL(JLAYER  ) &
          + PD_WALL(JLAYER+1)/PTC_WALL(JLAYER+1) )
     ZHC_D_WALL   (JLAYER) = PHC_WALL(JLAYER) * PD_WALL (JLAYER)
  end do
  !
  ZMTC_O_D_WALL(IWALL_LAYER) = 2. * PTC_WALL(IWALL_LAYER) &
       / PD_WALL (IWALL_LAYER)
  ZMTC_O_D_WALL(IWALL_LAYER) = 1./(  1./ZMTC_O_D_WALL(IWALL_LAYER)    &
       + 1./(XCPD*PRHOA*PAC_BLD)      )
  !
  ZHC_D_WALL   (IWALL_LAYER) = PHC_WALL(IWALL_LAYER) &
       * PD_WALL (IWALL_LAYER)
  !
  !-----------------------------------------------------------------------------
  !
  !*      2.    Preliminaries
  !             -------------
  !
  !*      2.1    snow-free surface fraction
  !              --------------------------
  !
  !
  !*      2.3    flux properties
  !              ---------------
  !
  ZSAC_T  =                  PAC_ROAD &
       + PWALL_O_ROAD * PAC_WALL &
       +                PAC_TOP 
  ZSAC_Q  = PDELT_ROAD * PAC_ROAD_WAT + PAC_TOP


  ZRHO_AC_R  = PRHOA * PAC_ROAD
  ZRHO_AC_W  = PRHOA * PAC_WALL
  ZRHO_ACF_R = PRHOA * PAC_ROAD 
  ZRHO_ACF_R_WAT = PRHOA * PAC_ROAD_WAT
  !
  !*      2.4    qsat, dqsat/dTs, and humidity for roads
  !              ---------------------------------------
  !

  ZDQSAT_ROAD=(rslifp(PPS,PTS_ROAD))

  !-----------------------------------------------------------------------------
  !
  !*      3.     LW properties
  !              -------------
  !ok
  call URBAN_LW_COEF(PEMIS_ROAD, PSVF_ROAD, PEMIS_WALL, PSVF_WALL,   &
       ZLW_W_TO_W, ZLW_R_TO_W, ZLW_W_TO_R, ZLW_R_TO_R, &
       ZLW_S_TO_W, ZLW_S_TO_R)

  !
  !-----------------------------------------------------------------------------
  !
  !*      4.    Inside wall layer coefficients
  !             ------------------------------
  !
  ILAYER=1
  !
  ZA(ILAYER) =   0.

  ZB(ILAYER) =   ZHC_D_WALL(IWALL_LAYER) / PTSTEP                          &
       + ZIMPL * (   ZMTC_O_D_WALL(IWALL_LAYER  )                          &
       + ZMTC_O_D_WALL(IWALL_LAYER-1)                          &
       )

  ZC(ILAYER) =                                                             &
       + ZIMPL * ( - ZMTC_O_D_WALL(IWALL_LAYER-1)                          &
       )
  !
  ZY(ILAYER) =   ZHC_D_WALL(IWALL_LAYER) / PTSTEP                          &
       * PT_WALL(IWALL_LAYER)      &
       + ZMTC_O_D_WALL(IWALL_LAYER) * PTI_BLD                  &
       + ZEXPL * ( - ZMTC_O_D_WALL(IWALL_LAYER  )                          &
       * PT_WALL(IWALL_LAYER  )                          &
       - ZMTC_O_D_WALL(IWALL_LAYER-1)                          &
       * PT_WALL(IWALL_LAYER  )                          &
       + ZMTC_O_D_WALL(IWALL_LAYER-1)                          &
       * PT_WALL(IWALL_LAYER-1)                          &
       )
  !
  !-----------------------------------------------------------------------------
  !
  !*      5.     Other wall layers coefficients
  !              ------------------------------
  !
  do JLAYER=2,IWALL_LAYER-1

     ILAYER=IWALL_LAYER-JLAYER+1

     ZA(ILAYER) =                                                       &
          ZIMPL * ( - ZMTC_O_D_WALL(JLAYER  )                           &
          )

     ZB(ILAYER) =   ZHC_D_WALL(JLAYER)/PTSTEP                           &
          + ZIMPL * (   ZMTC_O_D_WALL(JLAYER  )                           &
          + ZMTC_O_D_WALL(JLAYER-1)                           &
          )

     ZC(ILAYER) =                                                       &
          ZIMPL * ( - ZMTC_O_D_WALL(JLAYER-1)                           &
          )
     !
     ZY(ILAYER) =   ZHC_D_WALL(JLAYER)/PTSTEP * PT_WALL(JLAYER)         &
          + ZEXPL * (    ZMTC_O_D_WALL(JLAYER  ) * PT_WALL(JLAYER+1)   &
          - ZMTC_O_D_WALL(JLAYER  ) * PT_WALL(JLAYER  )   &
          - ZMTC_O_D_WALL(JLAYER-1) * PT_WALL(JLAYER  )   &
          + ZMTC_O_D_WALL(JLAYER-1) * PT_WALL(JLAYER-1)   &
          )
  end do
  !
  !-----------------------------------------------------------------------------
  !
  !*      6.     Surface wall layer coefficients
  !              -------------------------------
  !
  ILAYER=IWALL_LAYER
  !
  !
  ZA(ILAYER) =                                                           &
       ZIMPL * ( - ZMTC_O_D_WALL(1)                                      &
       )

  ZB(ILAYER) =   ZHC_D_WALL(1)/PTSTEP                                    &
       + ZIMPL * ( - 4. *PTS_WALL**3 * ZLW_W_TO_W                    &
       + ZRHO_AC_W * XCPD                                      &
       * (1.-PAC_WALL*PWALL_O_ROAD/ZSAC_T )          &
       + ZMTC_O_D_WALL(1)                                      &
       )

  ZC(ILAYER) =                                                           &
       ZIMPL * ( - 4.*PTS_ROAD**3 * ZLW_R_TO_W                     &
       - ZRHO_AC_W * XCPD * PAC_ROAD/ZSAC_T                &
       )


  ZY(ILAYER) =   ZHC_D_WALL(1)/PTSTEP*PTS_WALL                           &
       + PABS_SW_WALL                                        &
       + PLW_RAD        * ZLW_S_TO_W                         &
       + PTS_WALL   **4 * ZLW_W_TO_W                         &
       + PTS_ROAD   **4 * ZLW_R_TO_W                         &
       + ZRHO_AC_W * XCPD * PTA                              &
       * PAC_TOP / ZSAC_T              &
       + PAC_WALL / ZSAC_T                                   &
       * (   PH_TRAFFIC   / (1.-PBLD ))      &
       + ZIMPL * (                                                       &
       - 4.*PTS_WALL**4 * ZLW_W_TO_W                         &
       - 4.*PTS_ROAD**4 * ZLW_R_TO_W                         &
       )                                                       &
       + ZEXPL * (   ZMTC_O_D_WALL(1) * PT_WALL(2)                       &
       - ZMTC_O_D_WALL(1) * PT_WALL(1)                       &
       - ZRHO_AC_W * XCPD * (1.-PAC_WALL*PWALL_O_ROAD        &
       /ZSAC_T     ) &
       * PTS_WALL                             &
       + ZRHO_AC_W  * XCPD * PAC_ROAD /ZSAC_T           &
       * PTS_ROAD                             &
       )
  !
  !-----------------------------------------------------------------------------
  !
  !*      7.     Surface road layer coefficients
  !              -------------------------------
  !
  ILAYER=IWALL_LAYER+1
  !
  !
  ZA(ILAYER) =                                                           &
       ZIMPL * ( - 4.*PTS_WALL**3 * ZLW_W_TO_R                     &
       - ZRHO_ACF_R  * XCPD * PAC_WALL                         &
       * PWALL_O_ROAD  / ZSAC_T                        &
       )

  ZB(ILAYER) =   ZHC_D_ROAD(1)/PTSTEP                                    &
       + ZIMPL * ( - 4.*PTS_ROAD**3 * ZLW_R_TO_R                    &
       + ZRHO_ACF_R * XCPD * (1.-PAC_ROAD/ZSAC_T)          &
       + ZRHO_ACF_R_WAT * XLVTT * PDELT_ROAD * ZDQSAT_ROAD         &
       * (1.-PAC_ROAD_WAT *PDELT_ROAD /ZSAC_Q )        &
       + ZMTC_O_D_ROAD(1)                                      &
       )

  ZC(ILAYER) =                                                           &
       ZIMPL * ( - ZMTC_O_D_ROAD(1)                                      &
       )

  ZY(ILAYER) =   ZHC_D_ROAD(1)/PTSTEP*PT_ROAD(1)                         &
       + PABS_SW_ROAD 				       &
       + PLW_RAD	  * ZLW_S_TO_R  		       &
       + PTS_WALL   **4 * ZLW_W_TO_R  		       &
       + PTS_ROAD   **4 * ZLW_R_TO_R  		       &
       + ZRHO_ACF_R * XCPD * PTA                             &
       * PAC_TOP / ZSAC_T             &
       + PAC_ROAD  / ZSAC_T                                  &
       * (   PH_TRAFFIC  / (1.-PBLD))             &
       - ZRHO_ACF_R_WAT  * XLVTT * PDELT_ROAD                    &
       * (  PQSAT_ROAD                               &
       -(   PQSAT_ROAD      *PDELT_ROAD          &
       *PAC_ROAD_WAT        &
       + PQA             *PAC_TOP             &
       ) / ZSAC_Q                              &
       )                                           &
       + PAC_ROAD_WAT * PDELT_ROAD/ ZSAC_Q                   &
       * (   PLE_TRAFFIC  / (1.-PBLD))               &
       + ZIMPL * ( - 4.*PTS_WALL**4 * ZLW_W_TO_R                         &
       - 4.*PTS_ROAD**4 * ZLW_R_TO_R                         &
       + ZRHO_ACF_R_WAT * XLVTT * PDELT_ROAD                    &
       * (1.-PDELT_ROAD *PAC_ROAD_WAT          &
       /ZSAC_Q          ) &
       * ZDQSAT_ROAD  * PTS_ROAD               &
       )                                                       &
       + ZEXPL * (   ZRHO_ACF_R * XCPD * PAC_WALL * PWALL_O_ROAD         &
       * PTS_WALL / ZSAC_T                      &
       - ZRHO_ACF_R * XCPD * PTS_ROAD                        &
       * ( 1. - PAC_ROAD  /ZSAC_T  )         &
       - ZMTC_O_D_ROAD(1) * PT_ROAD(1)                       &
       + ZMTC_O_D_ROAD(1) * PT_ROAD(2)                       &
       )
  !
  !-----------------------------------------------------------------------------
  !
  !*      8.     Other road layers coefficients
  !              ------------------------------
  !
  do JLAYER=2,IROAD_LAYER-1

     ILAYER=IWALL_LAYER+JLAYER

     ZA(ILAYER) =                                                     &
          ZIMPL * ( - ZMTC_O_D_ROAD(JLAYER-1)                         &
          )

     ZB(ILAYER) =   ZHC_D_ROAD(JLAYER)/PTSTEP                         &
          + ZIMPL * (   ZMTC_O_D_ROAD(JLAYER-1)                         &
          + ZMTC_O_D_ROAD(JLAYER  )                         &
          )

     ZC(ILAYER) =                                                     &
          ZIMPL * ( - ZMTC_O_D_ROAD(JLAYER  )                         &
          )
     !
     ZY(ILAYER) =   ZHC_D_ROAD(JLAYER)/PTSTEP*PT_ROAD(JLAYER)         &
          + ZEXPL * (   ZMTC_O_D_ROAD(JLAYER-1) * PT_ROAD(JLAYER-1)   &
          - ZMTC_O_D_ROAD(JLAYER-1) * PT_ROAD(JLAYER  )   &
          - ZMTC_O_D_ROAD(JLAYER  ) * PT_ROAD(JLAYER  )   &
          + ZMTC_O_D_ROAD(JLAYER  ) * PT_ROAD(JLAYER+1)   &
          )
  end do
  !
  !-----------------------------------------------------------------------------
  !
  !*      9.     Inside road layer coefficients
  !              ------------------------------
  !
  ILAYER=IWALL_LAYER+IROAD_LAYER
  !
  ZA(ILAYER) =                                                      &
       ZIMPL * ( - ZMTC_O_D_ROAD(IROAD_LAYER-1)                     &
       )

  ZB(ILAYER) =   ZHC_D_ROAD(IROAD_LAYER) / PTSTEP                   &
       + ZIMPL * (   ZMTC_O_D_ROAD(IROAD_LAYER-1)                     &
       )

  ZC(ILAYER) =   0.
  !
  ZY(ILAYER) =   ZHC_D_ROAD(IROAD_LAYER) / PTSTEP                   &
       * PT_ROAD(IROAD_LAYER)   &
       + ZEXPL * (   ZMTC_O_D_ROAD(IROAD_LAYER-1)                     &
       * PT_ROAD(IROAD_LAYER-1)                     &
       - ZMTC_O_D_ROAD(IROAD_LAYER-1)                     &
       * PT_ROAD(IROAD_LAYER  )                     &
       )
  !
  !-----------------------------------------------------------------------------
  !
  !*     10.     Tri-diagonal system resolution
  !              ------------------------------
  !

  call TRID(ZX,ZA,ZB,ZC,ZY,6)
  !
  do JLAYER=1,IWALL_LAYER
     ILAYER=IWALL_LAYER-JLAYER+1
     PT_WALL(JLAYER) = ZX(ILAYER)
  end do
  !
  do JLAYER=1,IROAD_LAYER
     ILAYER=IWALL_LAYER+JLAYER
     PT_ROAD(JLAYER) = ZX(ILAYER)
  end do
  !
  !-----------------------------------------------------------------------------
  !
  !*     11.     special temperatures
  !              --------------------
  !
  !*     11.1    Road and wall surface temperatures
  !              ----------------------------------
  !
  PTS_ROAD = PT_ROAD(1)
  PTS_WALL = PT_WALL(1)
  !
  !*     11.2    Road deep temperature
  !              ---------------------
  !
  PTI_ROAD = PT_ROAD(IROAD_LAYER)
  !
  !-----------------------------------------------------------------------------
  !
  !*     12.    Road and wall absorbed infra-red radiation on snow-free surfaces
  !             ----------------------------------------------------------------
  !
  PABS_LW_ROAD =  ZLW_S_TO_R*PLW_RAD                 &
       + ZLW_R_TO_R*PTS_ROAD   **4          &
       + ZLW_W_TO_R*PTS_WALL   **4
  !
  PABS_LW_WALL =  ZLW_S_TO_W*PLW_RAD                 &
       + ZLW_W_TO_W*PTS_WALL   **4          &
       + ZLW_R_TO_W*PTS_ROAD   **4
  !
  !-----------------------------------------------------------------------------
  !
  !*     13.    Air canyon temperature at time t+dt
  !             -----------------------------------
  !
  PT_CANYON = (  PTS_ROAD    * PAC_ROAD                   &
       + PTS_WALL    * PAC_WALL * PWALL_O_ROAD    &
       + PTA         * PAC_TOP                    &
       + PH_TRAFFIC  / (1.-PBLD)/ PRHOA / XCPD ) &
       / ZSAC_T
  !-----------------------------------------------------------------------------
  !
  !*     14.     canyon air specific humidities
  !              ------------------------------
  !
  !*     14.1    New saturated specified humidity near the road surface
  !              ------------------------------------------------------
  !
  PQSAT_ROAD=rslif(PPS,PTS_ROAD)
  !
  !*     14.2    Canyon air specific humidity
  !              ----------------------------
  !
  PQ_PCANYON=PQ_CANYON
  PQ_CANYON = (  PQSAT_ROAD      * PAC_ROAD_WAT * PDELT_ROAD    & !linha 2676
       + PQA             * PAC_TOP                      )&
       !             + PLE_TRAFFIC  / (1.-PBLD) / PRHOA / XLVTT)     &
  / ZSAC_Q
  !
  !-----------------------------------------------------------------------------
  !
  return
end subroutine ROAD_WALL_LAYER_E_BUDGET


!   ##########################################################################
subroutine URBAN_LW_COEF(PEMIS_ROAD, PSVF_ROAD, PEMIS_WALL, PSVF_WALL,   &
     PLW_W_TO_W, PLW_R_TO_W, PLW_W_TO_R, PLW_R_TO_R, &
     PLW_S_TO_W, PLW_S_TO_R)
  !   ##########################################################################
  !
  !!****  *URBAN_LW_COEF*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the coefficients before each of the temperatures in the
  !     radiative budgets
  !         
  !     
  !!**  METHOD
  !     ------
  !
  ! without snow, the radiative budgets read:
  !
  !   Rn_w = abs_Rg_w 
  !  - sigma * emis_w                                                   * Ts_w**4 (t+dt)
  !  +         emis_w                       *      SVF_w                * LWR
  !  + sigma * emis_w * emis_r              *      SVF_w                * Ts_r**4 (t+dt)
  !  + sigma * emis_w * emis_w              * (1-2*SVF_w)               * Ts_w**4 (t+dt)
  !  +         emis_w            (1-emis_r) *      SVF_r  *      SVF_w  * LWR
  !  +         emis_w            (1-emis_w) *      SVF_w  * (1-2*SVF_w) * LWR
  !  + sigma * emis_w * emis_w * (1-emis_r) *      SVF_w  * (1-  SVF_r) * Ts_w**4 (t+dt)
  !  + sigma * emis_w * emis_w * (1-emis_w) * (1-2*SVF_w) * (1-2*SVF_w) * Ts_w**4 (t+dt)
  !  + sigma * emis_w * emis_r * (1-emis_w) *      SVF_w  * (1-2*SVF_w) * Ts_r**4 (t+dt)
  !
  !   Rn_r = abs_Rg_r
  !  - sigma * emis_r                                                   * Ts_r**4 (t+dt)
  !  +         emis_r                       *    SVF_r                  * LWR
  !  + sigma * emis_r * emis_w              * (1-SVF_r)                 * Ts_w**4 (t+dt)
  !  +         emis_r            (1-emis_w) * (1-SVF_r)   *      SVF_w  * LWR
  !  + sigma * emis_r * emis_w * (1-emis_w) * (1-SVF_r)   * (1-2*SVF_w) * Ts_w**4 (t+dt)
  !  + sigma * emis_r * emis_r * (1-emis_w) * (1-SVF_r)   *      SVF_w  * Ts_r**4 (t+dt)
  !
  !
  !
  !!      
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original    08/09/98 
  !!                  14/05/2002 Edmilson Freitas: changing from modules to 
  !!                             normal subroutines. Inclusion of the tebconst.h
  !!                             file, which has the constants used by TEB.!
  !!                             Elimination of declarations that are not needed.
  !-----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  use teb_vars_const
  implicit none
  !
  !*      0.1    declarations of arguments
  !
  !   INPUT VARIABLES
  !
  real :: PEMIS_ROAD  ! road emissivity
  real :: PSVF_ROAD   ! road sky view factor
  real :: PEMIS_WALL  ! wall emissivity
  real :: PSVF_WALL   ! wall sky view factor
  !
  !OUTPUT VARIABLES
  !
  real :: PLW_W_TO_W  ! L.W. interactions
  real :: PLW_R_TO_W  ! from first Temp.
  real :: PLW_W_TO_R  ! on second Temp.
  real :: PLW_R_TO_R  !
  real :: PLW_S_TO_W  ! idem. but from
  real :: PLW_S_TO_R  ! sky rad.
  !
  !*      0.2    declarations of local variables
  !
  !
  !
  PLW_W_TO_W =     XSTEFAN * PEMIS_WALL                      &
       * ( -   1.                                        &
       +       PEMIS_WALL  * (1.-2.*PSVF_WALL) &
       +       PEMIS_WALL  * (1.-   PSVF_ROAD) &
       * (1.-PEMIS_ROAD )*        PSVF_WALL  &
       +       PEMIS_WALL  * (1.-2.*PSVF_WALL) &
       * (1.-PEMIS_WALL) * (1.-2.*PSVF_WALL) &
       )
  !
  PLW_R_TO_W =   XSTEFAN * PEMIS_ROAD                        &
       * PEMIS_WALL                        &
       * PSVF_WALL                         &
       * ( 1. + (1.-PEMIS_WALL) * (1.-2.*PSVF_WALL))
  !
  PLW_S_TO_W =  PEMIS_WALL * PSVF_WALL                       &
       * (    1.                                      &
       + (1.-PEMIS_ROAD ) *        PSVF_ROAD       &
       + (1.-PEMIS_WALL) * (1.-2.*PSVF_WALL)      &
       )
  !
  !
  PLW_R_TO_R =            XSTEFAN * PEMIS_ROAD               &
       * (  -   1.                                    &
       +       PEMIS_ROAD * (1.-PEMIS_WALL)      &
       * (1.-PSVF_ROAD) * PSVF_WALL   &
       )
  !
  PLW_W_TO_R =   XSTEFAN * PEMIS_ROAD                        &
       * PEMIS_WALL                        &
       * (1. -  PSVF_ROAD)                 &
       * ( 1. + (1.-PEMIS_WALL) * (1.-2.*PSVF_WALL) )
  !
  PLW_S_TO_R =  PEMIS_ROAD* PSVF_ROAD                        &
       + PEMIS_ROAD*(1.-PEMIS_WALL)                   &
       *(1.-PSVF_ROAD)*PSVF_WALL
  !
  !
  !-----------------------------------------------------------------------------
  !
  return
end subroutine URBAN_LW_COEF
!   ##########################################################################
!   ##########################################################################
subroutine URBAN_HYDRO(PWS_ROOF_MAX,PWS_ROAD_MAX, PWS_ROOF, PWS_ROAD,  &
     PRR, PTSTEP, PBLD, PLE_ROOF, PLE_ROAD,          &
     PRUNOFF_ROOF,                                   &
     PRUNOFF_ROAD,                                   &
     PRUNOFF_TOWN                                    )
  !   ##########################################################################
  !
  !!****  *URBAN_HYDRO*  
  !!
  !!    PURPOSE
  !!    -------
  !
  !     Computes the evolution of prognostic water reservoirs
  !     of urbanized areas.
  !         
  !     
  !!**  METHOD
  !     ------
  !
  !
  !   The roof reservoir runoff goes directly into the road reservoir.
  !
  !   Runoff occurs for road reservoir (too much water), as well as drainage
  !   (evacuation system, typical time scale: 1 day)
  !
  !!      
  !!    AUTHOR
  !!    ------
  !!
  !!	V. Masson           * Meteo-France *
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!      Original    23/01/98 
  !-----------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  !
  use teb_vars_const
  implicit none
  !
  !*      0.1    declarations of arguments
  !
  ! INPUT/OUTPUT VARIABLES
  !
  real :: PWS_ROOF     ! roof water reservoir
  real :: PWS_ROAD     ! road water reservoir
  !
  !Only INPUT VARIABLES
  !
  real :: PWS_ROOF_MAX    ! maximum deepness of roof water reservoir
  real :: PWS_ROAD_MAX    ! maximum deepness of road water reservoir

  real :: PRR          ! rain rate
  real :: PTSTEP       ! time step
  real :: PBLD         ! fraction of buildings
  real :: PLE_ROOF     ! latent heat flux over roof
  real :: PLE_ROAD     ! latent heat flux over road
  !
  ! Only OUTPUT VARIABLES
  !
  real :: PRUNOFF_ROOF ! runoff (kg/m2/s)
  real :: PRUNOFF_ROAD ! runoff (kg/m2/s)
  real :: PRUNOFF_TOWN ! runoff (kg/m2/s)
  !
  !
  !-----------------------------------------------------------------------------
  !
  !*      1.     Roof reservoir evolution
  !              ------------------------
  !
  !
  !                                           evolution of the water reservoir
  !                                           (if we don't consider the runoff)
  !                                           PRR in kg/m2/s therefore PWS in mm
  !
  PWS_ROOF =  PWS_ROOF                                   &
       - PTSTEP * ( PLE_ROOF / XLVTT - PRR )
  !
  !                                           Ws_town must be positive
  !
  PWS_ROOF = max(0., PWS_ROOF)
  !
  !                                           if Ws_town > Ws_town_max,
  !                                           there is runoff
  !
  PWS_ROOF_MAX=1.

  PRUNOFF_ROOF = max(0., (PWS_ROOF - PWS_ROOF_MAX) /PTSTEP)  !divisao por ptstep foi retirada
  !
  PWS_ROOF = min(PWS_ROOF, PWS_ROOF_MAX)
  !
  !-----------------------------------------------------------------------------
  !
  !*      2.     Road reservoir evolution
  !              ------------------------
  !
  !
  !                                           evolution of the water reservoir
  !                                           (if we don't consider the runoff)
  !                                           PRR in kg/m2/s therefore PWS in mm
  !
  PWS_ROAD =  PWS_ROAD                                 &
       - PTSTEP * ( PLE_ROAD / XLVTT -  PRR )
  !
  !                                           Ws_town must be positive
  !
  PWS_ROAD = max(0., PWS_ROAD)
  !
  !                                           if Ws_town > Ws_town_max,
  !                                           there is runoff
  !
  !
  PWS_ROAD_MAX=1.
  PRUNOFF_ROAD = max(0., (PWS_ROAD - PWS_ROAD_MAX)/ PTSTEP  ) ! 
  !
  PWS_ROAD = min(PWS_ROAD, PWS_ROAD_MAX)
  !
  !-----------------------------------------------------------------------------
  !
  !*      3.     Area-averaged runoff
  !              --------------------
  !
  PRUNOFF_TOWN = PRUNOFF_ROAD * (1.-PBLD) + PRUNOFF_ROOF * PBLD
  !
  !-----------------------------------------------------------------------------
  !
  return
end subroutine URBAN_HYDRO

